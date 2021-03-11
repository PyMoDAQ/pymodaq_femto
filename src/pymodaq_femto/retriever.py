import sys
import subprocess

import logging
from pathlib import Path
import numpy as np
from PyQt5 import QtGui, QtWidgets, QtCore
from PyQt5.QtCore import Qt, QObject, pyqtSlot, QThread, pyqtSignal, QLocale
from PyQt5.QtGui import QIcon, QPixmap
from pyqtgraph.dockarea import Dock
from pyqtgraph.parametertree import Parameter, ParameterTree
from pymodaq.daq_utils import daq_utils as utils
from pymodaq.daq_utils import gui_utils as gutils
from pymodaq.daq_utils.h5modules import browse_data
from pymodaq.daq_utils.plotting.viewer2D.viewer2D_main import Viewer2D
from pymodaq.daq_utils.plotting.viewer1D.viewer1D_main import Viewer1D
from pymodaq.daq_utils.plotting.viewer0D.viewer0D_main import Viewer0D
from pymodaq_femto.simulation import Simulator, methods, MplCanvas, NavigationToolbar, MeshDataPlot, nlprocesses, \
    PulsePlot
from pymodaq_femto.result_graphics import RetrievalResultPlot
from collections import OrderedDict
from pypret import FourierTransform, Pulse, PNPS, lib, MeshData, random_gaussian
from pypret.frequencies import om2wl, wl2om, convert
import scipy
from scipy.fftpack import next_fast_len
from pymodaq.daq_utils.h5modules import H5BrowserUtil
from pyqtgraph.graphicsItems.GradientEditorItem import Gradients
from pymodaq_femto import _PNPS_CLASSES
from pypret.retrieval.retriever import _RETRIEVER_CLASSES

retriever_algos = list(_RETRIEVER_CLASSES.keys())

config = utils.load_config()
logger = utils.set_logger(utils.get_module_name(__file__))


class DataIn(OrderedDict):
    def __init__(self, name='', source='', trace_in=None, pulse_in=None, raw_spectrum=None, raw_trace=None, **kwargs):
        """class subclassing from OrderedDict defining data to be processed by the retriever, either experimental or
        simulated
        Parameters
        ----------
        name: (str) data identifier
        source: (str) either "simulated" or "experimental"
        trace_in: (MeshData) MeshData object as defined in pypret and containing the trace data
        pulse_in: (Pulse) Pulse object as defined in pypret and containing fundamental spectrum (at least)
        raw_spectrum: (dict) with data and axis as keys containing the spectrum and an Axis object
        raw_trace: (dict) with data axis and parameter_axis keys
        """
        if not isinstance(name, str):
            raise TypeError('name for the DataIn class should be a string')
        self['name'] = name
        if not isinstance(source, str):
            raise TypeError('source for the DataIn class should be a string')
        elif not ('simulated' in source or 'experimental' in source):
            raise ValueError('Invalid "source" for the DataIn class')
        self['source'] = source

        self['trace_in'] = trace_in
        self['pulse_in'] = pulse_in

        for k in kwargs:
            self[k] = kwargs[k]


def pulse_from_spectrum(wavelength, spectrum, pulse=None):
    """ Generates a pulse instance from a measured spectrum.
    """
    # scale to intensity over frequency, convert to amplitude and normalize
    spectrum = spectrum * wavelength * wavelength
    spectrum[spectrum < 0.0] = 0.0
    spectrum = np.sqrt(spectrum + 0.0j)
    spectrum /= spectrum.max()
    # calculate angular frequencies
    w = convert(wavelength, "wl", "om")
    if pulse is None:
        # create pulse parameters from the measured spectrum
        # choose center wavelength as the mean of the intensity
        w0 = lib.mean(w, lib.abs2(spectrum))
        # choose simulation grid that encompasses the measured spectrum
        dw = abs(np.mean(np.diff(w)))
        N = 4 * next_fast_len(int(abs(w[-1] - w[0]) / dw))
        ft = FourierTransform(N, dw=dw)
        pulse = Pulse(ft, w0, unit="om")
    # interpolate
    pulse.spectrum = scipy.interpolate.interp1d(w - pulse.w0, spectrum, bounds_error=False, fill_value=0.0)(pulse.w)
    return pulse


def preprocess(trace, signal_range=None, dark_signal_range=None):
    if dark_signal_range is not None:
        dark_signal = trace.copy()
        dark_signal.limit(dark_signal_range, axis=1)
        dark_signal = np.median(dark_signal.data, axis=1)
    if signal_range is not None:
        trace.limit(*signal_range)
    if dark_signal_range is not None:
        # subtract dark counts for every spectrum separately
        trace.data -= dark_signal[:, None]
    # normalize
    trace.normalize()
    return trace

# interpolate the measurement
def preprocess2(trace, pnps):
    if trace.units[1] == "m":
        # scaled in wavelength -> has to be corrected
        wavelength = trace.axes[1]
        frequency = convert(wavelength, "wl", "om")
        trace.scale(wavelength * wavelength)
        trace.normalize()
        trace.axes[1] = frequency
        trace.units[1] = "Hz"
    trace.interpolate(axis2=pnps.process_w)


def mask(x, y, where, **kwargs):
    y = scipy.interpolate.interp1d(x[~where], y[~where], **kwargs)(x)
    return y


class Retriever(QObject):
    """
    Main class initializing a DAQ_Scan module with its dashboard and scanning control panel
    """
    status_signal = pyqtSignal(str)

    params_in = [
        {'title': 'Data Info', 'name': 'data_in_info', 'type': 'group', 'children': [
            {'title': 'Trace type:', 'name': 'trace_type', 'type': 'list', 'values': methods,
             'tip': 'Characterization technique used to obtain this trace'},
            {'title': 'NL process:', 'name': 'nlprocess', 'type': 'list',
             'values': nlprocesses,
             'tip': 'Non Linear process used in the experiment'},
            {'title': 'Trace Info', 'name': 'trace_in_info', 'type': 'group', 'children': [
                {'title': 'Wl0 (nm)', 'name': 'wl0', 'type': 'float', 'value': 0, 'readonly': True,
                 'tip': 'Central spectrum wavelength in nanometers'},
                {'title': 'FWHM (nm)', 'name': 'wl_fwhm', 'type': 'float', 'value': 0, 'readonly': True,
                 'tip': 'FWHM of the spectrum in nanometers'},
                {'title': 'Param Size', 'name': 'trace_param_size', 'type': 'int', 'value': 0, 'readonly': True},
                {'title': 'Wavelentgh Size', 'name': 'trace_wl_size', 'type': 'int', 'value': 0, 'readonly': True},
                {'title': 'Scaling (m)', 'name': 'wl_scaling', 'type': 'float', 'value': 1e-9, 'readonly': False,
                 'tip': 'Scaling to go from the Trace wavelength values to wavelength in meters'},
                {'title': 'Scaling Parameter', 'name': 'param_scaling', 'type': 'float', 'value': 1e-15,
                 'readonly': False,
                 'tip': 'Scaling to go from the trace parameter values to delay in seconds, insertion in m (dscan) '
                        'or phase in rad (miips)'},

            ]},
            {'title': 'Spectrum Info', 'name': 'spectrum_in_info', 'type': 'group', 'children': [
                {'title': 'Wl0 (nm)', 'name': 'wl0', 'type': 'float', 'value': 0, 'readonly': True,
                 'tip': 'Central spectrum wavelength in nanometers'},
                {'title': 'FWHM (nm)', 'name': 'wl_fwhm', 'type': 'float', 'value': 0, 'readonly': True,
                 'tip': 'FWHM of the spectrum in nanometers'},
                {'title': 'Wavelength Size', 'name': 'spectrum_size', 'type': 'int', 'value': 0, 'readonly': True},
                {'title': 'Scaling (m)', 'name': 'wl_scaling', 'type': 'float', 'value': 1e-9, 'readonly': False,
                 'tip': 'Scaling to go from the spectrum wavelength values to wavelength in meters'},
            ]},
        ]},
        {'title': 'Processing', 'name': 'processing', 'type': 'group', 'children': [
            {'title': 'Grid settings:', 'name': 'grid_settings', 'type': 'group', 'children': [
                {'title': 'lambda0 (nm):', 'name': 'wl0', 'type': 'float', 'value': 750,
                 'tip': 'Central Wavelength of the Pulse spectrum and frequency grid'},
                {'title': 'Npoints:', 'name': 'npoints', 'type': 'list', 'values': [2 ** n for n in range(8, 16)],
                 'value': 1024,
                 'tip': 'Number of points for the temporal and Fourier Transform Grid'},
                {'title': 'Time resolution (fs):', 'name': 'time_resolution', 'type': 'float', 'value': 0.5,
                 'tip': 'Time spacing between 2 points in the time grid'},
            ]},
            {'title': 'Trace limits:', 'name': 'ROIselect', 'type': 'group', 'visible': False, 'children': [
                {'title': 'x0:', 'name': 'x0', 'type': 'int', 'value': 0, 'min': 0},
                {'title': 'y0:', 'name': 'y0', 'type': 'int', 'value': 0, 'min': 0},
                {'title': 'width:', 'name': 'width', 'type': 'int', 'value': 10, 'min': 1},
                {'title': 'height:', 'name': 'height', 'type': 'int', 'value': 10, 'min': 1},
            ]},
            {'title': 'Process trace', 'name': 'process_trace', 'type': 'action',
             'tip': 'Use one ROI to select a frequency area from the trace in order to remove the background'
                    ' and use ROISelect to reduce the area around the trace'},
            {'title': 'Process Spectrum', 'name': 'process_spectrum', 'type': 'action',
             'tip': 'Use ROIs to select frequency areas from the spectrum that should be removed'},
            {'title': 'Process Both', 'name': 'process_both', 'type': 'action',
             'tip': 'Process both the trace and the spectrum'},

        ]},
    ]

    params_retriever = [
        {'title': 'Algo type:', 'name': 'algo_type', 'type': 'list', 'values': retriever_algos,
         'tip': 'Retriever Algorithm'},
        {'title': 'Verbose Info:', 'name': 'verbose', 'type': 'bool', 'value': True,
         'tip': 'Display infos during retrieval'},
        {'title': 'Max iteration:', 'name': 'max_iter', 'type': 'int', 'value': 30,
         'tip': 'Max iteration for the algorithm'},
        {'title': 'Start Retrieval', 'name': 'start', 'type': 'action',
         'tip': 'Start the retrieval process'},
           ]

    def __init__(self, dockarea=None, dashboard=None):
        """

        Parameters
        ----------
        dockarea: (dockarea) instance of the modified pyqtgraph Dockarea (see daq_utils)
        dashboard: (DashBoard) instance of the pymodaq dashboard
        """
        QLocale.setDefault(QLocale(QLocale.English, QLocale.UnitedStates))
        logger.info('Initializing Retriever Extension')
        super().__init__()

        self.h5browse = H5BrowserUtil()

        self.dockarea = dockarea
        self.dashboard = dashboard
        self.mainwindow = self.dockarea.parent()

        self.settings_data_in = Parameter.create(name='dataIN_settings', type='group', children=self.params_in)
        self.settings_data_in.sigTreeStateChanged.connect(self.settings_in_changed)

        self.settings_retriever = Parameter.create(name='retriever_settings', type='group',
                                                   children=self.params_retriever)
        self.settings_retriever.sigTreeStateChanged.connect(self.settings_retriever_changed)

        self.setupUI()
        self.create_menu(self.mainwindow.menuBar())
        self.simulator = None
        self.data_in = None
        self.ft = None
        self.retriever = None
        self.pnps = None
        self.settings_data_in.child('processing', 'process_trace').sigActivated.connect(self.process_trace)
        self.settings_data_in.child('processing', 'process_spectrum').sigActivated.connect(self.process_spectrum)
        self.settings_data_in.child('processing', 'process_both').sigActivated.connect(self.process_both)
        self.settings_retriever.child('start').sigActivated.connect(self.start_retriever)

        self.viewer_trace_in.ROI_select_signal.connect(self.update_ROI)
        self.viewer_trace_in.ROIselect_action.triggered.connect(self.show_ROI)

    def create_menu(self, menubar):
        """
            Create the menubar object looking like :
        """
        menubar.clear()

        # %% create Settings menu
        self.main_menu = menubar.addMenu('Main')
        self.quit_action = self.main_menu.addAction('Quit')
        self.restart_action = self.main_menu.addAction('Restart')
        self.quit_action.triggered.connect(self.quit_fun)
        self.restart_action.triggered.connect(self.restart_fun)

        self.data_in_menu = menubar.addMenu('Data In')
        self.data_in_menu.addAction(self.load_trace_in_action)
        self.data_in_menu.addAction(self.gen_trace_in_action)

    def settings_in_changed(self, param, changes):
        for param, change, data in changes:
            path = self.settings_data_in.childPath(param)
            if change == 'childAdded':
                pass
            elif change == 'parent':
                pass
            elif change == 'value':
                if param.name() == 'method':
                    self.settings_data_in.child('algo',
                                                'nlprocess').setLimits(list(_PNPS_CLASSES[param.value()].keys()))

    def settings_retriever_changed(self, param, changes):
        for param, change, data in changes:
            path = self.settings_data_in.childPath(param)
            if change == 'childAdded':
                pass
            elif change == 'parent':
                pass
            elif change == 'value':
                if param.name() == 'method':
                    pass

    def quit_fun(self):
        """

        """
        try:
            if hasattr(self, 'mainwindow'):
                self.mainwindow.close()

        except Exception as e:
            logger.exception(str(e))

    def restart_fun(self, ask=False):
        ret = False
        mssg = QtWidgets.QMessageBox()
        if ask:
            mssg.setText('You have to restart the application to take the modifications into account!')
            mssg.setInformativeText("Do you want to restart?")
            mssg.setStandardButtons(mssg.Ok | mssg.Cancel)
            ret = mssg.exec()

        if ret == mssg.Ok or not ask:
            self.quit_fun()
            subprocess.call([sys.executable, __file__])

    def setupUI(self):
        self.ui = QObject()

        #  create main docks

        self.ui.dock_data_in = Dock('Data In')
        self.dockarea.addDock(self.ui.dock_data_in, 'top')

        self.ui.dock_processed = Dock('Processed Data')
        self.dockarea.addDock(self.ui.dock_processed, 'bottom', self.ui.dock_data_in)


        self.ui.dock_retriever = Dock('Retriever')
        self.dockarea.addDock(self.ui.dock_retriever, 'below', self.ui.dock_processed)

        self.ui.dock_retrieved_trace = Dock('Retrieved Trace')
        self.dockarea.addDock(self.ui.dock_retrieved_trace, 'below', self.ui.dock_retriever)

        self.ui.dock_retrieved_data = Dock('Retrieved Data')
        self.dockarea.addDock(self.ui.dock_retrieved_data, 'below', self.ui.dock_retrieved_trace)

        self.ui.dock_processed.raiseDock()

        # ######################################################
        #  setup data in dock
        self.data_in_toolbar = QtWidgets.QToolBar()
        self.load_trace_in_action = gutils.QAction(QIcon(QPixmap(":/icons/Icon_Library/Open_2D.png")),
                                                   'Load Experimental Trace')
        self.load_spectrum_in_action = gutils.QAction(QIcon(QPixmap(":/icons/Icon_Library/Open_1D.png")),
                                                      'Load Experimental Spectrum')
        self.gen_trace_in_action = gutils.QAction(QIcon(QPixmap(":/icons/Icon_Library/ini.png")),
                                           'Simulate Experimental Trace')
        self.load_from_simulation_action = gutils.QAction(QIcon(QPixmap(":/icons/Icon_Library/Open_sim.png")),
                                                  'Load Data from Simulation')

        self.load_trace_in_action.triggered.connect(self.load_trace_in)
        self.load_spectrum_in_action.triggered.connect(self.load_spectrum_in)
        self.gen_trace_in_action.triggered.connect(self.open_simulator)
        self.load_from_simulation_action.triggered.connect(self.load_from_simulator)

        self.data_in_toolbar.addAction(self.load_trace_in_action)
        self.data_in_toolbar.addAction(self.load_spectrum_in_action)
        self.data_in_toolbar.addAction(self.gen_trace_in_action)
        self.data_in_toolbar.addAction(self.load_from_simulation_action)
        data_in_splitter = QtWidgets.QSplitter()
        self.viewer_trace_in = Viewer2D()
        self.viewer_trace_in.ui.histogram_red.gradient.restoreState(Gradients['femto'])
        self.viewer_trace_in.aspect_ratio_action.click()
        self.viewer_spectrum_in = Viewer1D()
        self.settings_data_in_tree = ParameterTree()
        self.settings_data_in_tree.setMinimumWidth(300)

        data_in_splitter.addWidget(self.viewer_trace_in.parent)
        data_in_splitter.addWidget(self.viewer_spectrum_in.parent)
        data_in_splitter.addWidget(self.settings_data_in_tree)
        self.ui.dock_data_in.addWidget(self.data_in_toolbar)
        self.ui.dock_data_in.addWidget(data_in_splitter)

        self.settings_data_in_tree.setParameters(self.settings_data_in, showTop=False)
        self.settings_data_in.sigTreeStateChanged.connect(self.data_in_settings_changed)

        # #################################################
        # setup retriever dock
        retriever_widget = QtWidgets.QSplitter()
        self.viewer_live_trace = Viewer2D()
        self.viewer_live_trace.ui.histogram_red.gradient.restoreState(Gradients['femto'])
        self.viewer_error = Viewer0D()
        self.settings_retriever_tree = ParameterTree()
        self.settings_retriever_tree.setMinimumWidth(300)
        self.settings_retriever_tree.setParameters(self.settings_retriever, showTop=False)

        self.ui.dock_retriever.addWidget(retriever_widget)
        retriever_widget.addWidget(self.viewer_live_trace.parent)
        retriever_widget.addWidget(self.viewer_error.parent)
        retriever_widget.addWidget(self.settings_retriever_tree)

        #####################################
        # setup processed dock
        main_widget = QtWidgets.QWidget()
        main_widget.setLayout(QtWidgets.QHBoxLayout())
        trace_widget = QtWidgets.QWidget()
        pulse_widget = QtWidgets.QWidget()
        main_widget.layout().addWidget(trace_widget)
        main_widget.layout().addWidget(pulse_widget)
        self.ui.dock_processed.addWidget(main_widget)

        self.pulse_canvas = MplCanvas(pulse_widget, width=5, height=4, dpi=100)
        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        toolbar_pulse = NavigationToolbar(self.pulse_canvas, trace_widget)

        self.trace_canvas = MplCanvas(trace_widget, width=5, height=4, dpi=100)
        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        toolbar_trace = NavigationToolbar(self.trace_canvas, trace_widget)

        pulse_widget.setLayout(QtWidgets.QVBoxLayout())
        pulse_widget.layout().addWidget(toolbar_pulse)
        pulse_widget.layout().addWidget(self.pulse_canvas)
        trace_widget.setLayout(QtWidgets.QVBoxLayout())
        trace_widget.layout().addWidget(toolbar_trace)
        trace_widget.layout().addWidget(self.trace_canvas)

        ##################################################
        # setup retrieved trace dock



        ##################################################
        # setup retrievd data dock

    def open_simulator(self):
        simulator_widget = QtWidgets.QWidget()
        self.simulator = Simulator(simulator_widget)
        simulator_widget.setWindowTitle('PyMoDAQ Femto Simulator')
        simulator_widget.show()

    def load_from_simulator(self):
        if self.simulator is not None:
            data, axis, parameter_axis = self.simulator.trace_exp(Npts=512)
            spectrum_axis, spectrum = self.simulator.spectrum_exp(Npts=512)
            if self.data_in is None:
                self.data_in = DataIn(source='simulated')

            self.data_in.update(dict(source='simulated',
                                     raw_spectrum=dict(data=spectrum, axis=spectrum_axis),
                                     raw_trace=dict(data=data, axis=axis, parameter_axis=parameter_axis)))

            self.display_data_in()
            self.update_spectrum_info(self.data_in['raw_spectrum'])
            self.update_trace_info(self.data_in['raw_trace'])

    def update_spectrum_info(self, raw_spectrum):
            wl0, fwhm = utils.my_moment(raw_spectrum['axis']['data'], raw_spectrum['data'])
            self.settings_data_in.child('data_in_info', 'spectrum_in_info', 'wl0').setValue(wl0 * 1e9)
            self.settings_data_in.child('data_in_info', 'spectrum_in_info', 'wl_fwhm').setValue(fwhm * 1e9)
            self.settings_data_in.child('data_in_info', 'spectrum_in_info',
                                        'spectrum_size').setValue(len(raw_spectrum['data']))

    def update_trace_info(self, raw_trace):
        wl0, fwhm = utils.my_moment(raw_trace['axis']['data'], np.sum(raw_trace['data'], 0))
        self.settings_data_in.child('data_in_info', 'trace_in_info', 'wl0').setValue(wl0 * 1e9)
        self.settings_data_in.child('data_in_info', 'trace_in_info', 'wl_fwhm').setValue(fwhm * 1e9)

        self.settings_data_in.child('data_in_info', 'trace_in_info', 'trace_param_size').setValue(
            len(raw_trace['parameter_axis']))
        self.settings_data_in.child('data_in_info', 'trace_in_info',
                                    'trace_wl_size').setValue(len(raw_trace['axis']['data']))

        self.settings_data_in.child('processing', 'grid_settings', 'wl0').setValue(wl0 * 1e9)
        self.settings_data_in.child('processing', 'grid_settings',
                                    'npoints').setValue(len(raw_trace['parameter_axis']['data']))

        method = self.settings_data_in.child('data_in_info', 'trace_type').value()
        if not (method == 'dscan' or method == 'miips'):
            self.settings_data_in.child('processing', 'grid_settings',
                                        'time_resolution').setValue(np.mean(
                np.diff(raw_trace['parameter_axis']['data'])) * 1e15)

    def generate_ft_grid(self):
        wl0 = self.settings_data_in.child('processing', 'grid_settings', 'wl0').value() * 1e-9
        Npts = self.settings_data_in.child('processing', 'grid_settings', 'npoints').value()
        dt = self.settings_data_in.child('processing', 'grid_settings', 'time_resolution').value() * 1e-15
        self.ft = FourierTransform(Npts, dt, w0=wl2om(-wl0 - 300e-9))


    def process_trace(self):
        trace_in = self.get_trace_in()
        ## substract dark range of the spectrometer (if any)
        if len(self.viewer_trace_in.roi_manager.ROIs) > 0:
            roi = [self.viewer_trace_in.roi_manager.ROIs.keys()][0]
            pos = self.viewer_trace_in.roi_manager.ROIs[roi].pos()
            width, height = self.viewer_trace_in.roi_manager.ROIs[roi].size()
            xlim_pxls = np.array([pos.x(), pos.y()+width])
            ylim_pxls = np.array([pos.x(), pos.y() + height])
            xlim, ylim = self.viewer_trace_in.scale_axis(xlim_pxls, ylim_pxls)
            trace_in = preprocess(trace_in, signal_range=None, dark_signal_range=tuple(xlim))

        if self.viewer_trace_in.ROIselect_action.isChecked():
            x0 = self.settings_data_in.child('processing', 'ROIselect', 'x0').value()
            y0 = self.settings_data_in.child('processing', 'ROIselect', 'y0').value()
            width = self.settings_data_in.child('processing', 'ROIselect', 'width').value()
            height = self.settings_data_in.child('processing', 'ROIselect', 'height').value()
            xlim_pxls = np.array([x0, x0+width])
            ylim_pxls = np.array([y0, y0+height])
            xlim, ylim = self.viewer_trace_in.scale_axis(xlim_pxls, ylim_pxls)
            trace_in = preprocess(trace_in, signal_range=(tuple(ylim), tuple(xlim)))

        self.data_in['trace_in'] = trace_in
        preprocess2(self.data_in['trace_in'], self.pnps)

        self.trace_canvas.figure.clf()
        MeshDataPlot(trace_in, self.trace_canvas.figure)
        self.trace_canvas.draw()

    def process_both(self):
        self.process_trace()
        self.process_spectrum()

    def process_spectrum(self):
        self.generate_ft_grid()
        method = self.settings_data_in.child('data_in_info', 'trace_type').value()
        nlprocess = self.settings_data_in.child('data_in_info', 'nlprocess').value()
        wl0 = self.settings_data_in.child('processing', 'grid_settings', 'wl0').value() * 1e-9
        spectrum = self.data_in['raw_spectrum']['data']
        wavelength = self.data_in['raw_spectrum']['axis']['data']
        self.data_in['pulse_in'] = Pulse(self.ft, wl0)

        for roi in self.viewer_spectrum_in.roi_manager.ROIs:
            range = self.viewer_spectrum_in.roi_manager.ROIs[roi].pos()
            spectrum = mask(wavelength, spectrum, (range[0] <= wavelength) & (wavelength <= range[1]))

        self.data_in['pulse_in'] = pulse_from_spectrum(wavelength, spectrum, pulse=self.data_in['pulse_in'])
        self.pnps = PNPS(self.data_in['pulse_in'], method, nlprocess)

        self.pulse_canvas.figure.clf()
        PulsePlot(self.data_in['pulse_in'], self.pulse_canvas.figure)
        self.pulse_canvas.draw()

    def start_retriever(self):
        retriever_cls = _RETRIEVER_CLASSES[self.settings_retriever.child('algo_type').value()]
        verbose = self.settings_retriever.child('verbose').value()
        max_iter = self.settings_retriever.child('max_iter').value()
        preprocess2(self.data_in['trace_in'], self.pnps)
        if self.pnps is not None:
            self.retriever = retriever_cls(self.pnps, logging=True, verbose=verbose, maxiter=max_iter)
            random_gaussian(self.data_in['pulse_in'], 1000e-15, phase_max=0.1)
            # now retrieve from the synthetic trace simulated above
            self.retriever.retrieve(self.data_in['trace_in'], self.data_in['pulse_in'].spectrum, weights=None)

            result = self.retriever.result()
            fundamental = self.data_in['raw_spectrum']['data']
            wavelength = self.data_in['raw_spectrum']['axis']['data']
            fundamental *= (wavelength * wavelength)

            spec = self.data_in['pulse_in'].spectral_intensity
            spec = scipy.interpolate.interp1d(self.data_in['pulse_in'].wl, spec,
                                              bounds_error=False,
                                              fill_value=0.0)(wavelength)
            fundamental *= lib.best_scale(fundamental, spec)
            print("spectrum error", "%e" % lib.nrms(fundamental, spec))

            # do the retrieval plot
            RetrievalResultPlot(result, fundamental=fundamental,
                                fundamental_wavelength=wavelength,
                                oversampling=8, phase_blanking=True,
                                phase_blanking_threshold=0.01, limit=True)
        else:
            return
    @pyqtSlot(QtCore.QRectF)
    def update_ROI(self, rect=QtCore.QRectF(0, 0, 1, 1)):
        self.settings_data_in.child('processing', 'ROIselect', 'x0').setValue(int(rect.x()))
        self.settings_data_in.child('processing', 'ROIselect', 'y0').setValue(int(rect.y()))
        self.settings_data_in.child('processing', 'ROIselect', 'width').setValue(max([1, int(rect.width())]))
        self.settings_data_in.child('processing', 'ROIselect', 'height').setValue(max([1, int(rect.height())]))

    def show_ROI(self):
        self.settings_data_in.child('processing', 'ROIselect').setOpts(
            visible=self.viewer_trace_in.ROIselect_action.isChecked())
        pos = self.viewer_trace_in.ui.ROIselect.pos()
        size = self.viewer_trace_in.ui.ROIselect.size()
        self.update_ROI(QtCore.QRectF(pos[0], pos[1], size[0], size[1]))

    def get_trace_in(self):
        method = self.settings_data_in.child('data_in_info', 'trace_type').value()
        if method == 'dscan':
            label = 'Insertion'
            unit = 'm'
        elif method == 'miips':
            label = 'Phase'
            unit = 'rad'
        else:
            label = 'Delay'
            unit = 's'

        self.data_in['trace_in'] = MeshData(self.data_in['raw_trace']['data'],
                                            self.data_in['raw_trace']['parameter_axis']['data'],
                                            self.data_in['raw_trace']['axis']['data'],
                                            labels=[label, "wavelength"], units=[unit, "m"])

        return self.data_in['trace_in']

    def get_pulse_in(self):

        self.data_in['pulse_in'] = pulse_from_spectrum(self.data_in['raw_spectrum']['axis']['data'],
                                                       self.data_in['raw_spectrum']['data'])


    def get_axes_from_trace_node(self, fname, node_path):
        h5file = self.h5browse.open_file(fname)
        data, axes, nav_axes, is_spread = self.h5browse.get_h5_data(node_path)
        self.h5browse.close_file()
        return axes['x_axis'], axes['nav_00']

    def load_trace_in(self, fname=None, node_path=None):
        try:
            if fname is not None and node_path is not None:
                h5file = self.h5browse.open_file(fname)
                data, axes, nav_axes, is_spread = self.h5browse.get_h5_data(node_path)
                self.h5browse.close_file()
            else:
                data, fname, node_path = browse_data(ret_all=True, message='Select the node corresponding to the'
                                                                           'Charaterization Trace')

            if fname != '':
                wl, parameter_axis = self.get_axes_from_trace_node(fname, node_path)
                if self.data_in is None:
                    self.data_in = DataIn(source='experimental')

                scaling_parameter = self.settings_data_in.child('data_in_info',
                                                                'trace_in_info', 'param_scaling').value()
                scaling_wl = self.settings_data_in.child('data_in_info', 'trace_in_info', 'wl_scaling').value()

                wl['units'] = 'm'
                wl['data'] *= scaling_wl

                parameter_axis['data'] *= scaling_parameter
                parameter_axis['units'] = 'p.u.'

                self.data_in.update(dict(raw_trace={'data': data, 'axis': wl, 'parameter_axis': parameter_axis}))


                #
                # trace_in = MeshData(data, parameter_axis['data'] * scaling_parameter, wl['data'] * scaling_wl,
                #                          labels=[parameter_axis['label'], wl['label']],
                #                          units=['s', 'm'])
                self.update_trace_info(self.data_in['raw_trace'])

                # dt = np.mean(np.diff(parameter_axis['data'])) * scaling_parameter
                # wl0 = self.settings_data_in.child('data_in_info', 'trace_in_info', 'wl0').value() * 1e-9
                # ft = FourierTransform(len(parameter_axis['data']), dt, w0=wl2om(-wl0 - 300e-9))
                # self.data_in = DataIn(name='data_in', source='experimental', trace_in=trace_in,
                #                       pulse_in=Pulse(ft, self.settings_data_in.child('data_in_info',
                #                                                                        'trace_in_info',
                #                                                                        'wl0').value() * scaling_wl))

                self.display_trace_in()

        except Exception as e:
            logger.exception(str(e))

    def load_spectrum_in(self, fname=None, node_path=None):
        if fname is not None and node_path is not None:
            h5file = self.h5browse.open_file(fname)
            data, axes, nav_axes, is_spread = self.h5browse.get_h5_data(node_path)
            self.h5browse.close_file()

        else:
            data, fname, node_path = browse_data(ret_all=True, message='Select the node corresponding to the'
                                                                   'Fundamental Spectrum')
            if fname != '':
                h5file = self.h5browse.open_file(fname)
                data, axes, nav_axes, is_spread = self.h5browse.get_h5_data(node_path)
                self.h5browse.close_file()
            else:
                return

        if self.data_in is None:
            self.data_in = DataIn(source='experimental')
        self.data_in.update(dict(raw_spectrum={'data': data, 'axis': axes['x_axis']}))

        self.update_spectrum_info(self.data_in['raw_spectrum'])
        self.display_spectrum_in()

    def display_trace_in(self):
        self.viewer_trace_in.setImage(self.data_in['raw_trace']['data'])
        self.viewer_trace_in.x_axis = self.data_in['raw_trace']['axis']
        self.viewer_trace_in.y_axis = self.data_in['raw_trace']['parameter_axis']
        
    def display_spectrum_in(self):
        self.viewer_spectrum_in.show_data([self.data_in['raw_spectrum']['data']],
                                          x_axis=self.data_in['raw_spectrum']['axis'],
                                          labels=['Spectrum'])

    def display_data_in(self):
        self.display_trace_in()
        self.display_spectrum_in()

    def data_in_settings_changed(self, param, changes):
        for param, change, data in changes:
            if change == 'childAdded':
                pass

            elif change == 'value':
                pass

            elif change == 'parent':
                pass


def main():
    from pymodaq.daq_utils.daq_utils import get_set_preset_path
    app = QtWidgets.QApplication(sys.argv)
    win = QtWidgets.QMainWindow()
    area = gutils.DockArea()
    win.setCentralWidget(area)
    win.resize(1000, 500)
    win.setWindowTitle('PyMoDAQ Retriever')

    prog = Retriever(dashboard=None, dockarea=area)
    win.show()
    prog.load_trace_in(fname='C:\\Data\\2021\\20210309\\Dataset_20210309_000\\Dataset_20210309_000.h5',
                        node_path='/Raw_datas/Scan000/Detector000/Data1D/Ch000/Data')
    prog.load_spectrum_in(fname='C:\\Users\\weber\\Desktop\\pulse.h5',
                        node_path='/Raw_datas/Detector000/Data1D/Ch000/Data')

    sys.exit(app.exec_())


if __name__ == '__main__':
    main()


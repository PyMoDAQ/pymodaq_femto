import sys
import subprocess

import logging
from pathlib import Path
import numpy as np
from types import SimpleNamespace
import io
from contextlib import redirect_stdout
from PyQt5 import QtGui, QtWidgets, QtCore
from PyQt5.QtCore import Qt, QObject, pyqtSlot, QThread, pyqtSignal, QLocale
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtGui import QTextCursor
from pyqtgraph.dockarea import Dock
from pyqtgraph.parametertree import Parameter, ParameterTree
from pymodaq.daq_utils import daq_utils as utils
from pymodaq.daq_utils import gui_utils as gutils
from pymodaq.daq_utils.parameter import utils as putils, ioxml
from pymodaq.daq_utils.h5modules import browse_data
from pymodaq.daq_utils.plotting.viewer2D.viewer2D_main import Viewer2D
from pymodaq.daq_utils.plotting.viewer1D.viewer1D_main import Viewer1D
from pymodaq.daq_utils.plotting.viewer0D.viewer0D_main import Viewer0D
from pymodaq.daq_utils.managers.roi_manager import LinearROI
from pymodaq_femto.graphics import RetrievalResultPlot, MplCanvas, NavigationToolbar, MeshDataPlot, PulsePlot
from pymodaq_femto.simulation import Simulator, methods, nlprocesses, materials
from collections import OrderedDict
from pypret import FourierTransform, Pulse, PNPS, lib, MeshData, random_gaussian
from pypret.frequencies import om2wl, wl2om, convert
import scipy
from scipy.fftpack import next_fast_len
from pymodaq.daq_utils.h5modules import H5BrowserUtil, H5Saver
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
        dark_signal.limit(dark_signal_range, axes=1)
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

def substract_linear_phase(pulse):
    phase = np.unwrap(np.angle(pulse.spectrum))
    z = np.polyfit(pulse.w, phase, 1)
    pulse.spectrum *= np.exp(- 1j * np.poly1d(z)(pulse.w))
    return pulse
    
    
def mask(x, y, where, **kwargs):
    y = scipy.interpolate.interp1d(x[~where], y[~where], **kwargs)(x)
    return y
params_simul = Simulator.params
params_algo = utils.find_dict_in_list_from_key_val(params_simul, 'name', 'algo')

class Retriever(QObject):
    """
    Main class initializing a DAQ_Scan module with its dashboard and scanning control panel
    """
    status_signal = pyqtSignal(str)
    retriever_signal = pyqtSignal(str)
    params_in = [params_algo,
        {'title': 'Data Info', 'name': 'data_in_info', 'type': 'group', 'children': [
            {'title': 'Loaded file:', 'name': 'loaded_file', 'type': 'text', 'value': "", 'readonly': True,
             'tip': 'Loaded trace file'},
            {'title': 'Loaded node:', 'name': 'loaded_node', 'type': 'str', 'value': "", 'readonly': True,
             'tip': 'Loaded node within trace file'},
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
            {'title': 'Trace limits:', 'name': 'ROIselect', 'type': 'group', 'visible': True, 'children': [
                {'title': 'Crop Trace?:', 'name': 'crop_trace', 'type': 'bool', 'value': False},
                {'title': 'x0:', 'name': 'x0', 'type': 'int', 'value': 0, 'min': 0},
                {'title': 'y0:', 'name': 'y0', 'type': 'int', 'value': 0, 'min': 0},
                {'title': 'width:', 'name': 'width', 'type': 'int', 'value': 10, 'min': 1},
                {'title': 'height:', 'name': 'height', 'type': 'int', 'value': 10, 'min': 1},
            ]},
            {'title': 'Substract background:', 'name': 'linearselect', 'type': 'group', 'visible': True, 'children': [
                {'title': 'Substract?:', 'name': 'dosubstract', 'type': 'bool', 'value': False},
                {'title': 'wl0:', 'name': 'wl0', 'type': 'float', 'value': 0.},
                {'title': 'wl1:', 'name': 'wl1', 'type': 'float', 'value': 10.},
            ]},
            {'title': 'Process Spectrum', 'name': 'process_spectrum', 'type': 'action',
             'tip': 'Use ROIs to select frequency areas from the spectrum that should be removed'},
            {'title': 'Process trace', 'name': 'process_trace', 'type': 'action',
             'tip': 'Use one ROI to select a frequency area from the trace in order to remove the background'
                    ' and use ROISelect to reduce the area around the trace'},

            {'title': 'Process Both', 'name': 'process_both', 'type': 'action',
             'tip': 'Process both the trace and the spectrum'},
        ]},
        {'title': 'Retrieving', 'name': 'retrieving', 'type': 'group', 'children': [
            {'title': 'Algo type:', 'name': 'algo_type', 'type': 'list', 'values': retriever_algos,
             'tip': 'Retriever Algorithm'},
            {'title': 'Verbose Info:', 'name': 'verbose', 'type': 'bool', 'value': True,
             'tip': 'Display infos during retrieval'},
            {'title': 'Max iteration:', 'name': 'max_iter', 'type': 'int', 'value': 30,
             'tip': 'Max iteration for the algorithm'},
            {'title': 'Initial Pulse Guess', 'name': 'pulse_guess', 'type': 'group', 'children': [
                {'title': 'FWHM (fs):', 'name': 'fwhm', 'type': 'float', 'value': 100,
                 'tip': 'Guess of the pulse duration (used as a starting point)'},
                {'title': 'Phase amp. (rad):', 'name': 'phase_amp', 'type': 'float', 'value': 0.1,
                 'tip': 'Amplitude of the random phase applied to the initial guess'},
            ]},

            {'title': 'Start Retrieval', 'name': 'start', 'type': 'action',
             'tip': 'Start the retrieval process'},
            {'title': 'Stop Retrieval', 'name': 'stop', 'type': 'action',
             'tip': 'Stop the retrieval process'},
           ]},
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

        self.settings = Parameter.create(name='dataIN_settings', type='group', children=self.params_in)
        self.settings.sigTreeStateChanged.connect(self.settings_changed)

        self.setupUI()
        self.create_menu(self.mainwindow.menuBar())
        self.simulator = None
        self.data_in = None
        self.ft = None
        self.retriever = None
        self.pnps = None
        self.retriever_thread = None
        self.save_file_pathname = None
        self.settings.child('processing', 'process_trace').sigActivated.connect(self.process_trace)
        self.settings.child('processing', 'process_spectrum').sigActivated.connect(self.process_spectrum)
        self.settings.child('processing', 'process_both').sigActivated.connect(self.process_both)
        self.settings.child('retrieving', 'start').sigActivated.connect(self.start_retriever)
        self.settings.child('retrieving', 'stop').sigActivated.connect(self.stop_retriever)

        self.viewer_trace_in.ROI_select_signal.connect(self.update_ROI)
        self.viewer_trace_in.ROIselect_action.triggered.connect(self.show_ROI)

        self.settings.child('algo', 'miips_parameter').hide()
        self.settings.child('algo', 'dscan_parameter').hide()

    def save_data(self, save_file_pathname=None):
        try:
            if save_file_pathname is None:
                save_file_pathname = gutils.select_file(start_path=self.save_file_pathname, save=True,
                                                        ext='h5')  # see daq_utils
            h5saver = H5Saver(save_type='custom')
            h5saver.init_file(update_h5=True, custom_naming=False, addhoc_file_path=save_file_pathname,
                              raw_group_name='PyMoDAQFemtoAnalysis')

            settings_str = b'<All_settings>' + ioxml.parameter_to_xml_string(self.settings)
            settings_str += b'</All_settings>'
            # TODO
            h5saver.set_attr(h5saver.raw_group, 'settings', settings_str)
            data_in_group = h5saver.get_set_group(h5saver.raw_group, "DataIn")
            trace_group = h5saver.get_set_group(data_in_group, 'NLTrace')
            spectrum_group = h5saver.get_set_group(data_in_group, 'FunSpectrum')
            h5saver.add_data(trace_group, self.data_in['raw_trace'], scan_type='')
            h5saver.add_data(spectrum_group, self.data_in['raw_spectrum'], scan_type='')

        except Exception as e:
            pass

        h5saver.close_file()

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
        self.data_in_menu.addAction(self.load_spectrum_in_action)
        self.data_in_menu.addSeparator()
        self.data_in_menu.addAction(self.gen_trace_in_action)
        self.data_in_menu.addAction(self.load_from_simulation_action)

    def settings_changed(self, param, changes):
        for param, change, data in changes:
            path = self.settings.childPath(param)
            if change == 'childAdded':
                pass
            elif change == 'parent':
                pass
            elif change == 'value':
                if param.name() == 'method':
                    self.settings.child('algo', 'nlprocess').setLimits(list(_PNPS_CLASSES[param.value()].keys()))

                    if param.value() == 'miips':
                        self.settings.child('algo', 'alpha').show()
                        self.settings.child('algo', 'gamma').show()
                    else:
                        self.settings.child('algo', 'alpha').hide()
                        self.settings.child('algo', 'gamma').hide()
                    if param.value() == 'dscan':
                        self.settings.child('algo', 'material').show()
                    else:
                        self.settings.child('algo', 'material').hide()


                elif param.name() in putils.iter_children(self.settings.child('processing', 'ROIselect'),
                                                          []) and 'ROIselect' in param.parent().name():  # to be sure
                    # a param named 'y0' for instance will not collide with the y0 from the ROI
                    try:
                        self.viewer_trace_in.ROI_select_signal.disconnect(self.update_ROI)
                    except Exception as e:
                        pass
                    if self.settings.child('processing', 'ROIselect', 'crop_trace').value():
                        if not self.viewer_trace_in.ROIselect_action.isChecked():
                            self.viewer_trace_in.ROIselect_action.trigger()
                            QtWidgets.QApplication.processEvents()
                        self.viewer_trace_in.ui.ROIselect.setPos(
                            self.settings.child('processing', 'ROIselect', 'x0').value(),
                            self.settings.child('processing', 'ROIselect', 'y0').value())
                        self.viewer_trace_in.ui.ROIselect.setSize(
                            [self.settings.child('processing', 'ROIselect', 'width').value(),
                             self.settings.child('processing', 'ROIselect', 'height').value()])
                        self.viewer_trace_in.ROI_select_signal.connect(self.update_ROI)
                    else:
                        if self.viewer_trace_in.ROIselect_action.isChecked():
                            self.viewer_trace_in.ROIselect_action.trigger()
                elif param.name() in putils.iter_children(self.settings.child('processing', 'linearselect'),
                                                          []) and 'linearselect' in param.parent().name():  # to be sure
                    # a param named 'y0' for instance will not collide with the y0 from the ROI
                    try:
                        self.linear_region.sigRegionChangeFinished.disconnect(self.update_linear)
                    except Exception as e:
                        pass
                    self.linear_region.setVisible(self.settings.child('processing',
                                                                              'linearselect', 'dosubstract').value())
                    pos_real = np.array([self.settings.child('processing', 'linearselect', 'wl0').value(),
                                         self.settings.child('processing',
                                                                     'linearselect', 'wl1').value()]) * 1e-9

                    pos_pxl, y = self.viewer_trace_in.unscale_axis(np.array(pos_real), np.array([0, 1]))
                    self.linear_region.setPos(pos_pxl)
                    self.linear_region.sigRegionChangeFinished.connect(self.update_linear)
                    
                elif param.name() == 'method':
                    self.settings.child('algo', 'nlprocess').setLimits(list(_PNPS_CLASSES[param.value()].keys()))

                    if param.value() == 'miips':
                        self.settings.child('algo', 'alpha').show()
                        self.settings.child('algo', 'gamma').show()
                        self.settings.child('algo', 'miips_parameter').show()
                    else:
                        self.settings.child('algo', 'alpha').hide()
                        self.settings.child('algo', 'gamma').hide()
                        self.settings.child('algo', 'miips_parameter').hide()

                    if param.value() == 'dscan':
                        self.settings.child('algo', 'material').show()
                        self.settings.child('algo', 'dscan_parameter').show()
                    else:
                        self.settings.child('algo', 'material').hide()
                        self.settings.child('algo', 'dscan_parameter').hide()


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

        self.ui.dock_settings = Dock('Settings')
        self.dockarea.addDock(self.ui.dock_settings, 'top')

        self.ui.dock_data_in = Dock('Data In')
        self.dockarea.addDock(self.ui.dock_data_in, 'left', self.ui.dock_settings)

        self.ui.dock_processed = Dock('Processed Data')
        self.dockarea.addDock(self.ui.dock_processed, 'below', self.ui.dock_data_in)


        self.ui.dock_retriever = Dock('Retriever')
        self.dockarea.addDock(self.ui.dock_retriever, 'below', self.ui.dock_processed)


        self.ui.dock_retrieved_data = Dock('Retrieved Data')
        self.dockarea.addDock(self.ui.dock_retrieved_data, 'below', self.ui.dock_retriever)

        self.ui.dock_processed.raiseDock()


        # ######################################################
        #  setup settings in dock
        self.settings_tree = ParameterTree()
        self.settings_tree.setMinimumWidth(300)
        self.ui.dock_settings.addWidget(self.settings_tree)
        self.ui.dock_settings.setStretch(0.5)

        # setup toolbar
        self.toolbar = QtWidgets.QToolBar()
        self.mainwindow.addToolBar(self.toolbar)
        self.load_trace_in_action = gutils.QAction(QIcon(QPixmap(":/icons/Icon_Library/Open_2D.png")),
                                                   'Load Experimental Trace')
        self.load_spectrum_in_action = gutils.QAction(QIcon(QPixmap(":/icons/Icon_Library/Open_1D.png")),
                                                      'Load Experimental Spectrum')
        self.gen_trace_in_action = gutils.QAction(QIcon(QPixmap(":/icons/Icon_Library/ini.png")),
                                           'Simulate Experimental Trace')
        self.load_from_simulation_action = gutils.QAction(QIcon(QPixmap(":/icons/Icon_Library/Open_sim.png")),
                                                  'Load Data from Simulation')

        self.save_data_action = gutils.QAction(QIcon(QPixmap(":/icons/Icon_Library/Save.png")),
                                                  'Save Data')

        self.load_trace_in_action.triggered.connect(self.load_trace_in)
        self.load_spectrum_in_action.triggered.connect(self.load_spectrum_in)
        self.gen_trace_in_action.triggered.connect(self.open_simulator)
        self.load_from_simulation_action.triggered.connect(self.load_from_simulator)
        self.save_data_action.triggered.connect(self.save_data)

        self.toolbar.addAction(self.load_trace_in_action)
        self.toolbar.addAction(self.load_spectrum_in_action)
        self.toolbar.addAction(self.gen_trace_in_action)
        self.toolbar.addAction(self.load_from_simulation_action)
        self.toolbar.addSeparator()
        self.toolbar.addAction(self.save_data_action)
        
        
        # ######################################################
        #  setup data in dock
        
        data_in_splitter = QtWidgets.QSplitter()
        self.viewer_trace_in = Viewer2D()
        self.viewer_trace_in.ui.histogram_red.gradient.restoreState(Gradients['femto'])
        self.viewer_trace_in.aspect_ratio_action.click()

        pos = self.viewer_trace_in.roi_manager.viewer_widget.plotItem.vb.viewRange()[0]
        self.linear_region = LinearROI(index=0, pos=pos)
        self.linear_region.setZValue(-10)
        self.linear_region.setOpacity(0.5)
        self.linear_region.setBrush([255, 0, 0, 0])
        self.viewer_trace_in.roi_manager.viewer_widget.plotItem.addItem(self.linear_region)
        self.linear_region.sigRegionChangeFinished.connect(self.update_linear)
        self.linear_region.setVisible(False)

        self.viewer_spectrum_in = Viewer1D()

        data_in_splitter.addWidget(self.viewer_trace_in.parent)
        data_in_splitter.addWidget(self.viewer_spectrum_in.parent)
        self.ui.dock_data_in.addWidget(data_in_splitter)

        self.settings_tree.setParameters(self.settings, showTop=False)

        # #################################################
        # setup retriever dock
        retriever_widget = QtWidgets.QSplitter()
        self.viewer_live_trace = Viewer2D()
        self.viewer_live_trace.ui.histogram_red.gradient.restoreState(Gradients['femto_error'])
        self.viewer_live_trace.aspect_ratio_action.trigger()
        #self.viewer_live_trace.auto_levels_action.trigger()
        self.viewer_live_time = Viewer1D()
        self.viewer_live_lambda = Viewer1D()
        self.info_widget = QtWidgets.QTextEdit()
        self.info_widget.setReadOnly(True)

        vsplitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.ui.dock_retriever.addWidget(retriever_widget)
        retriever_widget.addWidget(self.viewer_live_trace.parent)
        vsplitter.addWidget(self.viewer_live_time.parent)
        vsplitter.addWidget(self.viewer_live_lambda.parent)
        retriever_widget.addWidget(vsplitter)
        retriever_widget.addWidget(self.info_widget)

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
        # setup retrievd data dock
        data_widget = QtWidgets.QWidget()
        data_widget.setLayout(QtWidgets.QVBoxLayout())
        self.data_canvas = MplCanvas(data_widget, width=5, height=4, dpi=100)
        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        toolbar_data = NavigationToolbar(self.data_canvas, data_widget)
        data_widget.layout().addWidget(toolbar_data)
        data_widget.layout().addWidget(self.data_canvas)
        self.ui.dock_retrieved_data.addWidget(data_widget)

        self.ui.dock_data_in.raiseDock()

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

            for child in putils.iter_children_params(self.settings.child('algo'), childlist=[]):
                path = ['algo']
                path.extend(self.settings.child('algo').childPath(child))
                self.settings.child(*path).setValue(self.simulator.settings.child(*path).value())
            self.settings.child('data_in_info', 'loaded_file').setValue('Simulation')

    def update_spectrum_info(self, raw_spectrum):
        wl0, fwhm = utils.my_moment(raw_spectrum['x_axis']['data'], raw_spectrum['data'])
        self.settings.child('data_in_info', 'spectrum_in_info', 'wl0').setValue(wl0 * 1e9)
        self.settings.child('data_in_info', 'spectrum_in_info', 'wl_fwhm').setValue(fwhm * 1e9)
        self.settings.child('data_in_info', 'spectrum_in_info',
                                    'spectrum_size').setValue(len(raw_spectrum['data']))

        self.settings.child('processing', 'grid_settings', 'wl0').setValue(wl0 * 1e9)

    def update_trace_info(self, raw_trace):
        wl0, fwhm = utils.my_moment(raw_trace['x_axis']['data'], np.sum(raw_trace['data'], 0))
        self.settings.child('data_in_info', 'trace_in_info', 'wl0').setValue(wl0 * 1e9)
        self.settings.child('data_in_info', 'trace_in_info', 'wl_fwhm').setValue(fwhm * 1e9)

        self.settings.child('data_in_info', 'trace_in_info', 'trace_param_size').setValue(
            len(raw_trace['y_axis']))
        self.settings.child('data_in_info', 'trace_in_info',
                                    'trace_wl_size').setValue(len(raw_trace['x_axis']['data']))


        self.settings.child('processing', 'grid_settings',
                                    'npoints').setValue(len(raw_trace['y_axis']['data']))

        method = self.settings.child('algo', 'method').value()
        if not (method == 'dscan' or method == 'miips'):
            self.settings.child('processing', 'grid_settings',
                                        'time_resolution').setValue(np.mean(
                np.diff(raw_trace['y_axis']['data'])) * 1e15)

    def generate_ft_grid(self):
        wl0 = self.settings.child('processing', 'grid_settings', 'wl0').value() * 1e-9
        Npts = self.settings.child('processing', 'grid_settings', 'npoints').value()
        dt = self.settings.child('processing', 'grid_settings', 'time_resolution').value() * 1e-15
        self.ft = FourierTransform(Npts, dt, w0=wl2om(-wl0 - 300e-9))


    def process_trace(self):
        self.ui.dock_processed.raiseDock()
        if self.pnps is None:
            logger.info('PNPS is not yet defined, process the spectrum first')
            return
        trace_in = self.get_trace_in()
        # TODO
        # ## substract bright spots range (if any)
        # if len(self.viewer_trace_in.roi_manager.ROIs) > 0:
        #     roi = [self.viewer_trace_in.roi_manager.ROIs.keys()][0]
        #     pos = self.viewer_trace_in.roi_manager.ROIs[roi].pos()
        #     width, height = self.viewer_trace_in.roi_manager.ROIs[roi].size()
        #     xlim_pxls = np.array([pos.x(), pos.y()+width])
        #     ylim_pxls = np.array([pos.x(), pos.y() + height])
        #     xlim, ylim = self.viewer_trace_in.scale_axis(xlim_pxls, ylim_pxls)
        #     trace_in = preprocess(trace_in, signal_range=None, bright_signal_range=tuple(xlim)) # bright_signal_range doesn't exists yet'

        if self.settings.child('processing', 'linearselect', 'dosubstract').value():
            xlim = np.array((self.settings.child('processing', 'linearselect', 'wl0').value(),
                    self.settings.child('processing', 'linearselect', 'wl1').value())) * 1e-9
            trace_in = preprocess(trace_in, signal_range=None, dark_signal_range=tuple(xlim))

        if self.settings.child('processing', 'ROIselect', 'crop_trace').value():
            x0 = self.settings.child('processing', 'ROIselect', 'x0').value()
            y0 = self.settings.child('processing', 'ROIselect', 'y0').value()
            width = self.settings.child('processing', 'ROIselect', 'width').value()
            height = self.settings.child('processing', 'ROIselect', 'height').value()
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
        self.process_spectrum()
        self.process_trace()


    def process_spectrum(self):

        self.ui.dock_processed.raiseDock()

        self.generate_ft_grid()
        method = self.settings.child('algo', 'method').value()
        nlprocess = self.settings.child('algo', 'nlprocess').value()
        wl0 = self.settings.child('data_in_info', 'trace_in_info', 'wl0').value() * 1e-9
        spectrum = self.data_in['raw_spectrum']['data']
        wavelength = self.data_in['raw_spectrum']['x_axis']['data']

        if 'shg' in nlprocess:
            wl0real = 2 * wl0
        elif 'thg' in nlprocess:
            wl0real = 3 * wl0
        else:
            wl0real = wl0

        self.data_in['pulse_in'] = Pulse(self.ft, wl0real)

        for roi in self.viewer_spectrum_in.roi_manager.ROIs:
            range = self.viewer_spectrum_in.roi_manager.ROIs[roi].pos()
            spectrum = mask(wavelength, spectrum, (range[0] <= wavelength) & (wavelength <= range[1]))

        self.data_in['pulse_in'] = pulse_from_spectrum(wavelength, spectrum, pulse=self.data_in['pulse_in'])
       #self.pnps = PNPS(self.data_in['pulse_in'], method, nlprocess)

        if method == 'dscan':
            material = materials[self.settings.child('algo', 'material').value()]
            self.pnps = PNPS(self.data_in['pulse_in'], method, nlprocess, material=material)
            parameter = utils.linspace_step(self.settings.child('algo', 'dscan_parameter', 'min').value(),
                                      self.settings.child('algo', 'dscan_parameter', 'max').value(),
                                      self.settings.child('algo', 'dscan_parameter', 'step').value())
            parameter *= 1e-3
        elif method == 'miips':
            alpha = self.settings.child('algo', 'alpha').value()
            gamma = self.settings.child('algo', 'gamma').value()
            self.pnps = PNPS(self.data_in['pulse_in'], method, nlprocess, alpha=alpha, gamma=gamma)
            parameter = utils.linspace_step(self.settings.child('algo', 'miips_parameter', 'min').value(),
                                      self.settings.child('algo', 'miips_parameter', 'max').value(),
                                      self.settings.child('algo', 'miips_parameter', 'step').value())
        else:
            self.pnps = PNPS(self.data_in['pulse_in'], method, nlprocess)


        self.pulse_canvas.figure.clf()
        PulsePlot(self.data_in['pulse_in'], self.pulse_canvas.figure)
        self.pulse_canvas.draw()

    def start_retriever(self):
        self.ui.dock_retriever.raiseDock()
        self.info_widget.clear()
        # mandatory to deal with multithreads
        if self.retriever_thread is not None:
            if self.retriever_thread.isRunning():
                self.retriever_thread.terminate()
                while not self.retriever_thread.isFinished():
                    QThread.msleep(100)
                self.retriever_thread = None

        self.retriever_thread = QThread()
        retriever = RetrieverWorker(self.data_in, self.pnps, self.settings)
        retriever.moveToThread(self.retriever_thread)
        retriever.status_sig[str].connect(self.update_retriever_info)
        retriever.result_signal[SimpleNamespace].connect(self.display_results)
        retriever.callback_sig[list].connect(self.update_retriever)
        self.retriever_signal[str].connect(retriever.command_retriever)
        self.retriever_thread.retriever = retriever
        self.retriever_thread.start()

        self.retriever_signal.emit('start')

    def stop_retriever(self):
        self.retriever_signal.emit('stop')

    def update_retriever_info(self, info):
        self.info_widget.moveCursor(QTextCursor.End)
        self.info_widget.insertPlainText(info+'\n')
        self.info_widget.moveCursor(QTextCursor.End)

    @pyqtSlot(list)
    def update_retriever(self, args):
        max = 0.8*np.max([np.abs(np.max(args[0])), np.abs(np.min(args[0]))])
        self.viewer_live_trace.ui.histogram_red.setHistogramRange(-max, max)
        self.viewer_live_trace.ui.histogram_red.setLevels(-max, max)
        self.viewer_live_trace.setImage(args[0])
        self.viewer_live_trace.x_axis = utils.Axis(data=args[2], label='Time', unit='s')
        self.viewer_live_trace.y_axis = utils.Axis(data=args[1], label='Frequency', unit='m')

        self.data_in['pulse_in'].spectrum = args[3]
        #self.data_in['pulse_in'] = substract_linear_phase(self.data_in['pulse_in'])
        self.viewer_live_time.show_data([np.abs(self.data_in['pulse_in'].field)**2],
                                        x_axis=utils.Axis(data=self.data_in['pulse_in'].t, label='Time', unit='s'),
                                        labels=['Temporal Intensity'])
        self.viewer_live_lambda.show_data([np.abs(self.data_in['pulse_in'].spectrum)**2],
                                        x_axis=utils.Axis(data=self.data_in['pulse_in'].wl, label='Wavelength',
                                                          unit='m'),
                                        labels=['Spectral Intensity'])

    @pyqtSlot(SimpleNamespace)
    def display_results(self, result):
        self.result = result
        self.ui.dock_retrieved_data.raiseDock()

        self.data_in['pulse_in'].spectrum = result.pulse_retrieved
        fundamental = self.data_in['raw_spectrum']['data']
        wavelength = self.data_in['raw_spectrum']['x_axis']['data']
        fundamental *= (wavelength * wavelength)
        spec = self.data_in['pulse_in'].spectral_intensity
        spec = scipy.interpolate.interp1d(self.data_in['pulse_in'].wl, spec,
                                          bounds_error=False,
                                          fill_value=0.0)(wavelength)

        fundamental *= lib.best_scale(fundamental, spec)
        print("spectrum error", "%e" % lib.nrms(fundamental, spec))

        # do the retrieval plot

        self.data_canvas.figure.clf()
        RetrievalResultPlot(result, fig=self.data_canvas.figure, fundamental=fundamental,
                            fundamental_wavelength=wavelength,
                            oversampling=8, phase_blanking=True,
                            phase_blanking_threshold=0.01, limit=True)
        self.data_canvas.draw()

    @pyqtSlot(QtCore.QRectF)
    def update_ROI(self, rect=QtCore.QRectF(0, 0, 1, 1)):
        self.settings.child('processing', 'ROIselect', 'x0').setValue(int(rect.x()))
        self.settings.child('processing', 'ROIselect', 'y0').setValue(int(rect.y()))
        self.settings.child('processing', 'ROIselect', 'width').setValue(max([1, int(rect.width())]))
        self.settings.child('processing', 'ROIselect', 'height').setValue(max([1, int(rect.height())]))

    def update_linear(self, linear_roi):
        pos = linear_roi.pos()
        pos_real, y = self.viewer_trace_in.scale_axis(np.array(pos), np.array([0, 1]))
        pos_real *= 1e9
        self.settings.child('processing', 'linearselect', 'wl0').setValue(pos_real[0])
        self.settings.child('processing', 'linearselect', 'wl1').setValue(pos_real[1])

    def show_ROI(self):
        # self.settings.child('processing', 'ROIselect').setOpts(
        #     visible=self.viewer_trace_in.ROIselect_action.isChecked())
        data = self.data_in['raw_trace']['data']
        axes = [np.arange(0, data.shape[0]), np.arange(0, data.shape[1])]
        axes_index = list(range(data.ndim))
        marginals = lib.marginals(data)
        limits = []
        for index in axes_index:
            limit = lib.limit(axes[index], marginals[index],
                              threshold=1e-2, padding=0.25)
            limits.append(limit)

        self.viewer_trace_in.ui.ROIselect.setPos((limits[1][0], limits[0][0]))
        self.viewer_trace_in.ui.ROIselect.setSize((limits[1][1]-limits[1][0], limits[0][1] - limits[0][0]))

        self.linear_region.setPos(limits[1])

        pos = self.viewer_trace_in.ui.ROIselect.pos()
        size = self.viewer_trace_in.ui.ROIselect.size()
        self.update_ROI(QtCore.QRectF(pos[0], pos[1], size[0], size[1]))

    def get_trace_in(self):
        method = self.settings.child('algo', 'method').value()
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
                                            self.data_in['raw_trace']['y_axis']['data'],
                                            self.data_in['raw_trace']['x_axis']['data'],
                                            labels=[label, "wavelength"], units=[unit, "m"])

        return self.data_in['trace_in']

    def get_pulse_in(self):

        self.data_in['pulse_in'] = pulse_from_spectrum(self.data_in['raw_spectrum']['x_axis']['data'],
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
                self.save_file_pathname = fname
                self.settings.child('data_in_info', 'loaded_file').setValue(fname)
                self.settings.child('data_in_info', 'loaded_node').setValue(node_path)
                wl, parameter_axis = self.get_axes_from_trace_node(fname, node_path)
                if self.data_in is None:
                    self.data_in = DataIn(source='experimental')

                scaling_parameter = self.settings.child('data_in_info',
                                                                'trace_in_info', 'param_scaling').value()
                scaling_wl = self.settings.child('data_in_info', 'trace_in_info', 'wl_scaling').value()

                wl['units'] = 'm'
                wl['data'] *= scaling_wl

                parameter_axis['data'] *= scaling_parameter
                parameter_axis['units'] = 'p.u.'

                self.data_in.update(dict(raw_trace={'data': data, 'x_axis': wl, 'y_axis': parameter_axis},
                                         file_path=fname,
                                         node_path=node_path))


                #
                # trace_in = MeshData(data, parameter_axis['data'] * scaling_parameter, wl['data'] * scaling_wl,
                #                          labels=[parameter_axis['label'], wl['label']],
                #                          units=['s', 'm'])
                self.update_trace_info(self.data_in['raw_trace'])

                # dt = np.mean(np.diff(parameter_axis['data'])) * scaling_parameter
                # wl0 = self.settings.child('data_in_info', 'trace_in_info', 'wl0').value() * 1e-9
                # ft = FourierTransform(len(parameter_axis['data']), dt, w0=wl2om(-wl0 - 300e-9))
                # self.data_in = DataIn(name='data_in', source='experimental', trace_in=trace_in,
                #                       pulse_in=Pulse(ft, self.settings.child('data_in_info',
                #                                                                        'trace_in_info',
                #                                                                        'wl0').value() * scaling_wl))

                self.display_trace_in()
                self.viewer_trace_in.ROIselect_action.trigger()

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
        self.data_in.update(dict(raw_spectrum={'data': data, 'x_axis': axes['x_axis']}))

        self.update_spectrum_info(self.data_in['raw_spectrum'])
        self.display_spectrum_in()

    def display_trace_in(self):
        self.viewer_trace_in.setImage(self.data_in['raw_trace']['data'])
        self.viewer_trace_in.x_axis = self.data_in['raw_trace']['x_axis']
        self.viewer_trace_in.y_axis = self.data_in['raw_trace']['y_axis']
        
    def display_spectrum_in(self):
        self.viewer_spectrum_in.show_data([self.data_in['raw_spectrum']['data']],
                                          x_axis=self.data_in['raw_spectrum']['x_axis'],
                                          labels=['Spectrum'])

    def display_data_in(self):
        self.display_trace_in()
        self.display_spectrum_in()


class RetrieverWorker(QObject):
    result_signal = pyqtSignal(SimpleNamespace)
    status_sig = pyqtSignal(str)
    callback_sig = pyqtSignal(list)

    def __init__(self, data_in, pnps, settings):
        super().__init__()
        self.settings = settings
        self.data_in = data_in
        self.pnps = pnps
        self.retirever = None

    # def send_callback(self, pnps):
    #     self.callback_sig.emit([pnps.Tmn, [pnps.parameter, pnps.process_w], pnps.pulse.field, pnps.field.t])

    @pyqtSlot(str)
    def command_retriever(self, command):
        if command == 'start':
            self.start_retriever()
        elif command == 'stop':
            self.stop_retriever()

    def start_retriever(self):
        retriever_cls = _RETRIEVER_CLASSES[self.settings.child('retrieving', 'algo_type').value()]
        verbose = self.settings.child('retrieving', 'verbose').value()
        max_iter = self.settings.child('retrieving', 'max_iter').value()
        fwhm = self.settings.child('retrieving', 'pulse_guess', 'fwhm').value()
        amplitude = self.settings.child('retrieving', 'pulse_guess', 'phase_amp').value()


        preprocess2(self.data_in['trace_in'], self.pnps)


        self.retriever = retriever_cls(self.pnps, logging=True, verbose=verbose, maxiter=max_iter,
                                       status_sig=self.status_sig, callback=self.callback_sig.emit,
                                       step_command=QtWidgets.QApplication.processEvents)
        random_gaussian(self.data_in['pulse_in'], fwhm*1e-15, phase_max=amplitude)
        # now retrieve from the synthetic trace simulated above
        self.retriever.retrieve(self.data_in['trace_in'], self.data_in['pulse_in'].spectrum, weights=None)

        self.result_signal.emit(self.retriever.result())

    def stop_retriever(self):
        self.retriever._retrieval_state.running = False


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
    prog.load_trace_in(fname='C:\\Data\\2021\\20210315\\Dataset_20210315_001\\Dataset_20210315_001.h5',
                        node_path='/Raw_datas/Scan001/Detector000/Data1D/Ch000/Data')
    prog.load_spectrum_in(fname='C:\\Users\\weber\\Desktop\\pulse.h5',
                        node_path='/Raw_datas/Detector000/Data1D/Ch000/Data')
    prog.save_data('C:\\Users\\weber\\Desktop\\pulse_analysis.h5')
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()


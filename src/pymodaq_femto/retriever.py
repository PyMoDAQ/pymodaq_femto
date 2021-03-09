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
from pymodaq_femto.simulation import Simulator
from collections import OrderedDict
from pypret import FourierTransform, Pulse, PNPS, PulsePlot, lib, MeshData
from pypret.frequencies import om2wl, wl2om, convert
import scipy
from scipy.fftpack import next_fast_len
from pymodaq.daq_utils.h5modules import H5BrowserUtil
from types import SimpleNamespace
from pyqtgraph.graphicsItems.GradientEditorItem import Gradients

config = utils.load_config()
logger = utils.set_logger(utils.get_module_name(__file__))

class DataIn(OrderedDict):
    def __init__(self, name='', source='', trace_in=None, pulse_in=None, **kwargs):
        """class subclassing from OrderedDict defining data to be processed by the retriever, either experimental or
        simulated
        Parameters
        ----------
        name: (str) data identifier
        source: (str) either "simulated" or "experimental"
        trace_in: (MeshData) MeshData object as defined in pypret and containing the trace data
        pulse_in: (Pulse) Pulse object as defined in pypret and containing fundamental spectrum (at least)
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



class RetrieverUI(QObject):
    """
    Main class initializing a DAQ_Scan module with its dashboard and scanning control panel
    """
    status_signal = pyqtSignal(str)

    params_in = [
        {'title': 'Data Info', 'name': 'data_in_info', 'type': 'group', 'children': [
            {'title': 'Trace Info', 'name': 'trace_in_info', 'type': 'group', 'children': [
                {'title': 'Wl0 (nm)', 'name': 'wl0', 'type': 'float', 'value': 0, 'readonly': True,
                 'tip': 'Central spectrum wavelength in nanometers'},
                {'title': 'FWHM (nm)', 'name': 'wl_fwhm', 'type': 'float', 'value': 0, 'readonly': True,
                 'tip': 'FWHM of the spectrum in nanometers'},
                {'title': 'Param Size', 'name': 'trace_param_size', 'type': 'int', 'value': 0, 'readonly': True},
                {'title': 'Wavelentgh Size', 'name': 'trace_wl_size', 'type': 'int', 'value': 0, 'readonly': True},
                {'title': 'Scaling (m)', 'name': 'wl_scaling', 'type': 'float', 'value': 1e-9, 'readonly': False,
                 'tip': 'Scaling to go from the Trace wavelength values to wavelength in meters'},
                {'title': 'Scaling (s)', 'name': 'param_scaling', 'type': 'float', 'value': 1e-15, 'readonly': False,
                 'tip': 'Scaling to go from the trace parameter values to delay in seconds'},
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

        self.setupUI()
        self.create_menu(self.mainwindow.menuBar())



        self.simulator = None
        self.data_in = None

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

        self.ui.dock_retriever = Dock('Retriever')
        self.dockarea.addDock(self.ui.dock_retriever, 'bottom', self.ui.dock_data_in)

        self.ui.dock_retrieved_trace = Dock('Retrieved Trace')
        self.dockarea.addDock(self.ui.dock_retrieved_trace, 'below', self.ui.dock_retriever)

        self.ui.dock_retrieved_data = Dock('Retrieved Data')
        self.dockarea.addDock(self.ui.dock_retrieved_data, 'below', self.ui.dock_retrieved_trace)

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
        self.ui.dock_retriever.addWidget(retriever_widget)
        retriever_widget.addWidget(self.viewer_live_trace.parent)
        retriever_widget.addWidget(self.viewer_error.parent)
        retriever_widget.addWidget(self.settings_retriever_tree)

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
            self.data_in = DataIn(name='data_in', source='simulated', trace_in=self.simulator.trace_exp(),
                                  pulse_in=self.simulator.pulse)
            self.display_data_in()
            self.update_spectrum_info(self.data_in['pulse_in'])
            self.update_trace_info(self.data_in['trace_in'])

    def update_spectrum_info(self, pulse):
            wl0, fwhm = utils.my_moment(pulse.wl, pulse.spectrum)
            self.settings_data_in.child('data_in_info', 'spectrum_in_info', 'wl0').setValue(wl0 * 1e9)
            self.settings_data_in.child('data_in_info', 'spectrum_in_info', 'wl_fwhm').setValue(fwhm * 1e9)
            self.settings_data_in.child('data_in_info', 'spectrum_in_info',
                                        'spectrum_size').setValue(len(pulse.spectrum))

    def update_trace_info(self, md):
        wl0, fwhm = utils.my_moment(md.axes[1], np.sum(md.data, 0))
        self.settings_data_in.child('data_in_info', 'trace_in_info', 'wl0').setValue(wl0 * 1e9)
        self.settings_data_in.child('data_in_info', 'trace_in_info', 'wl_fwhm').setValue(fwhm * 1e9)

        self.settings_data_in.child('data_in_info', 'trace_in_info', 'trace_param_size').setValue(
            len(md.axes[0]))
        self.settings_data_in.child('data_in_info', 'trace_in_info', 'trace_wl_size').setValue(len(md.axes[1]))

    def get_axes_from_trace_node(self, fname, node_path):
        h5file = self.h5browse.open_file(fname)
        data, axes, nav_axes, is_spread = self.h5browse.get_h5_data(node_path)
        self.h5browse.close_file()
        return axes['x_axis'], axes['nav_00']

    def load_trace_in(self):
        try:
            data, fname, node_path = browse_data(ret_all=True, message='Select the node corresponding to the'
                                                                       'Charaterization Trace')
            if fname is not None:
                wl, parameter_axis = self.get_axes_from_trace_node(fname, node_path)
                scaling_parameter = self.settings_data_in.child('data_in_info',
                                                                'trace_in_info', 'param_scaling').value()
                scaling_wl = self.settings_data_in.child('data_in_info', 'trace_in_info', 'wl_scaling').value()

                trace_in = MeshData(data, parameter_axis['data'] * scaling_parameter, wl['data'] * scaling_wl,
                                         labels=[parameter_axis['label'], wl['label']],
                                         units=['s', 'm'])
                self.update_trace_info(trace_in)

                dt = np.mean(np.diff(parameter_axis['data'])) * scaling_parameter
                wl0 = self.settings_data_in.child('data_in_info', 'trace_in_info', 'wl0').value() * 1e-9
                ft = FourierTransform(len(parameter_axis['data']), dt, w0=wl2om(-wl0 - 300e-9))
                self.data_in = DataIn(name='data_in', source='experimental', trace_in=trace_in,
                                      pulse_in=Pulse(ft, self.settings_data_in.child('data_in_info',
                                                                                       'trace_in_info',
                                                                                       'wl0').value() * scaling_wl))

                self.display_trace_in()

        except Exception as e:
            logger.exception(str(e))


    def load_spectrum_in(self):
        data, fname, node_path = browse_data(ret_all=True, message='Select the node corresponding to the'
                                                                   'Fundamental Spectrum')
        if fname is not None:
            h5file = self.h5browse.open_file(fname)
            data, axes, nav_axes, is_spread = self.h5browse.get_h5_data(node_path)
            self.h5browse.close_file()

            self.data_in['pulse_in'] = pulse_from_spectrum(axes['x_axis']['data'], data, pulse=self.data_in['pulse_in'])
            self.update_spectrum_info(self.data_in['pulse_in'])
            self.display_pulse_in()

    def display_trace_in(self):
        self.viewer_trace_in.setImage(self.data_in['trace_in'].data)
        self.viewer_trace_in.x_axis = utils.Axis(data=(self.data_in['trace_in'].axes[1]),
                                                 units=self.data_in['trace_in'].units[1],
                                                 label=self.data_in['trace_in'].labels[1])
        self.viewer_trace_in.y_axis = utils.Axis(data=(self.data_in['trace_in'].axes[0]),
                                                 units=self.data_in['trace_in'].units[0],
                                                 label=self.data_in['trace_in'].labels[0])
    def display_pulse_in(self):
        self.viewer_spectrum_in.show_data([utils.normalize(lib.abs2(self.data_in['pulse_in'].spectrum))],
                                      x_axis=utils.Axis(data=self.data_in['pulse_in'].wl, units='m',
                                                        label='Wavelength'),
                                          labels=['Spectrum'])

    def display_data_in(self):
        self.display_trace_in()
        self.display_pulse_in()



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

    prog = RetrieverUI(dashboard=None, dockarea=area)
    win.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

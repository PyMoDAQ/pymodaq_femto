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
from pymodaq.daq_utils.h5modules import H5Backend
config = utils.load_config()
logger = utils.set_logger(utils.get_module_name(__file__))


class DataIn(OrderedDict):
    def __init__(self, name='', source='', pnps=None, **kwargs):
        """class subclassing from OrderedDict defining data to be processed by the retriever, either experimental or
        simulated
        Parameters
        ----------
        name: (str) data identifier
        source: (str) either "simulated" or "experimental"
        pnps: (Pnps) Pnps object as defined in pypret and containing data for retrieval
        """
        if not isinstance(name, str):
            raise TypeError('name for the DataIn class should be a string')
        self['name'] = name
        if not isinstance(source, str):
            raise TypeError('source for the DataIn class should be a string')
        elif not ('simulated' in source or 'experimental' in source):
            raise ValueError('Invalid "source" for the DataIn class')
        self['source'] = source

        self['pnps'] = pnps

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
    pulse.spectrum = scipy.interpolate.interp1d(w - pulse.w0, spectrum,
                                                bounds_error=False,
                                                fill_value=0.0)(pulse.w)
    return pulse



class RetrieverUI(QObject):
    """
    Main class initializing a DAQ_Scan module with its dashboard and scanning control panel
    """
    status_signal = pyqtSignal(str)

    params = []

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

        self.h5backend = H5Backend()

        self.dockarea = dockarea
        self.dashboard = dashboard
        self.mainwindow = self.dockarea.parent()

        self.setupUI()
        self.create_menu(self.mainwindow.menuBar())

        self.simulator = None

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
        self.dockarea.addDock(self.ui.dock_retriever, 'below', self.ui.dock_data_in)

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
        self.viewer_spectrum_in = Viewer1D()
        self.settings_data_in_tree = ParameterTree()
        self.settings_data_in_tree.setMinimumWidth(300)

        data_in_splitter.addWidget(self.viewer_trace_in.parent)
        data_in_splitter.addWidget(self.viewer_spectrum_in.parent)
        data_in_splitter.addWidget(self.settings_data_in_tree)
        self.ui.dock_data_in.addWidget(self.data_in_toolbar)
        self.ui.dock_data_in.addWidget(data_in_splitter)

        params = [

            ]
        self.settings_data_in = Parameter.create(name='dataIN_settings', type='group', children=params)
        self.settings_data_in_tree.setParameters(self.settings_data_in, showTop=False)
        self.settings_data_in.sigTreeStateChanged.connect(self.data_in_settings_changed)

        # #################################################
        # setup retriever dock
        retriever_widget = QtWidgets.QSplitter()
        self.viewer_live_trace = Viewer2D()
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
            self.data_in = DataIn(name='data_in', source='simulated', pnps=self.simulator.pnps)
            self.display_data_in()

    def get_axes_from_trace_node(self, fname, node_path):
        h5file = self.h5backend.open_file(fname)
        self.h5backend.get_children(node_path.parent)


    def load_trace_in(self):
        try:
            data, fname, node_path = browse_data(ret_all=True, message='Select the node corresponding to the'
                                                                       'Charaterization Trace')
            M, N = data.shape
            
            dparameter = self.get_axes_from_trace_node(fname, node_path)
            trace = MeshData(data)
            self.data_in = DataIn(name='data_in', source='experimental')

        except Exception as e:
            logger.exception(str(e))


    def load_spectrum_in(self):
        data, fname, node_path = browse_data(ret_all=True, message='Select the node corresponding to the'
                                                                   'Charaterization Trace')

    def display_data_in(self):
        max_pnps = np.max(self.data_in['pnps'].Tmn)

        self.viewer_trace_in.setImage(utils.normalize(self.data_in['pnps'].trace.data/max_pnps))
        self.viewer_trace_in.x_axis = utils.Axis(data=(self.data_in['pnps'].w+self.data_in['pnps'].w0),
                                                       units='Hz', label='Frequency')
        self.viewer_trace_in.y_axis = utils.Axis(data=self.data_in['pnps'].parameter,
                                                 units=self.data_in['pnps'].parameter_unit,
                                                 label=utils.capitalize(self.data_in['pnps'].parameter_name))
        self.viewer_spectrum_in.show_data([utils.normalize(lib.abs2(self.data_in['pnps'].spectrum))],
                                          x_axis=utils.Axis(
                                              data=convert(self.data_in['pnps'].w + self.data_in['pnps'].w0,
                                                           "om", "wl"),
                                              units='nm', label='Wavelength'),
                                          labels=['Spectrum'])


    def data_in_settings_changed(self, param, changes):
        for param, change, data in changes:
            path = self.settings.childPath(param)
            if path is not None:
                childName = '.'.join(path)
            else:
                childName = param.name()
            if change == 'childAdded':
                pass

            elif change == 'value':
                if param.name() == 'scan_average':
                    self.show_average_dock(param.value() > 1)

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

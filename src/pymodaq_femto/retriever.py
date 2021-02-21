import sys
import subprocess

import logging
from pathlib import Path

from PyQt5 import QtGui, QtWidgets, QtCore
from PyQt5.QtCore import Qt, QObject, pyqtSlot, QThread, pyqtSignal, QLocale
from PyQt5.QtGui import QIcon, QPixmap
from pyqtgraph.dockarea import Dock
from pyqtgraph.parametertree import Parameter, ParameterTree
from pymodaq.daq_utils import daq_utils as utils
from pymodaq.daq_utils import gui_utils as gutils
from pymodaq.daq_utils.plotting.viewer2D.viewer2D_main import Viewer2D
from pymodaq.daq_utils.plotting.viewer1D.viewer1D_main import Viewer1D
from pymodaq.daq_utils.plotting.viewer0D.viewer0D_main import Viewer0D

config = utils.load_config()
logger = utils.set_logger(utils.get_module_name(__file__))

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


        self.dockarea = dockarea
        self.dashboard = dashboard
        self.mainwindow = self.dockarea.parent()

        self.setupUI()
        self.create_menu(self.mainwindow.menuBar())

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
        self.data_in_menu.addAction(self.load_trace_in)
        self.data_in_menu.addAction(self.gen_trace_in)



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
        self.load_trace_in = gutils.QAction(QIcon(QPixmap(":/icons/Icon_Library/Open.png")), 'Load Experimental Trace')
        self.gen_trace_in = gutils.QAction(QIcon(QPixmap(":/icons/Icon_Library/ini.png")),
                                           'Simulate Experimental Trace')

        self.data_in_toolbar.addAction(self.load_trace_in)
        self.data_in_toolbar.addAction(self.gen_trace_in)
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
            {'title': 'Trace Simulation:', 'name': 'trace_simulation', 'type': 'group', 'children': [
                {'title': 'Pulse Source:', 'name': 'pulse_source', 'type': 'list', 'values': ['Simulated', 'From File'],},
                {'title': 'Show Pulse:', 'name': 'show_pulse', 'type': 'bool_push', 'value': False, },
                {'title': 'Show trace:', 'name': 'show_trace', 'type': 'bool_push', 'value': False, },
                {'title': 'Pulse Settings:', 'name': 'pulse_settings', 'type': 'group', 'children': [
                    {'title': 'FWHM (fs):', 'name': 'fwhm_time', 'type': 'float', 'value': 5,
                     'tip': 'Fourier Limited Pulse duration in femtoseconds'},
                    {'title': 'GDD (fs2):', 'name': 'GDD_time', 'type': 'float', 'value': 245,
                     'tip': 'Group Delay Dispersion in femtosecond square'},
                    {'title': 'TOD (fs3):', 'name': 'TOD_time', 'type': 'float', 'value': 100,
                     'tip': 'Third Order Dispersion in femtosecond cube'},
                    {'title': 'Data File:', 'name': 'data_file_path', 'type': 'browsepath', 'filetype': True, 'visible': False,
                     'value': str(Path(__file__).parent.parent.parent.joinpath('data/spectral_data.csv')),
                     'tip': 'Path to a CSV file containing in columns: wavelength(nm), Normalized Sprectral Intensity and phase'
                            ' in radians'},
                    ]},
                {'title': 'Algorithm Options:', 'name': 'algo', 'type': 'group', 'children': [
                    {'title': 'Method:', 'name': 'method', 'type': 'list', 'values': ['frog', 'tdp', 'dscan', 'miips', 'ifrog'],
                     'tip': 'Characterization Method'},
                    {'title': 'NL process:', 'name': 'nlprocess', 'type': 'list', 'values': ['shg', 'thg', 'sd', 'pg', 'tg'],
                     'tip': 'Non Linear process used in the experiment'},
                ]},
                {'title': 'Grid settings:', 'name': 'grid_settings', 'type': 'group', 'children': [
                    {'title': 'Central Wavelength (nm):', 'name': 'wl0', 'type': 'float', 'value': 750,
                     'tip': 'Central Wavelength of the Pulse spectrum and frequency grid'},
                    {'title': 'Npoints:', 'name': 'npoints', 'type': 'list', 'values': [2**n for n in range(8,16)], 'value':512,
                     'tip': 'Number of points for the temporal and Fourier Transform Grid'},
                    {'title': 'Time resolution (fs):', 'name': 'time_resolution', 'type': 'float', 'value': 0.5,
                     'tip': 'Time spacing between 2 points in the time grid'},
                    ]},
                ]}
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

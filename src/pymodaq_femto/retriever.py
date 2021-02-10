import sys
import datetime
import subprocess
import pickle
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

        retriever_widget = QtWidgets.QSplitter()
        self.viewer_live_trace = Viewer2D()
        self.viewer_error = Viewer0D()
        self.settings_retriever_tree = ParameterTree()
        self.settings_retriever_tree.setMinimumWidth(300)
        self.ui.dock_retriever.addWidget(retriever_widget)
        retriever_widget.addWidget(self.viewer_live_trace.parent)
        retriever_widget.addWidget(self.viewer_error.parent)
        retriever_widget.addWidget(self.settings_retriever_tree)


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

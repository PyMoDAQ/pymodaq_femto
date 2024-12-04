from qtpy.QtCore import QObject
from qtpy import QtWidgets
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from pymodaq_femto.simulator import Simulator
from pyqtgraph.parametertree import ParameterTree

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

class SimulatorGUI(QObject):

    def __init__(self, parent=None):
        super().__init__()
        if parent is None:
            parent = QtWidgets.QWidget()
        self.parent = parent
        self.simulator = Simulator()

        self.setupUI()

    def setupUI(self):
        prog = Simulator()
        tree = ParameterTree()
        tree.setParameters(prog.settings, showTop=False)
        tree.setMaximumWidth(300)

        mplotlib_widget = QtWidgets.QWidget()
        sc = MplCanvas(mplotlib_widget, width=5, height=4, dpi=100)

        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        toolbar = NavigationToolbar(sc, mplotlib_widget)

        self.parent.setLayout(QtWidgets.QHBoxLayout())
        self.parent.layout().addWidget(tree)

        mplotlib_widget.setLayout(QtWidgets.QVBoxLayout())
        mplotlib_widget.layout().addWidget(toolbar)
        mplotlib_widget.layout().addWidget(sc)
        self.parent.layout().addWidget(mplotlib_widget)



def main():
    import sys

    app = QtWidgets.QApplication(sys.argv)
    win = QtWidgets.QWidget()
    win.setWindowTitle('PyMoDAQ Femto Simulator')
    prog = SimulatorGUI(win)

    win.show()

    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

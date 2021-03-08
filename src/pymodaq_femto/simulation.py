from PyQt5.QtCore import QObject
from PyQt5 import QtWidgets
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
from matplotlib.ticker import EngFormatter
from pathlib import Path
from pyqtgraph.parametertree import Parameter, ParameterTree
from pymodaq.daq_utils.parameter import pymodaq_ptypes
from pypret.frequencies import om2wl, wl2om, convert
from pypret import FourierTransform, Pulse, PNPS, PulsePlot, lib, MeshData
from pypret.graphics import plot_complex, plot_meshdata
from scipy.interpolate import interp2d
import numpy as np
from pymodaq.daq_utils.daq_utils import gauss1D, my_moment
from pypret.pnps import _PNPS_CLASSES

methods = list(_PNPS_CLASSES.keys())
methods.pop(methods.index('dscan'))
methods.sort()
nlprocesses = list(_PNPS_CLASSES[methods[0]].keys())

def normalize(x):
    x = x - np.min(x)
    x = x / np.max(x)
    return x

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(fig)


class PulsePlot:

    def __init__(self, pulse, fig=None, plot=True, **kwargs):
        self.pulse = pulse
        self.fig = fig
        if plot:
            self.plot(**kwargs)

    def plot(self, xaxis='wavelength', yaxis='intensity', limit=True,
             oversampling=False, phase_blanking=False,
             phase_blanking_threshold=1e-3, show=True):
        pulse = self.pulse

        if self.fig is None:
            fig, axs = plt.subplots(1, 2)
        else:
            axs = self.fig.subplots(nrows=1, ncols=2)

        fig = self.fig
        ax1, ax2 = axs.flat
        ax12 = ax1.twinx()
        ax22 = ax2.twinx()

        if oversampling:
            t = np.linspace(pulse.t[0], pulse.t[-1], pulse.N * oversampling)
            field = pulse.field_at(t)
        else:
            t = pulse.t
            field = pulse.field

        # time domain
        li11, li12, tamp, tpha = plot_complex(t, field, ax1, ax12, yaxis=yaxis,
                                              phase_blanking=phase_blanking, limit=limit,
                                              phase_blanking_threshold=phase_blanking_threshold)
        fx = EngFormatter(unit="s")
        ax1.xaxis.set_major_formatter(fx)
        ax1.set_title("time domain")
        ax1.set_xlabel("time")
        ax1.set_ylabel(yaxis)
        ax12.set_ylabel("phase (rad)")

        # frequency domain
        if oversampling:
            w = np.linspace(pulse.w[0], pulse.w[-1], pulse.N * oversampling)
            spectrum = pulse.spectrum_at(w)
        else:
            w = pulse.w
            spectrum = pulse.spectrum

        if xaxis == "wavelength":
            w = convert(w + pulse.w0, "om", "wl")
            unit = "m"
            label = "wavelength"
        elif xaxis == "frequency":
            w = w
            unit = " rad Hz"
            label = "frequency"

        li21, li22, samp, spha = plot_complex(w, spectrum, ax2, ax22, yaxis=yaxis,
                                              phase_blanking=phase_blanking, limit=limit,
                                              phase_blanking_threshold=phase_blanking_threshold)
        fx = EngFormatter(unit=unit)
        ax2.xaxis.set_major_formatter(fx)
        ax2.set_title("frequency domain")
        ax2.set_xlabel(label)
        ax2.set_ylabel(yaxis)
        ax22.set_ylabel("phase (rad)")

        self.fig = fig
        self.ax1, self.ax2 = ax1, ax2
        self.ax12, self.ax22 = ax12, ax22
        self.li11, self.li12, self.li21, self.li22 = li11, li12, li21, li22
        self.tamp, self.tpha = tamp, tpha
        self.samp, self.spha = samp, spha

        if show:
            fig.tight_layout()
        #     plt.show()


def plot_meshdata(ax, md, cmap="nipy_spectral", width_factor=6):
    x0, dx = my_moment(md.axes[1], np.sum(md.data, 0))
    y0, dy = my_moment(md.axes[0], np.sum(md.data, 1))
    x, y = lib.edges(md.axes[1]), lib.edges(md.axes[0])
    indx1 = np.argwhere(x >= x0-width_factor*dx)[0][0]
    indx2 = np.argwhere(x >= x0 + width_factor*dx)[0][0]
    indy1 = np.argwhere(y <= y0 + width_factor*dy)[0][0]
    indy2 = np.argwhere(y <= y0 - width_factor*dy)[0][0]

    im = ax.pcolormesh(x[indx1:indx2+1], y[indy1:indy2+1], normalize(md.data[indy1:indy2+1, indx1:indx2+1]), cmap=cmap)
    ax.set_xlabel(md.labels[1])
    ax.set_ylabel(md.labels[0])

    fx = EngFormatter(unit=md.units[1])
    ax.xaxis.set_major_formatter(fx)
    fy = EngFormatter(unit=md.units[0])
    ax.yaxis.set_major_formatter(fy)
    return im


class MeshDataPlot:

    def __init__(self, mesh_data, fig=None, plot=True, **kwargs):
        self.md = mesh_data
        self.fig = fig
        if plot:
            self.plot(**kwargs)

    def plot(self, show=True, **kwargs):
        md = self.md
        if self.fig is None:
            fig, ax = plt.subplots()
        else:
            fig = self.fig
            ax = self.fig.subplots(nrows=1, ncols=1)
            
        im = plot_meshdata(ax, md, "nipy_spectral", **kwargs)
        fig.colorbar(im, ax=ax)

        self.fig, self.ax = fig, ax
        self.im = im
        if show:
            fig.tight_layout()
            plt.show()

    def show(self):
        plt.show()

class Simulator(QObject):
    params = [
            {'title': 'Show Pulse', 'name': 'show_pulse', 'type': 'action', 'visible': False},
            {'title': 'Show Trace', 'name': 'show_trace', 'type': 'action', 'visible': False},
            {'title': 'Show both', 'name': 'show_plots', 'type': 'action', 'visible': False},
            {'title': 'Pulse Source:', 'name': 'pulse_source', 'type': 'list', 'values': ['Simulated', 'From File'],
             },

            {'title': 'Pulse Settings:', 'name': 'pulse_settings', 'type': 'group', 'children': [
                {'title': 'FWHM (fs):', 'name': 'fwhm_time', 'type': 'float', 'value': 5,
                 'tip': 'Fourier Limited Pulse duration in femtoseconds'},
                {'title': 'GDD (fs2):', 'name': 'GDD', 'type': 'float', 'value': 50,
                 'tip': 'Group Delay Dispersion in femtosecond square'},
                {'title': 'TOD (fs3):', 'name': 'TOD', 'type': 'float', 'value': 500,
                 'tip': 'Third Order Dispersion in femtosecond cube'},
                {'title': 'Data File:', 'name': 'data_file_path', 'type': 'browsepath', 'filetype': True,
                 'visible': False,
                 'value': str(Path(__file__).parent.parent.parent.joinpath('data/spectral_data.csv')),
                 'tip': 'Path to a CSV file containing in columns: wavelength(nm), Normalized Sprectral Intensity and phase'
                        ' in radians'},
            ]},
            {'title': 'Algorithm Options:', 'name': 'algo', 'type': 'group', 'children': [
                {'title': 'Method:', 'name': 'method', 'type': 'list',
                 'values': methods,
                 'tip': 'Characterization Method'},
                {'title': 'NL process:', 'name': 'nlprocess', 'type': 'list',
                 'values': nlprocesses,
                 'tip': 'Non Linear process used in the experiment'},
            ]},
            {'title': 'Grid settings:', 'name': 'grid_settings', 'type': 'group', 'children': [
                {'title': 'lambda0 (nm):', 'name': 'wl0', 'type': 'float', 'value': 750,
                 'tip': 'Central Wavelength of the Pulse spectrum and frequency grid'},
                {'title': 'Npoints:', 'name': 'npoints', 'type': 'list', 'values': [2 ** n for n in range(8, 16)],
                 'value': 1024,
                 'tip': 'Number of points for the temporal and Fourier Transform Grid'},
                {'title': 'Time resolution (fs):', 'name': 'time_resolution', 'type': 'float', 'value': 0.5,
                 'tip': 'Time spacing between 2 points in the time grid'},
            ]},
        ]

    def __init__(self, parent=None, show_ui=True):
        super().__init__()

        if parent is None:
            parent = QtWidgets.QWidget()

        self.parent = parent
        self.figs = []
        self.pnps = None
        self.max_pnps = 1
        self.pulse = None

        self.settings = Parameter.create(name='dataIN_settings', type='group', children=self.params)
        self.settings.sigTreeStateChanged.connect(self.settings_changed)

        if show_ui:
            self.setupUI()
            self.settings.child('show_plots').sigActivated.connect(self.show_pulse)
            self.settings.child('show_plots').show()
            self.settings.child('show_pulse').show()
            self.settings.child('show_trace').show()
            self.settings.child('show_plots').sigActivated.connect(self.show_trace)
            self.settings.child('show_trace').sigActivated.connect(self.show_trace)
            self.settings.child('show_pulse').sigActivated.connect(self.show_pulse)
        else:
            self.settings.child('show_plots').hide()
            self.settings.child('show_pulse').hide()
            self.settings.child('show_trace').hide()

        self.update_pulse()
        self.update_pnps()

    @property
    def trace_exp(self):
        """ Experimental trace on linear wavelength grid of the simulated trace

        Returns
        -------
        meshdata: (MeshData)
        """
        width_factor = 6
        w, y = lib.edges(self.pnps.trace.axes[1]), lib.edges(self.pnps.trace.axes[0])
        x0, dx = my_moment(self.pnps.trace.axes[1], np.sum(self.pnps.trace.data, 0))
        y0, dy = my_moment(self.pnps.trace.axes[0], np.sum(self.pnps.trace.data, 1))
        indx1 = np.argwhere(self.pnps.trace.axes[1] >= x0 - width_factor * dx)[0][0]
        indx2 = np.argwhere(self.pnps.trace.axes[1] >= x0 + width_factor * dx)[0][0]

        wl = self.pnps.wl[indx1:indx2]
        new_wl = np.linspace(np.min(wl), np.max(wl), 256)
        new_w = convert(new_wl, 'wl', 'om')
        new_data = interp2d(self.pnps.trace.axes[1][indx1:indx2], y, self.pnps.trace.data[indx1:indx2, :])(new_w, y)

        return MeshData(new_data, self.pnps.trace.axes[0], new_wl * 1e9,)

    @property
    def trace(self):
        return self.pnps.trace

    @property
    def parameter(self):
        return self.pnps.parameter


    def setupUI(self):
        self.settings_tree = ParameterTree()
        self.settings_tree.setParameters(self.settings, showTop=False)
        self.settings_tree.setMaximumWidth(300)

        mplotlib_widget = QtWidgets.QWidget()
        self.pulse_canvas = MplCanvas(mplotlib_widget, width=5, height=4, dpi=100)
        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        toolbar_pulse = NavigationToolbar(self.pulse_canvas, mplotlib_widget)

        self.trace_canvas = MplCanvas(mplotlib_widget, width=5, height=4, dpi=100)
        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        toolbar_trace = NavigationToolbar(self.trace_canvas, mplotlib_widget)

        self.parent.setLayout(QtWidgets.QHBoxLayout())
        self.parent.layout().addWidget(self.settings_tree)

        mplotlib_widget.setLayout(QtWidgets.QVBoxLayout())
        mplotlib_widget.layout().addWidget(toolbar_pulse)
        mplotlib_widget.layout().addWidget(self.pulse_canvas)
        mplotlib_widget.layout().addWidget(toolbar_trace)
        mplotlib_widget.layout().addWidget(self.trace_canvas)
        self.parent.layout().addWidget(mplotlib_widget)

        self.set_tight_layout(True)

    def settings_changed(self, param, changes):
        for param, change, data in changes:
            path = self.settings.childPath(param)
            if change == 'childAdded':
                pass
            elif change == 'parent':
                pass
            elif change == 'value':
                if param.name() == 'pulse_source':
                    for child in self.settings.child('pulse_settings').children():
                        if child.name() == 'data_file_path':
                            child.show(param.value() == 'From File')
                        else:
                            child.show(param.value() != 'From File')
                elif param.name() == 'method':
                    self.settings.child('algo', 'nlprocess').setLimits(list(_PNPS_CLASSES[param.value()].keys()))

    def set_tight_layout(self, tight=True):
        self.pulse_canvas.figure.set_tight_layout(tight)
        self.trace_canvas.figure.set_tight_layout(tight)

    def show_pulse(self):
        self.update_pulse()
        self.pulse_canvas.figure.clf()
        PulsePlot(self.pulse, self.pulse_canvas.figure)
        self.pulse_canvas.draw()

    def show_trace(self):
        self.update_pnps()
        self.trace_canvas.figure.clf()
        MeshDataPlot(self.pnps.trace, self.trace_canvas.figure)
        self.trace_canvas.draw()

    def update_grid(self):
        Nt = self.settings.child('grid_settings', 'npoints').value()
        dt = self.settings.child('grid_settings', 'time_resolution').value() * 1e-15
        wl0 = self.settings.child('grid_settings', 'wl0').value() * 1e-9
        self.ft = FourierTransform(Nt, dt=dt, w0=wl2om(-wl0 - 300e-9))

    def update_pnps(self):

        pulse = self.update_pulse()
        method = self.settings.child('algo', 'method').value()
        process = self.settings.child('algo', 'nlprocess').value()
        self.pnps = PNPS(pulse, method, process)
        parameter = np.linspace(self.ft.t[-1], self.ft.t[0], len(self.ft.t))
        self.pnps.calculate(pulse.spectrum, parameter)
        self.max_pnps = np.max(self.pnps.Tmn)
        return self.pnps

    def update_pulse(self):
        self.update_grid()
        wl0 = self.settings.child('grid_settings', 'wl0').value() * 1e-9
        w0 = convert(wl0, 'wl', 'om')
        pulse = Pulse(self.ft, wl0)

        if self.settings.child('pulse_source').value() == 'Simulated':
            fwhm = self.settings.child('pulse_settings', 'fwhm_time').value()
            domega = 4*np.log(2) / fwhm
            GDD = self.settings.child('pulse_settings', 'GDD').value()
            TOD = self.settings.child('pulse_settings', 'TOD').value()

            pulse.spectrum = gauss1D(pulse.w, x0=0., dx=domega * 1e15)  # x0=0 because the frequency axis is already
            # centered on wl0 (see Pulse(self.ft, wl0))
            #pulse.field = gauss1D(pulse.t, x0=0, dx=fwhm * 1e-15)
            pulse.spectrum = pulse.spectrum * np.exp(
                1j * (GDD * 1e-30) * ((pulse.w) ** 2) / 2 + 1j * (TOD * 1e-45) * (
                        (pulse.w) ** 3) / 6)

            # # recenter pulse in time domain
            # idx = np.argmax(pulse.intensity)
            # pulse.spectrum = pulse.spectrum * np.exp(-1j * pulse.t[idx] * (pulse.w - pulse.w0))

        else:
            data_path = self.settings.child('pulse_settings', 'data_file_path').value()
            data = np.genfromtxt(data_path, delimiter=',', skip_header=1)
            in_wl, in_int, in_phase = (data[:, i] for i in range(3))

            in_int = np.interp(pulse.wl, in_wl * 1e-9, np.maximum(0, in_int), left=0, right=0)
            in_phase = np.interp(pulse.wl, in_wl * 1e-9, in_phase, left=0, right=0)
            pulse.spectrum = in_int * np.exp(1j * in_phase)

        self.pulse = pulse
        return pulse




def main():
    import sys

    app = QtWidgets.QApplication(sys.argv)
    win = QtWidgets.QWidget()
    win.setWindowTitle('PyMoDAQ Femto Simulator')
    prog = Simulator(win, show_ui=True)

    win.show()

    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
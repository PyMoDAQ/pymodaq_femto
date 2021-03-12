import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import EngFormatter
from pypret import (Pulse, lib)
from pypret.frequencies import convert
from pypret.graphics import plot_complex, plot_meshdata


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


class RetrievalResultPlot:

    def __init__(self, retrieval_result, fig=None, plot=True, **kwargs):
        self.retrieval_result = retrieval_result
        self.fig = fig
        if plot:
            self.plot(**kwargs)

    def plot(self, xaxis='wavelength', yaxis='intensity', limit=True,
             oversampling=False, phase_blanking=False,
             phase_blanking_threshold=1e-3, show=True, fundamental=None,
             fundamental_wavelength=None):
        rr = self.retrieval_result
        # reconstruct a pulse from that
        pulse = Pulse(rr.pnps.ft, rr.pnps.w0, unit="om")
        if self.fig is None:
            fig = plt.figure(figsize=(30.0/2.54, 20.0/2.54))
        else:
            fig = self.fig
        # construct the figure

        gs1 = gridspec.GridSpec(2, 2, figure=fig)
        gs2 = gridspec.GridSpec(2, 6, figure=fig)
        ax1 = fig.add_subplot(gs1[0, 0])
        ax2 = fig.add_subplot(gs1[0, 1])
        ax3 = fig.add_subplot(gs2[1, :2])
        ax4 = fig.add_subplot(gs2[1, 2:4])
        ax5 = fig.add_subplot(gs2[1, 4:])
        ax12 = ax1.twinx()
        ax22 = ax2.twinx()

        # Plot in time domain
        pulse.spectrum = rr.pulse_retrieved  # the retrieved pulse
        if oversampling:
            t = np.linspace(pulse.t[0], pulse.t[-1], pulse.N * oversampling)
            field2 = pulse.field_at(t)
        else:
            t = pulse.t
            field2 = pulse.field
        field2 /= np.abs(field2).max()

        li11, li12, tamp2, tpha2 = plot_complex(t, field2, ax1, ax12, yaxis=yaxis,
                          phase_blanking=phase_blanking, limit=limit,
                          phase_blanking_threshold=phase_blanking_threshold)
        li11.set_linewidth(3.0)
        li11.set_color("#1f77b4")
        li11.set_alpha(0.6)
        li12.set_linewidth(3.0)
        li12.set_color("#ff7f0e")
        li12.set_alpha(0.6)

        fx = EngFormatter(unit="s")
        ax1.xaxis.set_major_formatter(fx)
        ax1.set_title("time domain")
        ax1.set_xlabel("time")
        ax1.set_ylabel(yaxis)
        ax12.set_ylabel("phase (rad)")
        ax1.legend([li11, li12], [yaxis, "phase"])

        # frequency domain        
        if oversampling:
            w = np.linspace(pulse.w[0], pulse.w[-1], pulse.N * oversampling)            
            spectrum2 = pulse.spectrum_at(w)
            pulse.spectrum = rr.pulse_retrieved
        else:
            w = pulse.w
            spectrum2 = rr.pulse_retrieved
        fund_w = convert(fundamental_wavelength, "wl", "om") - pulse.w0
        scale = np.abs(spectrum2).max()
        spectrum2 /= scale
        if fundamental is not None:
            fundamental /= scale*scale

        if xaxis == "wavelength":
            w = convert(w + pulse.w0, "om", "wl")
            fund_w = fundamental_wavelength
            unit = "m"
            label = "wavelength"
        elif xaxis == "frequency":
            unit = " rad Hz"
            label = "frequency"        

        # Plot in spectral domain
        li21, li22, samp2, spha2 = plot_complex(w, spectrum2, ax2, ax22, yaxis=yaxis,
                          phase_blanking=phase_blanking, limit=limit,
                          phase_blanking_threshold=phase_blanking_threshold)
        lines = [li21, li22]
        labels = ["intensity", "phase"]
        if fundamental is not None:
            li31, = ax2.plot(fund_w, fundamental, "r+", ms=4.0, mew=1.0, zorder=0)
            lines.append(li31)
            labels.append("measurement")
        li21.set_linewidth(3.0)
        li21.set_color("#1f77b4")
        li21.set_alpha(0.6)
        li22.set_linewidth(3.0)
        li22.set_color("#ff7f0e")
        li22.set_alpha(0.6)

        fx = EngFormatter(unit=unit)
        ax2.xaxis.set_major_formatter(fx)
        ax2.set_title("frequency domain")
        ax2.set_xlabel(label)
        ax2.set_ylabel(yaxis)
        ax22.set_ylabel("phase (rad)")
        ax2.legend(lines, labels)

        axes = [ax3, ax4, ax5]
        sc = 1.0 / rr.trace_input.max()
        traces = [rr.trace_input * sc, rr.trace_retrieved * sc,
                  (rr.trace_input - rr.trace_retrieved) * rr.weights * sc]
        titles = ["measured", "retrieved", "difference"]
        if np.any(rr.weights != 1.0):
            titles[-1] = "weighted difference"
        cmaps = ["nipy_spectral", "nipy_spectral", "RdBu"]
        md = rr.measurement
        for ax, trace, title, cmap in zip(axes, traces, titles, cmaps):
            x, y = lib.edges(rr.pnps.process_w), lib.edges(rr.parameter)
            im = ax.pcolormesh(x, y, trace, cmap=cmap)
            fig.colorbar(im, ax=ax)
            ax.set_xlabel(md.labels[1])
            ax.set_ylabel(md.labels[0])
            fx = EngFormatter(unit=md.units[1])
            ax.xaxis.set_major_formatter(fx)
            fy = EngFormatter(unit=md.units[0])
            ax.yaxis.set_major_formatter(fy)
            ax.set_title(title)
            ax.set_xlim(lib.limit(md.axes[1], md.marginals(axes=1)))

        ax1.grid()
        ax2.grid()
        self.fig = fig
        self.ax1, self.ax2 = ax1, ax2
        self.ax12, self.ax22 = ax12, ax22
        self.li11, self.li12, self.li21, self.li22 = li11, li12, li21, li22
        self.ax3, self.ax4, self.ax5 = ax3, ax4, ax5

        if show:
            #gs.tight_layout(fig)
            gs1.update(left=0.05, right=0.95, top=0.9, bottom=0.1,
                      hspace=0.25, wspace=0.3)
            gs2.update(left=0.1, right=0.95, top=0.9, bottom=0.1,
                      hspace=0.5, wspace=1.0)
            plt.show()
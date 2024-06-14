import matplotlib

matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import EngFormatter
from pypret import Pulse, lib
from pypret.frequencies import convert, om2wl, wl2om
from pypret.graphics import plot_complex, plot_meshdata
import math

def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper


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

    def plot(
        self,
        xaxis="wavelength",
        yaxis="intensity",
        limit=True,
        oversampling=False,
        phase_blanking=False,
        phase_blanking_threshold=1e-3,
        show=True,
    ):
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
        li11, li12, tamp, tpha = plot_complex(
            t,
            field,
            ax1,
            ax12,
            yaxis=yaxis,
            phase_blanking=phase_blanking,
            limit=limit,
            phase_blanking_threshold=phase_blanking_threshold,
        )
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

        li21, li22, samp, spha = plot_complex(
            w,
            spectrum,
            ax2,
            ax22,
            yaxis=yaxis,
            phase_blanking=phase_blanking,
            limit=limit,
            phase_blanking_threshold=phase_blanking_threshold,
        )
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


class PulsePropagationPlot(PulsePlot):
    def __init__(self, pulse, polynomial, fwhm=None, fig=None, plot=True, **kwargs):
        self.polynomial = polynomial
        self.pulse = pulse
        self.fig = fig
        self.fwhm = fwhm
        if plot:
            self.plot(**kwargs)

    def plot(
        self,
        xaxis="wavelength",
        yaxis="intensity",
        limit=True,
        oversampling=False,
        phase_blanking=False,
        phase_blanking_threshold=1e-3,
        show=True,
    ):

        super().plot(
            xaxis=xaxis,
            yaxis=yaxis,
            limit=limit,
            oversampling=oversampling,
            phase_blanking=phase_blanking,
            phase_blanking_threshold=phase_blanking_threshold,
            show=show,
        )

        # #Add fwhm plot
        if oversampling:
            self.t = np.linspace(
                self.pulse.t[0], self.pulse.t[-1], self.pulse.N * oversampling
            )
        else:
            self.t = self.fundamental.t
        self.intensity_fwhm = (
            lib.gaussian(
                self.t,
                self.t[np.argmax(self.tamp)],
                sigma=0.5 * (self.fwhm * 1e-15) / np.sqrt(2 * np.log(2.0)),
            )
            * self.tamp.max()
        )
        self.ax1.plot(self.t, self.intensity_fwhm, "r--", alpha=0.5)

        self.ax1.set_zorder(1)  # default zorder is 0 for ax1 and ax2
        self.ax1.patch.set_visible(False)  # prevents ax1 from hiding ax2

        # Add polynomial fit of phase
        if xaxis == "wavelength":
            w = convert(self.pulse.w + self.pulse.w0, "om", "wl")
        elif xaxis == "frequency":
            w = self.pulse.w

        phase_fit = np.poly1d(self.polynomial)(self.pulse.w)
        if yaxis == "intensity":
            amp = lib.abs2(self.pulse.spectrum)
        elif yaxis == "amplitude":
            amp = np.abs(self.pulse.spectrum)

        phase_fit -= lib.mean(phase_fit, amp * amp)
        self.ax22.plot(w, phase_fit, "--")

        if show:
            self.fig.tight_layout()


class MeshDataPlot:
    def __init__(self, mesh_data, fig=None, plot=True, limit=False, **kwargs):
        self.md = mesh_data
        self.fig = fig
        self.limit = limit
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
        if self.limit:
            ax.set_xlim(lib.limit(md.axes[1], md.marginals(axes=1)))

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

    def plot(
        self,
        xaxis="wavelength",
        yaxis="intensity",
        limit=True,
        oversampling=False,
        phase_blanking=False,
        phase_blanking_threshold=1e-3,
        show=True,
        fundamental=None,
        fundamental_wavelength=None,
        compare_fundamental=True
    ):
        rr = self.retrieval_result
        # reconstruct a pulse from that
        pulse = Pulse(rr.pnps.ft, rr.pnps.w0, unit="om")
        pulse.spectrum = rr.pulse_retrieved  # the retrieved pulse

        if rr.pnps.method =='dscan':
            # Look for optimal insertion
            propagated_pulse = Pulse(rr.pnps.ft, rr.pnps.w0, unit="om")

            precision = pulse.dt/10
            item = rr.pnps.material

            w1, w2 = sorted(wl2om(np.array(item._range)))
            w = propagated_pulse.w + propagated_pulse.w0
            valid = (w >= w1) & (w <= w2)

            w = w[valid]
            k = item.k(w, unit="om")
            k0 = item.k(propagated_pulse.w0, unit="om")
            k1 = item.k(propagated_pulse.w0 + propagated_pulse.ft.dw, unit="om")
            dk = (k1 - k0) / propagated_pulse.ft.dw

            # Add material dispersion without 0th and 1st Taylor orders (they don't change the pulse)
            kfull = np.zeros_like(propagated_pulse.w)
            kfull[valid] = k - k0 - dk * propagated_pulse.w[valid]

            shortest_fwhm = 1e50
            optimal_insertion = 0
            optimal_pulse = Pulse(rr.pnps.ft, rr.pnps.w0, unit="om")

            for insertion in rr.pnps.parameter:
                propagated_pulse.spectrum = pulse.spectrum
                propagated_pulse.spectrum *= np.exp(1j * kfull * insertion)
                fwhm = propagated_pulse.fwhm(precision)

                if fwhm < shortest_fwhm:
                    shortest_fwhm = fwhm
                    optimal_insertion = insertion
                    optimal_pulse.spectrum = propagated_pulse.spectrum

            pulse.spectrum = optimal_pulse.spectrum    # Use shortest pulse for display

        if self.fig is None:
            fig = plt.figure(figsize=(30.0 / 2.54, 20.0 / 2.54))
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
        if oversampling:
            t = np.linspace(pulse.t[0], pulse.t[-1], pulse.N * oversampling)
            field2 = pulse.field_at(t)
        else:
            t = pulse.t
            field2 = pulse.field
        field2 /= np.abs(field2).max()

        li11, li12, tamp2, tpha2 = plot_complex(
            t,
            field2,
            ax1,
            ax12,
            yaxis=yaxis,
            phase_blanking=phase_blanking,
            limit=limit,
            phase_blanking_threshold=phase_blanking_threshold,
        )
        li11.set_linewidth(3.0)
        li11.set_color("#1f77b4")
        li11.set_alpha(0.6)
        li12.set_linewidth(3.0)
        li12.set_color("#ff7f0e")
        li12.set_alpha(0.6)

        fx = EngFormatter(unit="s")
        ax1.xaxis.set_major_formatter(fx)
        if rr.pnps.method == "dscan":
            ax1.set_title("Temporal intensity at insertion = " + str(truncate(optimal_insertion*1000,4))
                          + 'mm [FWHM = ' + str(truncate(shortest_fwhm * 1e15, 4)) + " fs].")
        else:
            ax1.set_title("Temporal intensity")

        ax1.set_xlabel("time")
        ax1.set_ylabel(yaxis)
        ax12.set_ylabel("phase (rad)")
        ax1.legend([li11, li12], [yaxis, "phase"])

        # frequency domain
        if oversampling:
            w = np.linspace(pulse.w[0], pulse.w[-1], pulse.N * oversampling)
            spectrum2 = pulse.spectrum_at(w)
        else:
            w = pulse.w
            spectrum2 = rr.pulse_retrieved
        fund_w = convert(fundamental_wavelength, "wl", "om") - pulse.w0
        scale = np.abs(spectrum2).max()
        spectrum2 /= scale
        if fundamental is not None:
            fundamental /= np.abs(fundamental).max()

        if xaxis == "wavelength":
            w = convert(w + pulse.w0, "om", "wl")
            fund_w = fundamental_wavelength
            unit = "m"
            label = "wavelength"
        elif xaxis == "frequency":
            unit = " rad Hz"
            label = "frequency"

        # Plot in spectral domain
        li21, li22, samp2, spha2 = plot_complex(
            w,
            spectrum2,
            ax2,
            ax22,
            yaxis=yaxis,
            phase_blanking=phase_blanking,
            limit=limit,
            phase_blanking_threshold=phase_blanking_threshold,
        )
        lines = [li21, li22]
        labels = ["intensity", "phase"]
        if fundamental is not None and compare_fundamental:
            (li31,) = ax2.plot(fund_w, fundamental, "r+", ms=4.0, mew=1.0, zorder=0)
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
        traces = [
            rr.trace_input * sc,
            rr.trace_retrieved * sc,
            (rr.trace_input - rr.trace_retrieved) * rr.weights * sc,
        ]
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
            ax.set_xlim(
                lib.limit(md.axes[1], md.marginals(axes=1), threshold=1e-2, padding=0.1)
            )

        # if rr.pnps.method == 'dscan':
        #     wl_lim = ax3.get_xlim()
        #     d = (wl_lim[1]-wl_lim[0])/5
        #     ax3.annotate("", xy=(-0.5*d, optimal_insertion), xytext=(-1.5*d, optimal_insertion), xycoords='axes fraction',
        #     arrowprops=dict(arrowstyle="->"))

        ax1.grid()
        ax2.grid()
        self.fig = fig
        self.ax1, self.ax2 = ax1, ax2
        self.ax12, self.ax22 = ax12, ax22
        self.li11, self.li12, self.li21, self.li22 = li11, li12, li21, li22
        self.ax3, self.ax4, self.ax5 = ax3, ax4, ax5

        if show:
            # gs.tight_layout(fig)
            gs1.update(
                left=0.05, right=0.95, top=0.9, bottom=0.1, hspace=0.25, wspace=0.3
            )
            gs2.update(
                left=0.1, right=0.95, top=0.9, bottom=0.1, hspace=0.5, wspace=1.0
            )
            plt.show()

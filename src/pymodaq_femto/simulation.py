from PyQt5.QtCore import QObject
from PyQt5 import QtWidgets

from pathlib import Path
from pyqtgraph.parametertree import Parameter, ParameterTree
from pymodaq.daq_utils.parameter import pymodaq_ptypes
from pypret.frequencies import om2wl, wl2om, convert
from pypret import FourierTransform, Pulse, PNPS, lib, MeshData

import numpy as np
from pymodaq.daq_utils.daq_utils import gauss1D, my_moment, l2w, linspace_step, Axis, normalize
from pymodaq.daq_utils.array_manipulation import linspace_this_image, crop_vector_to_axis, crop_array_to_axis,\
    linspace_this_vect
from pypret.material import BK7
from pymodaq_femto.materials import FS_extended
from pymodaq_femto.graphics import MplCanvas, NavigationToolbar, MeshDataPlot, PulsePlot
from collections import OrderedDict
from pymodaq_femto import _PNPS_CLASSES



methods_tmp = list(_PNPS_CLASSES.keys())
methods_tmp.sort()
methods = ['frog']
methods.extend(methods_tmp)
nlprocesses = list(_PNPS_CLASSES[methods[0]].keys())
materials = OrderedDict(FS=FS_extended, BK7=BK7)



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
                {'title': 'Shaping type:', 'name': 'shaping_type', 'type': 'list', 'values': ['Taylor', 'Gaussian'],
                 },
                {'title': 'Npulses:', 'name': 'npulses', 'type': 'int', 'value': 1,
                 'tip': 'Number of pulse in a sequence'},
                {'title': 'Pulses separation:', 'name': 'delay_pulses', 'type': 'float', 'value': 100,
                 'tip': 'Delay between pulses in femtosecond', 'visible': False},
                {'title': 'Taylor Phase:', 'name': 'taylor_phase', 'type': 'group', 'children': [
                    {'title': 'Delay (fs):', 'name': 'GD', 'type': 'float', 'value': 0,
                     'tip': 'Group Delay in femtosecond'},
                    {'title': 'GDD (fs2):', 'name': 'GDD', 'type': 'float', 'value': 50,
                     'tip': 'Group Delay Dispersion in femtosecond square'},
                    {'title': 'TOD (fs3):', 'name': 'TOD', 'type': 'float', 'value': 500,
                     'tip': 'Third Order Dispersion in femtosecond cube'},
                ]},
                {'title': 'Gaussian Phase:', 'name': 'gaussian_phase', 'type': 'group', 'visible': False, 'children': [
                    {'title': 'Amplitude (rad):', 'name': 'gauss_amp', 'type': 'float', 'value': 6,
                     'tip': 'Amplitude of the gaussian phase in radian'},
                    {'title': 'dt (fs):', 'name': 'dtime', 'type': 'float', 'value': 10,
                     'tip': 'FWHM (in fs) of the gaussian temporal phase'},
                ]},

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
                {'title': 'Alpha (rad):', 'name': 'alpha', 'type': 'float', 'value': 1,
                 'tip': 'amplitude of the phase pattern (in rad)', 'visible': False},
                {'title': 'Gamma (Hz):', 'name': 'gamma', 'type': 'float', 'value': 10,
                 'tip': 'frequency of the phase pattern (in Hz)', 'visible': False},
                {'title': 'Material:', 'name': 'material', 'type': 'list',
                 'values': list(materials.keys()), 'visible': False,
                 'tip': 'Material used for the Dscan measurement'},
                {'title': 'Dscan Parameter Scan:', 'name': 'dscan_parameter', 'type': 'group', 'visible': False,
                 'children': [
                     {'title': 'Insertion min (mm):', 'name': 'min', 'type': 'float', 'value': -10.,
                      'tip': 'Minimum of the scanned parameter in mm'},
                     {'title': 'Insertion max (mm):', 'name': 'max', 'type': 'float', 'value': 10.,
                      'tip': 'Minimum of the scanned parameter in mm'},
                     {'title': 'Insertion step (mm):', 'name': 'step', 'type': 'float', 'value': 0.025,
                      'tip': 'Step size of the scanned parameter in mm'},
                 ]},
                {'title': 'MIIPS Parameter Scan:', 'name': 'miips_parameter', 'type': 'group', 'visible': False,
                 'children': [
                     {'title': 'Phase min (rad):', 'name': 'min', 'type': 'float', 'value': 0,
                      'tip': 'Minimum of the scanned parameter in radians'},
                     {'title': 'Phase max (rad):', 'name': 'max', 'type': 'float', 'value': 2 * np.pi,
                      'tip': 'Minimum of the scanned parameter in radian'},
                     {'title': 'Phase setp (rad):', 'name': 'step', 'type': 'float', 'value': 2 * np.pi / 100,
                      'tip': 'Step size of the scanned parameter in radians'},
                 ]},
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
            {'title': 'Plot settings:', 'name': 'plot_settings', 'type': 'group', 'children': [
                {'title': 'Units:', 'name': 'units', 'type': 'list', 'values': ['nm', 'Hz'],
                 'tip': 'Plot ad a function of the wavelength (in nm) or as a function of the angular frequency (in Hz)'},
                {'title': 'Autolimits?:', 'name': 'autolimits', 'type': 'bool', 'value': True,
                 'tip': 'Restrict the data plot to limits given from marginals and threshold'},
                {'title': 'Set Limits?:', 'name': 'setlimits', 'type': 'bool', 'value': False,
                 'tip': 'Restrict the data plot to limits given from marginals and threshold'},
                {'title': 'Autolimits Threshold:', 'name': 'autolim_thresh', 'type': 'float', 'value': 1e-2,
                 'tip': 'Threshold for the determination of the plotting limits'},
                {'title': 'Limit min:', 'name': 'limit_min', 'type': 'float', 'value': 500,
                 'tip': 'Min  value of the frequency axis for plotting (Hz or nm)', 'visible': False},
                {'title': 'Limit max:', 'name': 'limit_max', 'type': 'float', 'value': 1100,
                 'tip': 'Max  value of the frequency axis for plotting (Hz or nm)', 'visible': False},
                {'title': 'Npts:', 'name': 'Npts', 'type': 'list',
                 'values': [2 ** n for n in range(8, 16)], 'value': 512,
                 'tip': 'Number of points to display the frequency axis'},
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

                elif param.name() == 'autolimits':
                    if param.value():
                        self.settings.child('plot_settings', 'autolim_thresh').show()
                        self.settings.child('plot_settings', 'limit_min').hide()
                        self.settings.child('plot_settings', 'limit_max').hide()
                        self.settings.child('plot_settings', 'setlimits').setValue(False)

                elif param.name() == 'setlimits':
                    if param.value():
                        self.settings.child('plot_settings', 'autolim_thresh').hide()
                        self.settings.child('plot_settings', 'limit_min').show()
                        self.settings.child('plot_settings', 'limit_max').show()
                        self.settings.child('plot_settings', 'autolimits').setValue(False)

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

                elif param.name() == 'shaping_type':
                    if param.value() == 'Taylor':
                        self.settings.child('pulse_settings', 'taylor_phase').show()
                        self.settings.child('pulse_settings', 'gaussian_phase').hide()
                    elif param.value() == 'Gaussian':
                        self.settings.child('pulse_settings', 'gaussian_phase').show()
                        self.settings.child('pulse_settings', 'taylor_phase').hide()

                elif param.name() == 'npulses':
                    self.settings.child('pulse_settings', 'delay_pulses').show(param.value() > 1)


    def set_tight_layout(self, tight=True):
        self.pulse_canvas.figure.set_tight_layout(tight)
        self.trace_canvas.figure.set_tight_layout(tight)

    def show_pulse(self):
        self.update_pulse()
        self.pulse_canvas.figure.clf()
        if self.settings.child('plot_settings', 'units').value() == 'nm':
            PulsePlot(self.pulse, self.pulse_canvas.figure, xaxis='wavelength',
                      limit=self.settings.child('plot_settings', 'autolimits').value())
        else:
            PulsePlot(self.pulse, self.pulse_canvas.figure, xaxis='frequency',
                      limit=self.settings.child('plot_settings', 'autolimits').value())
        self.pulse_canvas.draw()

    def spectrum_exp(self, Npts=512, wl_lim=None):
        spectrum = normalize(lib.abs2(self.pulse.spectrum))
        wl = self.pulse.wl
        if wl_lim is not None:
            wl, spectrum = crop_vector_to_axis(wl, spectrum, wl_lim)
        wl_lin, spectrum_lin = linspace_this_vect(wl[::-1], spectrum[::-1], Npts)

        return Axis(data=wl_lin, label='Wavelength', units='m'), spectrum_lin

    def trace_exp(self, threshold=None, Npts=512, wl_lim=None):
        """ Experimental trace on linear wavelength grid of the simulated trace
        Parameters
        ----------
        threshold: (None or float)
        Npts: (int)
        wl_lim: (None or list of 2 floats)

        Returns
        -------
        meshdata: (MeshData)
        """
        md = self.pnps.trace.copy()
        md.normalize()
        md = self.get_trace_wl(md, Npts)
        md.axes[0] = md.axes[0][::-1]
        md.data = md.data[::-1, :]

        if threshold is not None:
            md.autolimit(threshold=threshold)
        elif wl_lim is not None:
            delay_c, wlc, trace_croped = crop_array_to_axis(md.axes[0], md.axes[1], md.data.T,
                                                   (np.min(md.axes[0]), np.max(md.axes[0]), wl_lim[0], wl_lim[1]))
            wl_lin, data_wl = linspace_this_image(wlc, trace_croped.T, axis=1,
                                                   Npts=Npts)
            md.data = data_wl
            md.axes[1] = wl_lin
            # md.data = trace_croped.T
            # md.axes[1] = wlc
        return md.data, Axis(data=md.axes[1], label=md.labels[1], units=md.units[1]),\
               Axis(data=md.axes[0], label=md.labels[0], units=md.units[0])

    def get_trace_wl(self, md, Npts=512):
        wl = l2w(md.axes[1] * 1e-15) * 1e-9
        wl = wl[::-1]
        md.data = md.data[:, ::-1]
        md.scale(1/wl**2)  # conversion has to be scaled by the Jacobian
        wl_lin, data_wl = linspace_this_image(wl, md.data, axis=1,
                                              Npts=Npts)

        md = MeshData(data_wl, *[md.axes[0], wl_lin], uncertainty=md.uncertainty,
                      labels=[md.labels[0], 'Wavelength'], units=[md.units[0], 'm'])
        md.normalize()
        return md

    def show_trace(self):
        self.update_pnps()
        self.trace_canvas.figure.clf()
        md = self.pnps.trace.copy()
        md.normalize()
        Npts = self.settings.child('plot_settings', 'Npts').value()
        if self.settings.child('plot_settings', 'units').value() == 'nm':
            md = self.get_trace_wl(md, Npts=Npts)

        if self.settings.child('plot_settings', 'autolimits').value():
            md.autolimit(threshold=self.settings.child('plot_settings', 'autolim_thresh').value())

        if self.settings.child('plot_settings', 'setlimits').value():
            lims = np.array([self.settings.child('plot_settings', 'limit_min').value(),
                             self.settings.child('plot_settings', 'limit_max').value()])
            if self.settings.child('plot_settings', 'units').value() == 'nm':
                lims *= 1e-9
            else:
                lims *= 1e15
            delay_c, xc, trace_croped = crop_array_to_axis(md.axes[0], md.axes[1], md.data.T,
                                                            (np.min(md.axes[0]), np.max(md.axes[0]), lims[0],
                                                             lims[1]))
            xlin, data_line = linspace_this_image(xc, trace_croped.T, axis=1, Npts=Npts)
            md.data = data_line
            md.axes[1] = xlin

        MeshDataPlot(md, self.trace_canvas.figure)
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

        if method == 'dscan':
            material = materials[self.settings.child('algo', 'material').value()]
            self.pnps = PNPS(pulse, method, process, material=material)
            parameter = linspace_step(self.settings.child('algo', 'dscan_parameter', 'min').value(),
                                      self.settings.child('algo', 'dscan_parameter', 'max').value(),
                                      self.settings.child('algo', 'dscan_parameter', 'step').value())
            parameter *= 1e-3
        elif method == 'miips':
            alpha = self.settings.child('algo', 'alpha').value()
            gamma = self.settings.child('algo', 'gamma').value()
            self.pnps = PNPS(pulse, method, process, alpha=alpha, gamma=gamma)
            parameter = linspace_step(self.settings.child('algo', 'miips_parameter', 'min').value(),
                                      self.settings.child('algo', 'miips_parameter', 'max').value(),
                                      self.settings.child('algo', 'miips_parameter', 'step').value())
        else:
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
            domega = 4 * np.log(2) / fwhm
            pulse.spectrum = gauss1D(pulse.w, x0=0., dx=domega * 1e15)  # x0=0 because the frequency axis is already

            if self.settings.child('pulse_settings', 'shaping_type').value() == 'Taylor':
                GD = self.settings.child('pulse_settings', 'taylor_phase','GD').value()
                GDD = self.settings.child('pulse_settings', 'taylor_phase', 'GDD').value()
                TOD = self.settings.child('pulse_settings', 'taylor_phase', 'TOD').value()
                phase = GD * 1e-15 * pulse.w +\
                        GDD * 1e-30 * pulse.w ** 2 / 2 +\
                        TOD * 1e-45 * pulse.w ** 3 / 6
                pulse.spectrum = pulse.spectrum * np.exp(1j * phase)
            elif self.settings.child('pulse_settings', 'shaping_type').value() == 'Gaussian':
                amp = self.settings.child('pulse_settings', 'gaussian_phase', 'gauss_amp').value()
                dtime = self.settings.child('pulse_settings', 'gaussian_phase', 'dtime').value() *1e-15
                phase = amp * gauss1D(pulse.t, 0, dtime)
                pulse.field = pulse.field * np.exp(1j * phase)

            Npulses = self.settings.child('pulse_settings', 'npulses').value()
            if Npulses > 1:
                delta_t = self.settings.child('pulse_settings', 'delay_pulses').value()
                spectrum = np.zeros_like(pulse.spectrum)
                for ind in range(Npulses):
                    spectrum += 1 / Npulses * pulse.spectrum * np.exp(1j * pulse.w * (-Npulses/2+ind) * delta_t * 1e-15)
                pulse.spectrum = spectrum

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

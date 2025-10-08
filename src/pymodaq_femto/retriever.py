import sys
import subprocess

import numpy as np
import warnings
import math, os
import inspect

from qtpy import QtWidgets, QtCore
from qtpy.QtCore import QObject, Slot, QThread, Signal, QLocale
from qtpy.QtGui import QIcon, QPixmap
from qtpy.QtGui import QTextCursor

from scipy.fftpack import next_fast_len
from scipy.interpolate import splrep, BSpline, interp1d
from collections import OrderedDict
from types import SimpleNamespace
from pathlib import Path

from pyqtgraph.dockarea import Dock
from pyqtgraph.parametertree import Parameter, ParameterTree

from pymodaq_utils import utils as utils
from pymodaq_utils.math_utils import my_moment, linspace_step
from pymodaq_utils.logger import set_logger, get_module_name, get_base_logger
from pymodaq_utils.config import Config

from pymodaq_gui.utils import DockArea
from pymodaq_gui.parameter import utils as putils, ioxml
from pymodaq_gui.h5modules.browsing import browse_data
from pymodaq_gui.h5modules.browsing import H5BrowserUtil, H5Browser
from pymodaq_gui.h5modules.saving import H5SaverLowLevel
from pymodaq_gui.utils.file_io import select_file
from pymodaq_gui.plotting.data_viewers.viewer1D import Viewer1D
from pymodaq_gui.plotting.data_viewers.viewer2D import Viewer2D
from pymodaq_gui.plotting.utils.plot_utils import RoiInfo
from pymodaq_gui.managers.action_manager import QAction
from pymodaq_gui.managers.roi_manager import LinearROI

from pymodaq_data.data import Axis, DataWithAxes, DataSource
from pymodaq_data.h5modules.data_saving import DataLoader, DataSaverLoader

from pypret import FourierTransform, Pulse, PNPS, lib, MeshData, random_gaussian
from pypret.frequencies import om2wl, wl2om, convert
from pypret.retrieval.retriever import _RETRIEVER_CLASSES


from pymodaq_femto.graphics import (
    RetrievalResultPlot,
    MplCanvas,
    NavigationToolbar,
    MeshDataPlot,
    PulsePlot,
    PulsePropagationPlot,
)
from pymodaq_femto.simulator import Simulator, methods, nlprocesses, materials, dscan_removed_message

from pymodaq_femto import _PNPS_CLASSES
import pymodaq_femto.materials

retriever_algos = list(_RETRIEVER_CLASSES.keys())

config = Config()
logger = set_logger(get_module_name(__file__))

materials_propagation = []
for item in inspect.getmembers(pymodaq_femto.materials):
    if isinstance(item[1], pymodaq_femto.materials.BaseMaterial):
        materials_propagation.append(item[1])


class DataIn(OrderedDict):
    def __init__(
        self,
        name="",
        source=DataSource.raw,
        trace_in=None,
        pulse_in=None,
        raw_spectrum=None,
        raw_trace=None,
        **kwargs
    ):
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
            raise TypeError("name for the DataIn class should be a string")
        self["name"] = name
        # if not isinstance(source, str):
        #     raise TypeError("source for the DataIn class should be a string")
        # elif not ("simulated" in source or "experimental" in source):
        #     raise ValueError('Invalid "source" for the DataIn class')
        self["source"] = source

        self["trace_in"] = trace_in
        self["pulse_in"] = pulse_in

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
    pulse.spectrum = interp1d(w - pulse.w0, spectrum, bounds_error=False, fill_value=0.0)(pulse.w)
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
        trace.labels[1] = "Frequency"
    trace.interpolate(axis2=pnps.process_w)


def substract_linear_phase(pulse):
    phase = np.unwrap(np.angle(pulse.spectrum))
    intensity = np.abs(pulse.spectrum)
    z = np.polyfit(pulse.w, phase, 1)
    pulse.spectrum *= np.exp(-1j * np.poly1d(z)(pulse.w))
    return pulse


def fit_pulse_phase(pulse, phase_blanking_threshold, order):
    phase = np.unwrap(np.angle(pulse.spectrum))
    amp = np.abs(pulse.spectrum)

    # delta = 100e-9
    x, phase = lib.mask_phase(pulse.w, amp, phase, phase_blanking_threshold)
    # fitarea = (pulse.wl > pulse.wl0-delta/2)&(pulse.wl < pulse.wl0+delta/2)
    # z = np.polyfit(pulse.w[fitarea], phase[fitarea], order)
    z = np.polyfit(x.compressed(), phase.compressed(), order)
    return z


def mask(x, y, where, **kwargs):
    y = interp1d(x[~where], y[~where], **kwargs)(x)
    return y


params_simul = Simulator.params
params_algo = utils.find_dict_in_list_from_key_val(params_simul, "name", "algo")


def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper


# method of pypret.Retriever which overwrites the _error_vector method for a modified one which uses
# an energy dependent weighting to account for unknown spectral response
def nonuniform_error_vector(self, Tmn, store=True):
    " Modified to allow for energy dependent weights"
    # rename
    rs = self._retrieval_state
    Tmn_meas = self.Tmn_meas
    # scaling factor
    w2 = self._weights * self._weights

    # mu is vector if spectral response is unknown
    mean_mu = np.sum(Tmn_meas * Tmn * w2) / np.sum(Tmn * Tmn * w2)
    mu = np.full(self.N, mean_mu)
    mask = (w2.sum(axis=0) > 0.0) & (  # weights equal to zero
            Tmn_meas.sum(axis=0) > 0.0
    )  # measurement is zero
    mu[mask] = (
            np.sum(Tmn_meas * Tmn * w2, axis=0)[mask] / np.sum(Tmn * Tmn * w2, axis=0)[mask]
    )
    # extend the edges of the response function
    idx1 = lib.find(mask, lambda x: x)
    mu[:idx1] = mu[idx1]
    idx2 = lib.find(mask, lambda x: x, n=-1)
    mu[idx2:] = mu[idx2]

    # store intermediate results in current retrieval state
    if store:
        rs.mu = mu
        rs.Tmn = Tmn
        rs.Smk = self.pnps.Smk
    return np.ravel((Tmn_meas - mu * Tmn) * self._weights)


# Optional modified retriever step calculation that keeps spectral intensity fixed
def retrieve_step_fix_spectrum(self, iteration, En):
    """ Perform a single COPRA step.

        Parameters
        ----------
        iteration : int
            The current iteration number - mainly for logging.
        En : 1d-array
            The current pulse spectrum.
        """
    # local rename
    ft = self.ft
    options = self.options
    pnps = self.pnps
    rs = self._retrieval_state
    Tmn_meas = self.Tmn_meas
    # current gradient -> last gradient
    rs.previous_max_gradient = rs.current_max_gradient
    rs.current_max_gradient = 0.0
    # switch iteration
    if rs.steps_since_improvement == 10:
        rs.mode == "local" if rs.mode == "global" else "global"
    # local iteration
    if rs.mode == "local":
        # running estimate for the trace
        Tmn = np.zeros((self.M, self.N))
        for m in np.random.permutation(np.arange(self.M)):
            p = self.parameter[m]
            Tmn[m, :] = pnps.calculate(En, p)
            Smk2 = self._project(Tmn_meas[m, :] / rs.mu, pnps.Smk)
            nablaZnm = pnps.gradient(Smk2, p)
            # calculate the step size
            Zm = lib.norm2(Smk2 - pnps.Smk)
            gradient_norm = lib.norm2(nablaZnm)
            if gradient_norm > rs.current_max_gradient:
                rs.current_max_gradient = gradient_norm
            gamma = Zm / max(rs.current_max_gradient, rs.previous_max_gradient)
            # update the spectrum
            En -= gamma * nablaZnm
            En = np.abs(self.initial_guess) * np.exp(1j * np.angle(En))
        # Tmn is only an approximation as En changed in the iteration!
        rs.approximate_error = True
        R = self._R(Tmn)  # updates rs.mu!!!
    # global iteration
    elif rs.mode == "global":
        Tmn = pnps.calculate(En, self.parameter)
        r = self._r(Tmn)
        R = self._Rr(r)  # updates rs.mu!!!
        rs.approximate_error = False
        # gradient descent w.r.t. Smk
        w2 = self._weights * self._weights
        gradrmk = (
                -4
                * ft.dt
                / (ft.dw * lib.twopi)
                * ft.backward(rs.mu * ft.forward(pnps.Smk) * (Tmn_meas - rs.mu * Tmn) * w2)
        )
        etar = options.alpha * r / lib.norm2(gradrmk)
        Smk2 = pnps.Smk - etar * gradrmk
        # gradient descent w.r.t. En
        nablaZn = pnps.gradient(Smk2, self.parameter).sum(axis=0)
        # calculate the step size
        Z = lib.norm2(Smk2 - pnps.Smk)
        etaz = options.alpha * Z / lib.norm2(nablaZn)
        # update the spectrum
        En -= etaz * nablaZn
        En = np.abs(self.initial_guess) * np.exp(1j * np.angle(En))
    return R, En


def popup_message(title, text):
    msg = QtWidgets.QMessageBox()
    msg.setWindowTitle(title)
    msg.setText(text)
    msg.setIcon(QtWidgets.QMessageBox.Warning)
    msg.exec_()


class Retriever(QObject):
    """
    Main class initializing a DAQ_Scan module with its dashboard and scanning control panel
    """

    status_signal = Signal(str)
    retriever_signal = Signal(str)
    params_in = [
        params_algo,
        {
            "title": "Data Info",
            "name": "data_in_info",
            "type": "group",
            "children": [
                {
                    "title": "Loaded file:",
                    "name": "loaded_file",
                    "type": "text",
                    "value": "",
                    "readonly": True,
                    "tip": "Loaded trace file",
                },
                {
                    "title": "Loaded node:",
                    "name": "loaded_node",
                    "type": "str",
                    "value": "",
                    "readonly": True,
                    "tip": "Loaded node within trace file",
                },
                {
                    "title": "Trace Info",
                    "name": "trace_in_info",
                    "type": "group",
                    "children": [
                        {
                            "title": "Wl0 (nm)",
                            "name": "wl0",
                            "type": "float",
                            "value": 0,
                            "readonly": True,
                            "tip": "Central spectrum wavelength in nanometers",
                        },
                        {
                            "title": "FWHM (nm)",
                            "name": "wl_fwhm",
                            "type": "float",
                            "value": 0,
                            "readonly": True,
                            "tip": "FWHM of the spectrum in nanometers",
                        },
                        {
                            "title": "Param Size",
                            "name": "trace_param_size",
                            "type": "int",
                            "value": 0,
                            "readonly": True,
                        },
                        {
                            "title": "Wavelength Size",
                            "name": "trace_wl_size",
                            "type": "int",
                            "value": 0,
                            "readonly": True,
                        },
                        {
                            "title": "Wavelength scaling",
                            "name": "wl_scaling",
                            "type": "float",
                            "value": 1,
                            "readonly": False,
                            "tip": "Scaling to go from the Trace wavelength values to wavelength in meters",
                        },
                        {
                            "title": "Parameter scaling",
                            "name": "param_scaling",
                            "type": "float",
                            "value": 1,
                            "readonly": False,
                            "tip": "Scaling to go from the trace parameter values to delay in seconds, insertion in m (dscan) "
                                   "or phase in rad (miips)",
                        },
                    ],
                },
                {
                    "title": "Spectrum Info",
                    "name": "spectrum_in_info",
                    "type": "group",
                    "children": [
                        {
                            "title": "Wl0 (nm)",
                            "name": "wl0",
                            "type": "float",
                            "value": 0,
                            "readonly": True,
                            "tip": "Central spectrum wavelength in nanometers",
                        },
                        {
                            "title": "FWHM (nm)",
                            "name": "wl_fwhm",
                            "type": "float",
                            "value": 0,
                            "readonly": True,
                            "tip": "FWHM of the spectrum in nanometers",
                        },
                        {
                            "title": "Fourier transform duration (fs)",
                            "name": "ftl",
                            "type": "float",
                            "value": 0,
                            "readonly": True,
                            "tip": "Fourier transform duration (process spectrum to get)",
                        },
                        {
                            "title": "Wavelength Size",
                            "name": "spectrum_size",
                            "type": "int",
                            "value": 0,
                            "readonly": True,
                        },
                        {
                            "title": "Wavelength scaling",
                            "name": "wl_scaling",
                            "type": "float",
                            "value": 1,
                            "readonly": False,
                            "tip": "Scaling to go from the spectrum wavelength values to wavelength in meters",
                        },
                    ],
                },
            ],
        },
        {
            "title": "Processing",
            "name": "processing",
            "type": "group",
            "children": [
                {
                    "title": "Grid settings:",
                    "name": "grid_settings",
                    "type": "group",
                    "children": [
                        {
                            "title": "lambda0 (nm):",
                            "name": "wl0",
                            "type": "float",
                            "value": 750,
                            "tip": "Central Wavelength of the Pulse spectrum and frequency grid",
                        },
                        {
                            "title": "Npoints:",
                            "name": "npoints",
                            "type": "list",
                            "limits": [2 ** n for n in range(8, 16)],
                            "value": 1024,
                            "tip": "Number of points for the temporal and Fourier Transform Grid",
                        },
                        {
                            "title": "Time resolution (fs):",
                            "name": "time_resolution",
                            "type": "float",
                            "value": 1.0,
                            "tip": "Time spacing between 2 points in the time grid",
                        },
                    ],
                },
                {
                    "title": "Trace limits:",
                    "name": "ROIselect",
                    "type": "group",
                    "visible": True,
                    "children": [
                        {
                            "title": "Crop Trace?:",
                            "name": "crop_trace",
                            "type": "bool",
                            "value": False,
                        },
                        {
                            "title": "x0:",
                            "name": "x0",
                            "type": "int",
                            "value": 0,
                            "min": 0,
                        },
                        {
                            "title": "y0:",
                            "name": "y0",
                            "type": "int",
                            "value": 0,
                            "min": 0,
                        },
                        {
                            "title": "width:",
                            "name": "width",
                            "type": "int",
                            "value": 10,
                            "min": 1,
                        },
                        {
                            "title": "height:",
                            "name": "height",
                            "type": "int",
                            "value": 10,
                            "min": 1,
                        },
                    ],
                },
                {
                    "title": "Substract trace background:",
                    "name": "linearselect",
                    "type": "group",
                    "visible": True,
                    "children": [
                        {
                            "title": "Substract?:",
                            "name": "dosubstract",
                            "type": "bool",
                            "value": False,
                        },
                        {"title": "wl0:", "name": "wl0", "type": "float", "value": 0.0},
                        {
                            "title": "wl1:",
                            "name": "wl1",
                            "type": "float",
                            "value": 10.0,
                        },
                    ],
                },
                {
                    "title": "Substract spectrum background:",
                    "name": "linearselect_spectrum",
                    "type": "group",
                    "visible": True,
                    "children": [
                        {
                            "title": "Substract?:",
                            "name": "dosubstract_spectrum",
                            "type": "bool",
                            "value": False,
                        },
                        {
                            "title": "wl0:",
                            "name": "wl0_s",
                            "type": "float",
                            "value": 0.0,
                        },
                        {
                            "title": "wl1:",
                            "name": "wl1_s",
                            "type": "float",
                            "value": 10.0,
                        },
                    ],
                },
                {
                    "title": "Process Spectrum",
                    "name": "process_spectrum",
                    "type": "action",
                    "tip": "Use ROIs to select frequency areas from the spectrum that should be removed",
                },
                {
                    "title": "Process trace",
                    "name": "process_trace",
                    "type": "action",
                    "tip": "Use one ROI to select a frequency area from the trace in order to remove the background"
                           " and use ROISelect to reduce the area around the trace",
                },
                {
                    "title": "Process Both",
                    "name": "process_both",
                    "type": "action",
                    "tip": "Process both the trace and the spectrum",
                },
            ],
        },
        {
            "title": "Retrieving",
            "name": "retrieving",
            "type": "group",
            "children": [
                {
                    "title": "Algo type:",
                    "name": "algo_type",
                    "type": "list",
                    "limits": retriever_algos,
                    "tip": "Retriever Algorithm",
                },
                {
                    "title": "Verbose Info:",
                    "name": "verbose",
                    "type": "bool",
                    "value": True,
                    "tip": "Display infos during retrieval",
                },
                {
                    "title": "Max iteration:",
                    "name": "max_iter",
                    "type": "int",
                    "value": 30,
                    "tip": "Max iteration for the algorithm",
                },
                {
                    "title": "Uniform spectral response:",
                    "name": "uniform_response",
                    "type": "bool",
                    "value": True,
                    "tip": "Assume uniform response of non-linear process. Turn off for real data.",
                },
                {
                    "title": "Keep spectral intensity fixed",
                    "name": "fix_spectrum",
                    "type": "bool",
                    "value": False,
                    "tip": "When true, only lets the phase evolve during the algorithm",
                },
                {
                    "title": "Initial guess:",
                    "name": "guess_type",
                    "type": "list",
                    "limits": ["Random gaussian", "Fundamental spectrum"],
                    "tip": "Retriever Algorithm",
                },
                {
                    "title": "Initial Pulse Guess",
                    "name": "pulse_guess",
                    "type": "group",
                    "visible": True,
                    "children": [
                        {
                            "title": "FWHM (fs):",
                            "name": "fwhm",
                            "type": "float",
                            "value": 5.0,
                            "tip": "Guess of the pulse duration (used as a starting point)",
                        },
                        {
                            "title": "Phase amp. (rad):",
                            "name": "phase_amp",
                            "type": "float",
                            "value": 0.1,
                            "tip": "Amplitude of the random phase applied to the initial guess",
                        },
                    ],
                },
                {
                    "title": "Start Retrieval",
                    "name": "start",
                    "type": "action",
                    "tip": "Start the retrieval process",
                },
                {
                    "title": "Stop Retrieval",
                    "name": "stop",
                    "type": "action",
                    "tip": "Stop the retrieval process",
                },
                {
                    "title": "Propagate result",
                    "name": "propagate",
                    "type": "action",
                    "tip": "Propagate the retrieved pulse",
                },
            ],
        },
    ]
    material_names = [material.name for material in materials_propagation]
    index = material_names.index("Air")
    material_names.insert(0, material_names.pop(index))
    air_first_material_names = material_names.copy()

    index = material_names.index("FS")
    material_names.insert(0, material_names.pop(index))

    prop_param = [
        {
            "title": "Materials",
            "name": "materials",
            "type": "group",
            "children": [
                {
                    "title": "Material 1:",
                    "name": "material1",
                    "type": "list",
                    "limits": air_first_material_names,
                    "readonly": False,
                    "tip": "First material",
                },
                {
                    "title": "Thickness (mm)",
                    "name": "thickness1",
                    "type": "float",
                    "value": 0,
                    "readonly": False,
                },
                {
                    "title": "Material 2:",
                    "name": "material2",
                    "type": "list",
                    "limits": material_names,
                    "readonly": False,
                    "tip": "Second material",
                },
                {
                    "title": "Thickness (mm)",
                    "name": "thickness2",
                    "type": "float",
                    "value": 0,
                    "readonly": False,
                },
                {
                    "title": "Plot oversampling",
                    "name": "prop_oversampling",
                    "type": "int",
                    "value": 2,
                    "readonly": False,
                },
                {
                    "title": "Resolution for FWHM calc (fs)",
                    "name": "dt_fwhm",
                    "type": "float",
                    "value": 0.5,
                    "readonly": False,
                },
                {
                    "title": "Phase masking threshold",
                    "name": "fit_threshold",
                    "type": "float",
                    "value": 0.1,
                    "readonly": False,
                },
            ],
        }
    ]
    pulse_prop = [
        {
            "title": "Pulse properties",
            "name": "pulse_prop",
            "type": "group",
            "children": [
                {
                    "title": "FWHM (fs)",
                    "name": "fwhm_meas",
                    "type": "float",
                    "value": 0.0,
                    "readonly": True,
                    "tip": "Full width at half maximum of propagated pulse",
                },
                # {'title': 'Fourier Limit (fs)', 'name': 'fwhm_ftl', 'type': 'float', 'limits': 0,
                # 'readonly': True,
                # 'tip': 'Full width at half maximum of fourier transformed pulse'},
                # {'title': 'Peak intensity compared to FTL (%)', 'name': 'ratio_main_pulse', 'type': 'float', 'limits': 0.0,
                # 'readonly': True,
                # 'tip': 'Peak intensity compared to the Fourier transform limited pulse'},
                {
                    "title": "GDD (fs^2)",
                    "name": "gdd",
                    "type": "float",
                    "value": 0.0,
                    "readonly": True,
                    "tip": "GDD",
                },
                {
                    "title": "TOD (fs3)",
                    "name": "tod",
                    "type": "float",
                    "value": 0.0,
                    "readonly": True,
                    "tip": "TOD",
                },
                {
                    "title": "FOD (fs4)",
                    "name": "fod",
                    "type": "float",
                    "value": 0.0,
                    "readonly": True,
                    "tip": "FOD",
                },
                {
                    "title": "Ratio in main pulse (%)",
                    "name": "ratio",
                    "type": "float",
                    "value": 0.0,
                    "readonly": True,
                    "tip": "Ratio",
                },
            ],
        }
    ]

    def __init__(self, dockarea=None, dashboard=None):
        """

        Parameters
        ----------
        dockarea: (dockarea) instance of the modified pyqtgraph Dockarea (see daq_utils)
        dashboard: (DashBoard) instance of the pymodaq dashboard
        """
        QLocale.setDefault(QLocale(QLocale.English, QLocale.UnitedStates))
        logger.info("Initializing Retriever Extension")
        super().__init__()

        self.h5utils = H5BrowserUtil()

        self.dockarea = dockarea
        self.dashboard = dashboard
        self.mainwindow = self.dockarea.parent()

        self.settings = Parameter.create(
            name="dataIN_settings", type="group", children=self.params_in
        )

        self.prop_settings = Parameter.create(
            name="propagation_settings", type="group", children=self.prop_param
        )
        self.pulse_settings = Parameter.create(
            name="pulse_settings", type="group", children=self.pulse_prop
        )

        self.resources_dir = os.path.abspath(os.path.dirname(__file__)) + "\\resources"

        self.setupUI()
        self.create_menu(self.mainwindow.menuBar())
        self.simulator = None
        self.data_in = None
        self.ft = None
        self.retriever = None
        self.pnps = None
        self.retriever_thread = None
        self.propagated_pulse = None
        self.result = None
        self.save_file_pathname = None
        self.state = []
        self.fake_fundamental = False

        self.settings.child("processing", "process_trace").sigActivated.connect(
            self.process_trace
        )
        self.settings.child("processing", "process_spectrum").sigActivated.connect(
            self.process_spectrum
        )
        self.settings.child("processing", "process_both").sigActivated.connect(
            self.process_both
        )
        self.settings.child("retrieving", "start").sigActivated.connect(
            self.start_retriever
        )
        self.settings.child("retrieving", "stop").sigActivated.connect(
            self.stop_retriever
        )
        self.settings.child("retrieving", "propagate").sigActivated.connect(
            self.propagate
        )

        self.settings.sigTreeStateChanged.connect(self.settings_changed)
        self.prop_settings.sigTreeStateChanged.connect(self.prop_settings_changed)

        self.viewer_trace_in.roi_select_signal.connect(self.update_ROI)
        self.viewer_trace_in.get_action('ROIselect').triggered.connect(self.show_ROI)

        self.settings.child("algo", "miips_parameter").hide()
        self.settings.child("algo", "dscan_parameter").hide()
        self.settings.child("algo", "alpha").hide()
        self.settings.child("algo", "gamma").hide()

    #################################
    # GUI Functions
    #################################
    def create_menu(self, menubar):
        """
            Create the menubar object looking like :
        """
        menubar.clear()

        # %% create Settings menu
        self.main_menu = menubar.addMenu("Main")
        self.quit_action = self.main_menu.addAction("Quit")
        self.restart_action = self.main_menu.addAction("Restart")
        self.quit_action.triggered.connect(self.quit_fun)
        self.restart_action.triggered.connect(self.restart_fun)

        self.data_in_menu = menubar.addMenu("Data In")
        self.data_in_menu.addAction(self.load_trace_in_action)
        self.data_in_menu.addAction(self.load_spectrum_in_action)
        self.data_in_menu.addSeparator()
        self.data_in_menu.addAction(self.gen_trace_in_action)
        self.data_in_menu.addAction(self.load_from_simulation_action)

        self.io_menu = menubar.addMenu("IO")
        self.io_menu.addAction(self.save_data_action)

    def setupUI(self):
        self.ui = QObject()

        #  create main docks

        self.ui.dock_settings = Dock("Settings")
        self.dockarea.addDock(self.ui.dock_settings, "top")

        self.ui.dock_data_in = Dock("Data In")
        self.dockarea.addDock(self.ui.dock_data_in, "left", self.ui.dock_settings)

        self.ui.dock_processed = Dock("Processed Data")
        self.dockarea.addDock(self.ui.dock_processed, "below", self.ui.dock_data_in)

        self.ui.dock_retriever = Dock("Retriever")
        self.dockarea.addDock(self.ui.dock_retriever, "below", self.ui.dock_processed)

        self.ui.dock_retrieved_data = Dock("Retrieved Data")
        self.dockarea.addDock(
            self.ui.dock_retrieved_data, "below", self.ui.dock_retriever
        )

        self.ui.dock_propagation = Dock("Propagation")
        self.dockarea.addDock(self.ui.dock_propagation, "below", self.ui.dock_retrieved_data)

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

        if self.dashboard is not None:
            if self.dashboard.scan_module is not None:
                self.load_last_scan_action = QAction(
                    QIcon(QPixmap(":/icons/Icon_Library/Open_2D.png")),
                    "Load last 2D scan",
                )
                self.toolbar.addAction(self.load_last_scan_action)
                self.toolbar.addSeparator()
                self.load_last_scan_action.triggered.connect(self.load_last_scan)

        self.load_trace_in_action = QAction("Load Experimental Trace")
        self.load_trace_in_action.set_icon("Open_2D")

        self.load_spectrum_in_action = QAction("Load Experimental Spectrum")
        self.load_spectrum_in_action.set_icon("Open_1D")

        self.gen_trace_in_action = QAction("Simulate Experimental Trace")
        self.gen_trace_in_action.set_icon("ini")

        self.load_from_simulation_action = QAction("Load Data from Simulation")
        self.load_from_simulation_action.set_icon("Open_sim")

        self.save_data_action = QAction("Save Data")
        self.save_data_action.set_icon("Save")

        self.save_settings_action = QAction(
            QIcon(QPixmap(os.path.join(self.resources_dir, 'save_settings.png'))),
            "Save current settings",
        )
        self.recall_settings_action = QAction(
            QIcon(QPixmap(os.path.join(self.resources_dir, 'load_settings.png'))),
            "Recall saved settings",
        )

        self.show_log_action = QAction("Show log file")
        self.show_log_action.set_icon("information2")

        self.load_trace_in_action.triggered.connect(self.load_trace_in)
        self.load_spectrum_in_action.triggered.connect(self.load_spectrum_in)
        self.gen_trace_in_action.triggered.connect(self.open_simulator)
        self.load_from_simulation_action.triggered.connect(self.load_from_simulator)
        self.save_data_action.triggered.connect(lambda: self.save_data(None))
        self.save_settings_action.triggered.connect(self.save_settings_to_file)
        self.recall_settings_action.triggered.connect(self.recall_settings_from_file)
        self.show_log_action.triggered.connect(self.show_log)

        self.toolbar.addAction(self.load_trace_in_action)
        self.toolbar.addAction(self.load_spectrum_in_action)
        self.toolbar.addSeparator()
        self.toolbar.addAction(self.gen_trace_in_action)
        self.toolbar.addAction(self.load_from_simulation_action)
        self.toolbar.addSeparator()
        self.toolbar.addAction(self.save_data_action)
        self.toolbar.addSeparator()
        self.toolbar.addAction(self.save_settings_action)
        self.toolbar.addAction(self.recall_settings_action)
        self.toolbar.addSeparator()
        self.toolbar.addAction(self.show_log_action)

        # ######################################################
        #  setup data in dock

        data_in_splitter = QtWidgets.QSplitter()
        self.viewer_trace_in = Viewer2D(QtWidgets.QWidget())
        self.viewer_trace_in.set_gradient('red', gradient="femto")  # Change default colormap

        for key in ['red', 'green', 'blue']:  # Hides all RGB controls (not needed for a trace)
            self.viewer_trace_in.get_action(key).setVisible(False)

        pos = self.viewer_trace_in.roi_manager.viewer_widget.plotItem.vb.viewRange()[0]
        self.linear_region = LinearROI(index=0, pos=pos)
        self.linear_region.setZValue(-10)
        self.linear_region.setOpacity(0.5)
        self.linear_region.setBrush([255, 0, 0, 0])
        self.viewer_trace_in.roi_manager.viewer_widget.plotItem.addItem(
            self.linear_region
        )
        self.linear_region.sigRegionChangeFinished.connect(self.update_linear)
        self.linear_region.setVisible(False)

        self.viewer_spectrum_in = Viewer1D()
        pos = self.viewer_spectrum_in.roi_manager.viewer_widget.plotItem.vb.viewRange()[
            0
        ]
        self.linear_region_spectrum = LinearROI(index=0, pos=pos)
        self.linear_region_spectrum.setZValue(-10)
        self.linear_region_spectrum.setOpacity(0.5)
        self.linear_region_spectrum.setBrush([255, 0, 0, 0])
        self.viewer_spectrum_in.roi_manager.viewer_widget.plotItem.addItem(
            self.linear_region_spectrum
        )
        self.linear_region_spectrum.sigRegionChangeFinished.connect(
            self.update_linear_spectrum
        )
        self.linear_region_spectrum.setVisible(False)

        data_in_splitter.addWidget(self.viewer_trace_in.parent)
        data_in_splitter.addWidget(self.viewer_spectrum_in.parent)
        self.ui.dock_data_in.addWidget(data_in_splitter)

        self.settings_tree.setParameters(self.settings, showTop=False)

        # #################################################
        # setup retriever dock
        retriever_widget = QtWidgets.QSplitter()
        self.viewer_live_trace = Viewer2D(QtWidgets.QWidget())
        self.viewer_live_trace.set_gradient('red', gradient="femto")
        for key in ['red', 'green', 'blue']:  # Hides all RGB controls (not needed for a trace)
            self.viewer_trace_in.get_action(key).setVisible(False)
        self.viewer_live_trace.get_action('aspect_ratio').setChecked(False)
        self.viewer_live_trace.get_action('aspect_ratio').trigger()
        self.viewer_live_trace.get_action('aspect_ratio').trigger()
        # self.viewer_trace_in.get_action('aspect_ratio').setVisible(False) #Disable aspect ratio

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

        ##################################################

        # setup propagation dock
        prop_widget = QtWidgets.QWidget()
        prop_widget.setLayout(QtWidgets.QVBoxLayout())

        param_widget = QtWidgets.QWidget()
        param_widget.setLayout(QtWidgets.QHBoxLayout())
        param_widget.setMinimumHeight(300)

        self.prop_tree = ParameterTree()
        self.pulse_tree = ParameterTree()

        param_widget.layout().addWidget(self.prop_tree, 1)
        param_widget.layout().addWidget(self.pulse_tree, 1)
        prop_widget.layout().addWidget(param_widget, 1)

        propagated_widget = QtWidgets.QWidget()
        prop_widget.layout().addWidget(propagated_widget, 5)
        self.ui.dock_propagation.addWidget(prop_widget)

        self.prop_canvas = MplCanvas(propagated_widget, width=5, height=4, dpi=100)
        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        toolbar_prop = NavigationToolbar(self.prop_canvas, propagated_widget)

        propagated_widget.setLayout(QtWidgets.QVBoxLayout())
        propagated_widget.layout().addWidget(toolbar_prop)
        propagated_widget.layout().addWidget(self.prop_canvas)
        self.prop_tree.setParameters(self.prop_settings, showTop=False)
        self.pulse_tree.setParameters(self.pulse_settings, showTop=False)

        self.ui.dock_data_in.raiseDock()

    def show_log(self):
        """Open the log file in the default text editor"""
        import webbrowser
        webbrowser.open(get_base_logger(logger).handlers[0].baseFilename)

    def settings_changed(self, param, changes):
        for param, change, data in changes:
            path = self.settings.childPath(param)
            if change == "childAdded":
                pass
            elif change == "parent":
                pass
            elif change == "value":
                if param.name() == "method":

                    # Reupdate trace params if method changes
                    if "trace_loaded" in self.state:
                        self.update_trace_info(self.data_in["raw_trace"])

                    self.settings.child("algo", "nlprocess").setLimits(
                        list(_PNPS_CLASSES[param.value()].keys())
                    )

                    if param.value() == "miips":
                        self.settings.child("algo", "alpha").show()
                        self.settings.child("algo", "gamma").show()
                    else:
                        self.settings.child("algo", "alpha").hide()
                        self.settings.child("algo", "gamma").hide()
                    if param.value() == "dscan":
                        dscan_removed_message()
                        self.settings.child("algo", "method").setValue('frog')
                    else:
                        self.settings.child("algo", "material").hide()

                elif (param.name() in putils.iter_children(self.settings.child("processing", "ROIselect"), [])
                      and "ROIselect" in param.parent().name()):  # to be sure
                    # a param named 'y0' for instance will not collide with the y0 from the ROI
                    try:
                        self.viewer_trace_in.roi_select_signal.disconnect(
                            self.update_ROI
                        )
                    except Exception as e:
                        pass
                    if self.settings.child("processing", "ROIselect", "crop_trace").value():
                        if not self.viewer_trace_in.get_action('ROIselect').isChecked():
                            self.viewer_trace_in.get_action('ROIselect').trigger()
                            QtWidgets.QApplication.processEvents()
                        self.viewer_trace_in.view.ROIselect.setPos(
                            self.settings.child(
                                "processing", "ROIselect", "x0"
                            ).value(),
                            self.settings.child(
                                "processing", "ROIselect", "y0"
                            ).value(),
                        )
                        self.viewer_trace_in.view.ROIselect.setSize(
                            [
                                self.settings.child(
                                    "processing", "ROIselect", "width"
                                ).value(),
                                self.settings.child(
                                    "processing", "ROIselect", "height"
                                ).value(),
                            ]
                        )
                        self.viewer_trace_in.roi_select_signal.connect(self.update_ROI)
                    else:
                        if self.viewer_trace_in.get_action('ROIselect').isChecked():
                            self.viewer_trace_in.get_action('ROIselect').trigger()
                elif (
                        param.name()
                        in putils.iter_children(
                    self.settings.child("processing", "linearselect"), []
                )
                        and "linearselect" in param.parent().name()
                ):  # to be sure
                    # a param named 'y0' for instance will not collide with the y0 from the ROI
                    try:
                        self.linear_region.sigRegionChangeFinished.disconnect(
                            self.update_linear
                        )
                    except Exception as e:
                        pass
                    self.linear_region.setVisible(
                        self.settings.child(
                            "processing", "linearselect", "dosubstract"
                        ).value()
                    )
                    pos_real = (
                            np.array(
                                [
                                    self.settings.child(
                                        "processing", "linearselect", "wl0"
                                    ).value(),
                                    self.settings.child(
                                        "processing", "linearselect", "wl1"
                                    ).value(),
                                ]
                            )
                            * 1e-9
                    )

                    pos_pxl, y = self.viewer_trace_in.view.unscale_axis(
                        np.array(pos_real), np.array([0, 1])
                    )
                    self.linear_region.setPos(pos_pxl)
                    self.linear_region.sigRegionChangeFinished.connect(
                        self.update_linear
                    )

                elif (
                        param.name()
                        in putils.iter_children(
                    self.settings.child("processing", "linearselect_spectrum"), []
                )
                        and "linearselect_spectrum" in param.parent().name()
                ):  # to be sure
                    # a param named 'y0' for instance will not collide with the y0 from the ROI
                    try:
                        self.linear_region_spectrum.sigRegionChangeFinished.disconnect(
                            self.update_linear_spectrum
                        )
                    except Exception as e:
                        pass
                    self.linear_region_spectrum.setVisible(
                        self.settings.child(
                            "processing",
                            "linearselect_spectrum",
                            "dosubstract_spectrum",
                        ).value()
                    )
                    pos_real = (
                            np.array(
                                [
                                    self.settings.child(
                                        "processing", "linearselect_spectrum", "wl0_s"
                                    ).value(),
                                    self.settings.child(
                                        "processing", "linearselect_spectrum", "wl1_s"
                                    ).value(),
                                ]
                            )
                            * 1e-9
                    )

                    # pos_pxl, y = self.viewer_spectrum_in.unscale_axis(np.array(pos_real), np.array([0, 1]))
                    self.linear_region_spectrum.setPos(pos_real)
                    self.linear_region_spectrum.sigRegionChangeFinished.connect(
                        self.update_linear_spectrum
                    )


                # If trace scalings are changed, rescale trace axes
                elif param.name() in ["param_scaling", "wl_scaling"] and param.parent().name() == "trace_in_info":
                    if "trace_loaded" in self.state:
                        self.load_trace_in(fname=self.data_in["file_path"], node_path=self.data_in["node_path"])

                # If spectrum scalings are changed, reload spectrum axis
                elif param.name() == "wl_scaling" and param.parent().name() == "spectrum_in_info":
                    if "spectrum_loaded" in self.state:
                        self.load_spectrum_in(fname=self.data_in["spectrum_file_path"],
                                              node_path=self.data_in["spectrum_node_path"])

                elif param.name() == "guess_type":
                    if param.value() == "Fundamental spectrum":
                        self.settings.child("retrieving", "pulse_guess").hide()
                    elif param.value() == "Random gaussian":
                        self.settings.child("retrieving", "pulse_guess").show()

                elif param.name() == "algo_type":
                    if param.value() == "copra":
                        self.settings.child("retrieving", "fix_spectrum").show()
                    else:
                        self.settings.child("retrieving", "fix_spectrum").hide()

    def restart_fun(self, ask=False):
        ret = False
        mssg = QtWidgets.QMessageBox()
        if ask:
            mssg.setText(
                "You have to restart the application to take the modifications into account!"
            )
            mssg.setInformativeText("Do you want to restart?")
            mssg.setStandardButtons(mssg.Ok | mssg.Cancel)
            ret = mssg.exec()

        if ret == mssg.Ok or not ask:
            self.quit_fun()
            subprocess.call([sys.executable, __file__])

    def quit_fun(self):
        """

        """
        try:
            if hasattr(self, "mainwindow"):
                self.mainwindow.close()

        except Exception as e:
            logger.exception(str(e))

    #################################
    # Experimental data loading
    #################################
    def load_spectrum_in(self, fname=None, node_path=None):

        # When spectrum file has already been selected, we reload it
        if fname is not None and node_path is not None:
            self.h5utils.open_file(fname, 'r+')
            dataloader = DataLoader(self.h5utils)
            spectrum = dataloader.load_data(node_path, with_bkg=True)

        # Otherwise open the dialog box
        else:
            spectrum, fname, node_path = browse_data(
                ret_all=True,
                message="Select the node corresponding to the" "Fundamental Spectrum",
            )

        if fname == "" or spectrum is None:  # User pressed cancel
            return

        else:  # User selected a node
            if self.data_in is None:
                self.data_in = DataIn(source=DataSource.raw)

            # Wavelength axis
            scaling_wl = self.settings.child(
                "data_in_info", "spectrum_in_info", "wl_scaling"
            ).value()

            nindex = spectrum.get_axis_indexes()
            if len(nindex) > 1:
                warnings.warn(
                    "The axes of the chosen spectrum have more than one dimension, defaulting to the first one.")

            wl_axis = spectrum.get_axis_from_index(nindex[0])[0] * scaling_wl

            # Store DataWithAxis: spectrum
            #      rescale axis and update labels
            spectrum.axes = [wl_axis]
            spectrum.labels = ['Fundamental Spectrum']
            self.data_in.update(
                dict(
                    raw_spectrum=spectrum,
                    spectrum_file_path=fname,
                    spectrum_node_path=node_path,
                )
            )

            wl_axis_data = wl_axis.get_data().astype("double")

            self.settings.child("processing", "linearselect_spectrum", "wl0_s").setValue(
                np.min(wl_axis_data * 1e9)
            )
            self.settings.child("processing", "linearselect_spectrum", "wl1_s").setValue(
                np.max(wl_axis_data * 1e9)
            )
            self.update_spectrum_info(self.data_in["raw_spectrum"])
            self.display_spectrum_in()
            self.fake_fundamental = False

    def update_spectrum_info(self, raw_spectrum):
        axis_index = raw_spectrum.get_axis_indexes()[0]
        wl_axis = raw_spectrum.get_axis_from_index(axis_index)[0].get_data()
        data = raw_spectrum.get_data_index()

        wl0, fwhm = my_moment(wl_axis, data)

        self.settings.child("data_in_info", "spectrum_in_info", "wl0").setValue(
            wl0 * 1e9
        )
        self.settings.child("data_in_info", "spectrum_in_info", "wl_fwhm").setValue(
            fwhm * 1e9
        )
        self.settings.child(
            "data_in_info", "spectrum_in_info", "spectrum_size"
        ).setValue(len(data))

        self.settings.child("processing", "grid_settings", "wl0").setValue(wl0 * 1e9)
        self.state.append("spectrum_loaded")

    def display_spectrum_in(self):
        self.viewer_spectrum_in.show_data(self.data_in["raw_spectrum"])

    def create_fake_fundamental(self):
        self.fake_fundamental = True
        if self.data_in is None:
            self.data_in = DataIn(source=DataSource.raw)

        # We use central wavelength of the trace
        trace_wl = self.get_trace_in().axes[1]

        if self.settings.child("processing", "ROIselect", "crop_trace").value():
            x0 = self.settings.child("processing", "ROIselect", "x0").value()
            w = self.settings.child("processing", "ROIselect", "width").value()
            wl0 = trace_wl[int(x0+w/2)]
        else:
            wl0 = self.settings["data_in_info", "trace_in_info", "wl0"] * 1e-9
        wl_fwhm = self.settings["data_in_info", "trace_in_info", "wl_fwhm"] * 1e-9
        nlprocess = self.settings.child("algo", "nlprocess").value()

        if "shg" in nlprocess:
            wlreal = 2 * trace_wl
            wl0 *= 2
        elif "thg" in nlprocess:
            wlreal = 3 * trace_wl
            wl0 *= 3
        else:
            wlreal = trace_wl

        spectrum = DataWithAxes("Fundamental Spectrum", source=DataSource.calculated,
                                data=[np.exp(-(wlreal - wl0) ** 2 / (2 * wl_fwhm ** 2))],
                                axes=[Axis(data=wlreal, label="Wavelength", units="m")])

        self.data_in.update(
            dict(
                raw_spectrum=spectrum,
                spectrum_file_path='',
                spectrum_node_path='',
            )
        )
        self.update_spectrum_info(self.data_in["raw_spectrum"])
        self.display_spectrum_in()

    def load_trace_in(self, fname=None, node_path=None):
        try:
            # When spectrum file has already been selected, we reload it
            if fname is not None and node_path is not None:
                self.h5utils.open_file(fname, 'r+')
                dataloader = DataLoader(self.h5utils)
                trace = dataloader.load_data(node_path, with_bkg=True)

            else:
                trace, fname, node_path = browse_data(
                    ret_all=True,
                    message="Select the node corresponding to the"
                            "Characterization Trace",
                )

            if fname != "":
                self.save_file_pathname = fname
                self.settings.child("data_in_info", "loaded_file").setValue(fname)
                self.settings.child("data_in_info", "loaded_node").setValue(node_path)
                self.set_data_in_exp(trace, fname, node_path)

        except Exception as e:
            logger.exception(str(e))

    def set_data_in_exp(self, trace, fname="", node_path=""):
        if self.data_in is None:
            self.data_in = DataIn(source="experimental")

        wl = trace.get_axis_from_index(trace.sig_indexes[0])[0]
        parameter_axis = trace.get_nav_axes_with_data()[0]

        scaling_parameter = self.settings.child(
            "data_in_info", "trace_in_info", "param_scaling"
        ).value()
        scaling_wl = self.settings.child(
            "data_in_info", "trace_in_info", "wl_scaling"
        ).value()

        # wl.units = "m"
        wl.data = wl.get_data() * scaling_wl

        parameter_axis.data = parameter_axis.get_data() * scaling_parameter
        # parameter_axis.units = "p.u."

        if trace.sig_indexes[0] == 0:
            trace.axes = [wl, parameter_axis]
        else:
            trace.axes = [parameter_axis, wl]

        trace.data[0] = trace.data[0].astype(float)
        self.data_in.update(
            dict(
                raw_trace=trace,
                file_path=fname,
                node_path=node_path,
                parameter_units = parameter_axis.units
            )
        )

        self.update_trace_info(self.data_in["raw_trace"])
        self.display_trace_in()

        if "trace_loaded" not in self.state:  # We don't clear the ROIs if this is not the first loaded trace
            self.viewer_trace_in.get_action('ROIselect').trigger()

    def update_trace_info(self, raw_trace):
        wl_axis = raw_trace.get_axis_from_index(raw_trace.sig_indexes[0])[0].get_data()
        param_axis = raw_trace.get_nav_axes_with_data()[0].get_data()
        data = raw_trace.get_data_index()

        wl0, fwhm = my_moment(wl_axis, np.sum(data, 0))

        self.settings.child("data_in_info", "trace_in_info", "wl0").setValue(wl0 * 1e9)
        self.settings.child("data_in_info", "trace_in_info", "wl_fwhm").setValue(
            fwhm * 1e9
        )

        self.settings.child("data_in_info", "trace_in_info", "trace_param_size").setValue(len(param_axis))
        self.settings.child("data_in_info", "trace_in_info", "trace_wl_size").setValue(len(wl_axis))

        self.settings.child("processing", "grid_settings", "npoints").setValue(
            next_fast_len(len(wl_axis))
        )

        method = self.settings.child("algo", "method").value()
        if method in ["dscan", "miips"]:
            tres = 1  # 1 fs by default
        else:
            tres = np.mean(np.diff(param_axis)) * 1e15
        self.settings.child("processing", "grid_settings", "time_resolution").setValue(tres)
        self.state.append("trace_loaded")

    def display_trace_in(self):
        self.ui.dock_data_in.raiseDock()
        self.viewer_trace_in.show_data(self.data_in['raw_trace'])
        self.viewer_trace_in.get_action('autolevels').trigger()  # Auto scale colormap
        self.viewer_trace_in.get_action('aspect_ratio').trigger()
        self.viewer_trace_in.get_action('aspect_ratio').setVisible(False) #Disable aspect ratio

        for key in ['red', 'green', 'blue']:  # Hides all RGB controls (not needed for a trace)
            # if not key == 'red': self.viewer_trace_in.get_action(key).trigger()
            self.viewer_trace_in.get_action(key).setVisible(False)
        if not self.viewer_trace_in.is_action_checked('histo'):
            self.viewer_trace_in.get_action('histo').trigger()

    def load_last_scan(self):
        try:
            viewer = self.dashboard.scan_module.ui.scan2D_graph

            data = self.dashboard.scan_module.scan_data_2D[0].T.copy()

            parameter_axis = Axis(
                data=viewer.x_axis.axis_data(data.shape[0]),
                label=viewer.x_axis.axis_label,
                units=viewer.x_axis.axis_units
            )
            wl = Axis(
                data=viewer.y_axis.axis_data(data.shape[1]),
                label=viewer.y_axis.axis_label,
                units=viewer.y_axis.axis_units
            )

            self.set_data_in_exp(data, wl, parameter_axis)

        except Exception as e:
            logger.exception(str(e))

    #################################
    # Loading from simulator
    #################################

    def open_simulator(self):
        simulator_widget = QtWidgets.QWidget()
        self.simulator = Simulator(simulator_widget)
        simulator_widget.setWindowTitle("PyMoDAQ Femto Simulator")
        simulator_widget.show()

    def load_from_simulator(self):
        if self.simulator is not None:
            data, axis, parameter_axis = self.simulator.trace_exp(Npts=512)
            spectrum_axis, spectrum_data = self.simulator.spectrum_exp(Npts=512)
            if self.data_in is None:
                self.data_in = DataIn(source="simulated")

            axis.index = 1
            trace = DataWithAxes('simulator_trace', source=DataSource.calculated, data=[data], axes=[axis, parameter_axis],
                                 nav_indexes=(0,))
            spectrum = DataWithAxes('simulator_spectrum', source=DataSource.calculated, data=[spectrum_data],
                                    axes=[spectrum_axis])
            self.data_in.update(
                dict(
                    source="simulated",
                    raw_trace=trace,
                    raw_spectrum=spectrum,
                    parameter_units=parameter_axis.units
                ))

            self.display_trace_in()
            self.display_spectrum_in()
            self.update_spectrum_info(self.data_in["raw_spectrum"])
            self.update_trace_info(self.data_in["raw_trace"])

            for child in putils.iter_children_params(
                    self.settings.child("algo"), childlist=[]
            ):
                path = ["algo"]
                path.extend(self.settings.child("algo").childPath(child))
                self.settings.child(*path).setValue(
                    self.simulator.settings.child(*path).value()
                )
            self.settings.child("data_in_info", "loaded_file").setValue("Simulation")

    #################################
    # ROI on DataIn
    #################################
    def show_ROI(self):
        # self.settings.child('processing', 'ROIselect').setOpts(
        #     visible=self.viewer_trace_in.view.ROIselect_action.isChecked())
        if "trace_loaded" in self.state:
            data = self.data_in["raw_trace"].get_data_index()
            axes = [np.arange(0, data.shape[0]), np.arange(0, data.shape[1])]
            axes_index = list(range(data.ndim))
            marginals = lib.marginals(data)
            limits = []
            for index in axes_index:
                limit = lib.limit(
                    axes[index], marginals[index].astype(float), threshold=1e-2, padding=0.25
                )
                limits.append(limit)

            self.viewer_trace_in.view.ROIselect.setPos((limits[0][0], limits[1][0]))
            self.viewer_trace_in.view.ROIselect.setSize(
                (limits[0][1] - limits[0][0], limits[1][1] - limits[1][0])
            )

            self.linear_region.setPos(limits[0])

            pos = self.viewer_trace_in.view.ROIselect.pos()
            size = self.viewer_trace_in.view.ROIselect.size()

        else:
            pos = (0, 0)
            size = (1, 1)

        self.update_ROI(RoiInfo(origin=pos, size=size))

    @Slot(RoiInfo)
    def update_ROI(self, roi_info):
        self.settings.child("processing", "ROIselect", "x0").setValue(int(roi_info.origin[1]))
        self.settings.child("processing", "ROIselect", "y0").setValue(int(roi_info.origin[0]))
        self.settings.child("processing", "ROIselect", "width").setValue(
            max([1, int(roi_info.size[1])])
        )
        self.settings.child("processing", "ROIselect", "height").setValue(
            max([1, int(roi_info.size[0])])
        )

    def update_linear(self, linear_roi):
        pos = linear_roi.pos()
        pos_real, y = self.viewer_trace_in.view.scale_axis(np.array(pos), np.array([0, 1]))
        pos_real *= 1e9
        self.settings.child("processing", "linearselect", "wl0").setValue(pos_real[0])
        self.settings.child("processing", "linearselect", "wl1").setValue(pos_real[1])

    def update_linear_spectrum(self, linear_roi):
        pos = linear_roi.pos()
        self.settings.child("processing", "linearselect_spectrum", "wl0_s").setValue(
            pos[0] * 1e9
        )
        self.settings.child("processing", "linearselect_spectrum", "wl1_s").setValue(
            pos[1] * 1e9
        )

    #################################
    # Processing
    #################################

    def process_both(self):
        ret = self.process_spectrum()
        if ret == 'Success': self.process_trace()

    def process_spectrum(self):
        if "trace_loaded" not in self.state:
            popup_message("Error", "Please load a trace first!")
            return

        if "spectrum_loaded" not in self.state:
            msg = QtWidgets.QMessageBox()
            msg.setWindowTitle('No fundamental spectrum is loaded')
            msg.setText('Proceed without a fundamental spectrum?')
            msg.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.Cancel)
            msg.setDefaultButton(QtWidgets.QMessageBox.Yes)
            msg.setIcon(QtWidgets.QMessageBox.Question)
            answer = msg.exec_()

            if answer == QtWidgets.QMessageBox.Yes:
                self.fake_fundamental = True
                self.create_fake_fundamental()
            else:
                self.fake_fundamental = False
                return 'Aborted'

        self.ui.dock_processed.raiseDock()

        self.generate_ft_grid()
        if len(np.unique(self.ft.w)) == 1:
            popup_message(
                "Error",
                "Frequency axis only has one point. Please check that i) the correct method and NL process are selected, ii) the grid settings in 'Processing' are correct. In particular, check that 'Time resolution (fs)' makes sense - typically it should be on the order of 1 fs for a standard femtosecond laser pulse.",
            )
            return 'Frequency Error'

        method = self.settings.child("algo", "method").value()
        nlprocess = self.settings.child("algo", "nlprocess").value()
        wl0 = self.settings.child("data_in_info", "trace_in_info", "wl0").value() * 1e-9
        spectrum = self.data_in["raw_spectrum"].data
        wavelength = self.data_in["raw_spectrum"].get_axis_from_index(0)[0].get_data()

        if "shg" in nlprocess:
            wl0real = 2 * wl0
        elif "thg" in nlprocess:
            wl0real = 3 * wl0
        else:
            wl0real = wl0

        self.data_in["pulse_in"] = Pulse(self.ft, wl0real)

        for roi in self.viewer_spectrum_in.roi_manager.ROIs:
            range = self.viewer_spectrum_in.roi_manager.ROIs[roi].pos()
            spectrum = mask(
                wavelength,
                spectrum,
                (range[0] <= wavelength) & (wavelength <= range[1]),
            )

        if self.settings.child(
                "processing", "linearselect_spectrum", "dosubstract_spectrum"
        ).value():
            x1 = (
                    self.settings.child(
                        "processing", "linearselect_spectrum", "wl0_s"
                    ).value()
                    * 1e-9
            )
            x2 = (
                    self.settings.child(
                        "processing", "linearselect_spectrum", "wl1_s"
                    ).value()
                    * 1e-9
            )

            idx1 = np.argmin(np.abs(wavelength - x1))
            idx2 = np.argmin(np.abs(wavelength - x2))

            if idx1 > idx2:
                idx1, idx2 = idx2, idx1

            spectrum -= np.mean(spectrum[idx1:idx2])

        self.data_in["pulse_in"] = pulse_from_spectrum(
            wavelength, spectrum, pulse=self.data_in["pulse_in"]
        )

        self.settings.child("data_in_info", "spectrum_in_info", "ftl").setValue(
            self.data_in["pulse_in"].fwhm(dt=0.1)
        )
        # self.pnps = PNPS(self.data_in['pulse_in'], method, nlprocess)

        if method == "dscan":
            material = materials[self.settings.child("algo", "material").value()]
            self.pnps = PNPS(
                self.data_in["pulse_in"], method, nlprocess, material=material
            )
            parameter = linspace_step(
                self.settings.child("algo", "dscan_parameter", "min").value(),
                self.settings.child("algo", "dscan_parameter", "max").value(),
                self.settings.child("algo", "dscan_parameter", "step").value(),
            )
            parameter *= 1e-3
        elif method == "miips":
            alpha = self.settings.child("algo", "alpha").value()
            gamma = self.settings.child("algo", "gamma").value()
            self.pnps = PNPS(
                self.data_in["pulse_in"], method, nlprocess, alpha=alpha, gamma=gamma
            )
            parameter = linspace_step(
                self.settings.child("algo", "miips_parameter", "min").value(),
                self.settings.child("algo", "miips_parameter", "max").value(),
                self.settings.child("algo", "miips_parameter", "step").value(),
            )
        else:
            self.pnps = PNPS(self.data_in["pulse_in"], method, nlprocess)

        self.state.append("spectrum_processed")
        self.pulse_canvas.figure.clf()

        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', UserWarning)
                PulsePlot(self.data_in["pulse_in"], self.pulse_canvas.figure)
        except ValueError:
            popup_message("Error",
                          "The wavelength axis of the processed spectrum seems to be wrong. Please check that i) the correct method and NL methods are selected, ii) the grid settings in 'Processing' are correct. In particular, check that 'Time resolution (fs)' makes sense - typically it is on the order of 1 fs for a standard femtosecond laser pulse.")
        self.pulse_canvas.draw()

        return 'Success'

    def process_trace(self):
        if "trace_loaded" not in self.state:
            popup_message("Error", "Please load a trace first!")
            return
        if "spectrum_processed" not in self.state:
            popup_message("Error", "Please process the spectrum first")
            return
        self.ui.dock_processed.raiseDock()

        if self.pnps is None:
            logger.info("PNPS is not yet defined, process the spectrum first")
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

        if self.settings.child("processing", "linearselect", "dosubstract").value():
            xlim = (
                    np.array(
                        (
                            self.settings.child(
                                "processing", "linearselect", "wl0"
                            ).value(),
                            self.settings.child(
                                "processing", "linearselect", "wl1"
                            ).value(),
                        )
                    )
                    * 1e-9
            )
            trace_in = preprocess(
                trace_in, signal_range=None, dark_signal_range=tuple(xlim)
            )

        if self.settings.child("processing", "ROIselect", "crop_trace").value():
            x0 = self.settings.child("processing", "ROIselect", "x0").value()
            y0 = self.settings.child("processing", "ROIselect", "y0").value()
            width = self.settings.child("processing", "ROIselect", "width").value()
            height = self.settings.child("processing", "ROIselect", "height").value()
            xlim_pxls = np.array([x0, x0 + width])
            ylim_pxls = np.array([y0, y0 + height])
            xlim, ylim = self.viewer_trace_in.view.scale_axis(xlim_pxls, ylim_pxls)
            trace_in = preprocess(trace_in, signal_range=(tuple(ylim), tuple(xlim)))

        self.data_in["trace_in"] = trace_in
        preprocess2(self.data_in["trace_in"], self.pnps)
        self.state.append("trace_processed")

        self.trace_canvas.figure.clf()
        MeshDataPlot(trace_in, self.trace_canvas.figure, limit=True)
        self.trace_canvas.draw()

    def generate_ft_grid(self):
        wl_center = self.settings.child("processing", "grid_settings", "wl0").value() * 1e-9
        Npts = self.settings.child("processing", "grid_settings", "npoints").value()
        dt = (
                self.settings.child(
                    "processing", "grid_settings", "time_resolution"
                ).value()
                * 1e-15
        )
        dw = np.pi / (0.5 * Npts * dt)

        self.ft = FourierTransform(Npts, dt, w0=wl2om(wl_center) - np.floor(Npts / 2) * dw)

    def get_trace_in(self):
        method = self.settings.child("algo", "method").value()
        if method == "dscan":
            label = "Insertion"
            unit = "m"
        elif method == "miips":
            label = "Phase"
            unit = "rad"
        else:
            label = "Delay"
            unit = "s"

        self.data_in["trace_in"] = MeshData(
            self.data_in["raw_trace"].get_data_index(),
            self.data_in["raw_trace"].get_nav_axes_with_data()[0].get_data(),
            self.data_in["raw_trace"].get_axis_from_index(self.data_in["raw_trace"].sig_indexes[0])[0].get_data(),
            labels=[label, "wavelength"],
            units=[unit, "m"],
        )

        return self.data_in["trace_in"]

    def get_pulse_in(self):
        self.data_in["pulse_in"] = pulse_from_spectrum(
            self.data_in["raw_spectrum"].get_axis_from_index(0)[0].get_data(),
            self.data_in["raw_spectrum"].get_data_index(),
        )

    def update_fwhm(self):
        precision = 1e-15 * self.prop_settings.child("materials", "dt_fwhm").value()
        try:
            fwhm = 1e15 * self.propagated_pulse.fwhm(precision)
            self.pulse_settings.child("pulse_prop", "fwhm_meas").setValue(
                truncate(fwhm, 4)
            )

        except ValueError:
            warnings.warn("FWHM is undefined.")
            self.pulse_settings.child("pulse_prop", "fwhm_meas").setValue(0)

    #################################
    # Retriever
    #################################

    def start_retriever(self):
        if "trace_processed" not in self.state:
            popup_message("Error", "Please process the trace first!")
            return
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

        self.retriever_signal.emit("start")
        self.state.append("retrieving")

    def stop_retriever(self):
        if "retrieving" not in self.state:
            popup_message("Error", "No retrieval running")
            return
        self.retriever_signal.emit("stop")

    def update_retriever_info(self, info):
        self.info_widget.moveCursor(QTextCursor.End)
        self.info_widget.insertPlainText(info + "\n")
        self.info_widget.moveCursor(QTextCursor.End)

    @Slot(list)
    def update_retriever(self, args):
        # max = 0.8 * np.max([np.abs(np.max(args[0])), np.abs(np.min(args[0]))])
        self.viewer_live_trace.set_gradient('red', gradient="femto")
        self.viewer_live_trace.show_data(DataWithAxes('retriever',
                                                      source=DataSource.calculated,
                                                      data=[args[0]],
                                                      axes=[
                                                          Axis(data=args[1], label="Parameter", units=self.data_in["parameter_units"], index=0),
                                                          Axis(data=args[2], label="Frequency", units="Hz", index=1)],
                                                      nav_indexes=(0,)))
        self.viewer_live_trace.get_action('autolevels').trigger()
        # self.viewer_live_trace.get_action('aspect_ratio').trigger()

        self.data_in["pulse_in"].spectrum = args[3]
        # self.data_in['pulse_in'] = substract_linear_phase(self.data_in['pulse_in'])
        self.viewer_live_time.show_data(DataWithAxes('retriever',
                                                     source=DataSource.calculated,
                                                     data=[np.abs(self.data_in["pulse_in"].field) ** 2],
                                                     axes=[Axis(data=self.data_in["pulse_in"].t, label="Time",
                                                                units="s")], ))
        self.viewer_live_lambda.show_data(DataWithAxes('retriever',
                                                       source=DataSource.calculated,
                                                       data=[self.data_in["pulse_in"].spectral_intensity],
                                                       axes=[Axis(data=convert(self.data_in["pulse_in"].w + self.data_in["pulse_in"].w0, "om", "wl"), label="Wavelength",
                                                                  units="m")], ))

    @Slot(SimpleNamespace)
    def display_results(self, result):
        self.result = result
        self.state.append("result_ok")
        self.ui.dock_retrieved_data.raiseDock()

        retrieved_pulse = self.data_in["pulse_in"].copy()
        retrieved_pulse.spectrum = result.pulse_retrieved
        fundamental = self.data_in["raw_spectrum"].get_data_index()
        wavelength = self.data_in["raw_spectrum"].get_axis_from_index(0)[0].get_data()
        # fundamental /= wavelength * wavelength
        spec = retrieved_pulse.spectral_intensity
        spec = interp1d(retrieved_pulse.wl, spec, bounds_error=False, fill_value=0.0)(wavelength)

        if not self.fake_fundamental:
            fundamental *= lib.best_scale(fundamental, spec)
            print("spectrum error", "%e" % lib.nrms(fundamental, spec))

        # do the retrieval plot

        self.data_canvas.figure.clf()
        RetrievalResultPlot(
            result,
            fig=self.data_canvas.figure,
            fundamental=fundamental,
            fundamental_wavelength=wavelength,
            oversampling=1,
            phase_blanking=True,
            phase_blanking_threshold=0.01,
            limit=True,
            compare_fundamental=not self.fake_fundamental
        )
        self.data_canvas.draw()

    #################################
    # Propagation
    #################################
    def prop_settings_changed(self, param, changes):
        for param, change, data in changes:
            path = self.settings.childPath(param)
            if change == "childAdded":
                pass
            elif change == "parent":
                pass
            elif change == "value":
                if param.name() in [
                    "material1",
                    "material2",
                    "thickness1",
                    "thickness2",
                    "prop_oversampling",
                    "fit_threshold",
                ]:
                    self.propagate()
                elif param.name() == "dt_fwhm":
                    self.update_fwhm()

    # def refine(self):
    #     if "result_ok" not in self.state:
    #         popup_message("Error", "Complete the retrieval first")
    #         return
    #
    #     self.propagated_pulse = Pulse(
    #         self.result.pnps.ft, self.result.pnps.w0, unit="om"
    #     )
    #     self.propagated_pulse.spectrum = self.result.pulse_retrieved
    #     tck = splrep(self.propagated_pulse.w, self.propagated_pulse.spectrum)
    #
    #     def fun(tck, pulse):
    #         mod_int = BSpline(*tck)(pulse.w)
    #         new_pulse = pulse.copy()
    #         new_pulse.spectrum = mod_int*np.abs(new_pulse.spectrum) * np.exp(1j* np.angle(new_pulse.spectrum))
    #
    #         return self.retriever.trace_error(new_pulse)
    #
    #     res = scipy.minimize(fun, tck, method='Nelder-Mead', tol=1e-6)
    #     print(res.x)
    #     mod_int = BSpline(*(res.x))(self.propagated_pulse.w)
    #     new_pulse = self.propagated_pulse.copy()
    #     new_pulse.spectrum = mod_int * np.abs(new_pulse.spectrum) * np.exp(1j * np.angle(new_pulse.spectrum))
    #
    #     self.data_canvas.figure.clf()
    #     RetrievalResultPlot(
    #         self.result,
    #         fig=self.data_canvas.figure,
    #         fundamental=fundamental,
    #         fundamental_wavelength=wavelength,
    #         oversampling=1,
    #         phase_blanking=True,
    #         phase_blanking_threshold=0.01,
    #         limit=True,
    #     )
    #     self.data_canvas.draw()

    def propagate(self):
        if "result_ok" not in self.state:
            popup_message("Error", "Complete the retrieval first")
            return
        self.ui.dock_propagation.raiseDock()

        self.propagated_pulse = Pulse(
            self.result.pnps.ft, self.result.pnps.w0, unit="om"
        )
        self.propagated_pulse.spectrum = self.result.pulse_retrieved

        material_list = []
        material_list.append(self.prop_settings.child("materials", "material1").value())
        material_list.append(self.prop_settings.child("materials", "material2").value())

        thickness_list = []
        thickness_list.append(
            self.prop_settings.child("materials", "thickness1").value()
        )
        thickness_list.append(
            self.prop_settings.child("materials", "thickness2").value()
        )

        for material, length in zip(material_list, thickness_list):
            item = getattr(pymodaq_femto.materials, material)
            w1, w2 = sorted(wl2om(np.array(item._range)))
            w = self.propagated_pulse.w + self.propagated_pulse.w0
            valid = (w >= w1) & (w <= w2)

            w = w[valid]
            k = item.k(w, unit="om")
            k0 = item.k(self.propagated_pulse.w0, unit="om")
            k1 = item.k(
                self.propagated_pulse.w0 + self.propagated_pulse.ft.dw, unit="om"
            )
            dk = (k1 - k0) / self.propagated_pulse.ft.dw

            # Add material dispersion without 0th and 1st Taylor orders (they don't change the pulse)
            kfull = np.zeros_like(self.propagated_pulse.w)
            kfull[valid] = k - k0 - dk * self.propagated_pulse.w[valid]
            self.propagated_pulse.spectrum *= np.exp(1j * kfull * 1e-3 * length)

        phasepoly = fit_pulse_phase(
            self.propagated_pulse,
            self.prop_settings.child("materials", "fit_threshold").value(),
            4,
        )
        self.propagated_pulse.spectrum *= np.exp(
            -1j * np.poly1d(phasepoly[-1])(self.propagated_pulse.w)
        )
        self.propagated_pulse.spectrum *= np.exp(
            -1j * np.poly1d(phasepoly[-2])(self.propagated_pulse.w)
        )
        self.pulse_settings.child("pulse_prop", "gdd").setValue(
            truncate(phasepoly[-3] * 1e30 * 2, 4)
        )
        self.pulse_settings.child("pulse_prop", "tod").setValue(
            truncate(phasepoly[-4] * 1e45 * 6, 4)
        )
        self.pulse_settings.child("pulse_prop", "fod").setValue(
            truncate(phasepoly[-5] * 1e60 * 24, 4)
        )

        self.update_fwhm()
        plot_oversampling = self.prop_settings.child(
            "materials", "prop_oversampling"
        ).value()

        self.prop_canvas.figure.clf()
        prop_plot = PulsePropagationPlot(
            self.propagated_pulse,
            phasepoly,
            fwhm=self.pulse_settings.child("pulse_prop", "fwhm_meas").value(),
            fig=self.prop_canvas.figure,
            oversampling=plot_oversampling,
            phase_blanking=True,
            phase_blanking_threshold=self.prop_settings.child(
                "materials", "fit_threshold"
            ).value(),
        )

        self.pulse_settings.child("pulse_prop", "ratio").setValue(
            prop_plot.intensity_fwhm.sum() / prop_plot.tamp.sum() * 100
        )
        self.prop_canvas.draw()

    #################################
    # Data saving
    #################################
    def save_data(self, save_file_pathname=None):
        try:
            if save_file_pathname is None:
                save_file_pathname = select_file(save=True, ext="h5", filter='h5 files (*.h5)')
            h5saver = H5SaverLowLevel(save_type="custom")
            h5saver.init_file(save_file_pathname,
                              raw_group_name='PyMoDAQFemto',
                              new_file=True
                              )
            datasaver = DataSaverLoader(h5saver)

            # DataIn
            data_in_group = h5saver.get_set_group(h5saver.raw_group, "DataIn")
            settings_str = b"<DataIn_settings>" + ioxml.parameter_to_xml_string(self.settings) + b"</DataIn_settings>"
            h5saver.set_attr(data_in_group, "settings", settings_str)

            if "trace_loaded" in self.state:
                trace_group = h5saver.get_set_group(data_in_group, "Trace")
                datasaver.add_data(trace_group, self.data_in["raw_trace"])

            if "spectrum_loaded" in self.state:
                spectrum_group = h5saver.get_set_group(data_in_group, "FunSpectrum")
                datasaver.add_data(spectrum_group, self.data_in["raw_spectrum"])

            # Retrieval result
            if self.result is not None:
                rr = self.result
                result_group = h5saver.get_set_group(h5saver.raw_group, "Retrieval Result")

                # Spectral Intensity
                spectrum_group = h5saver.get_set_group(result_group, "Spectral intensity")
                datasaver.add_data(spectrum_group, DataWithAxes("Spectral intensity", source=DataSource.calculated,
                                                                data=[lib.abs2(rr.pulse_retrieved)],
                                                                axes=[Axis(
                                                                    data=convert(rr.pnps.w + rr.pnps.w0, "om", "wl"),
                                                                    label="Wavelength",
                                                                    units="m")]))
                # Spectral phase
                phase_group = h5saver.get_set_group(result_group, "Spectral phase")
                datasaver.add_data(phase_group, DataWithAxes("Spectral intensity", source=DataSource.calculated,
                                                             data=[lib.phase(rr.pulse_retrieved)],
                                                             axes=[
                                                                 Axis(data=convert(rr.pnps.w + rr.pnps.w0, "om", "wl"),
                                                                      label="Wavelength",
                                                                      units="m")]))

                h5saver.set_attr(spectrum_group, "w0", rr.pnps.w0)
                h5saver.set_attr(spectrum_group, "Npts", rr.pnps.ft.N)
                h5saver.set_attr(phase_group, "w0", rr.pnps.w0)
                h5saver.set_attr(phase_group, "Npts", rr.pnps.ft.N)

                # Retrieved Trace
                trace_group_retrieved = h5saver.get_set_group(result_group, "Retrieved Trace")
                trace_retrieved = DataWithAxes("Retrieved Trace", source=DataSource.calculated,
                                               data=[rr.trace_retrieved],
                                               axes=[Axis(data=rr.parameter, label="Parameter", units="p.u.", index=0),
                                                     Axis(data=rr.pnps.process_w, label="Frequency", units="Hz",
                                                          index=1)],
                                               nav_indexes=(0,))
                datasaver.add_data(trace_group_retrieved, trace_retrieved)

            if self.propagated_pulse is not None:
                propag_group = h5saver.get_set_group(h5saver.raw_group, "Propagation")
                # Spectrum
                spectrum_prop_group = h5saver.get_set_group(propag_group, "Spectral intensity")
                datasaver.add_data(spectrum_prop_group, DataWithAxes("Spectral intensity", source=DataSource.calculated,
                                                                     data=[self.propagated_pulse.spectral_intensity],
                                                                     axes=[Axis(data=rr.pnps.ft.w, label="Frequency",
                                                                                units="Hz")]))

                phase_prop_group = h5saver.get_set_group(propag_group, "Spectral phase")
                datasaver.add_data(phase_prop_group, DataWithAxes("Spectral phase", source=DataSource.calculated,
                                                                  data=[self.propagated_pulse.spectral_phase],
                                                                  axes=[Axis(data=rr.pnps.ft.w, label="Frequency",
                                                                             units="Hz")]))
                # Time
                time_group = h5saver.get_set_group(propag_group, "Temporal intensity")
                t = np.linspace(self.propagated_pulse.t[0], self.propagated_pulse.t[-1], self.propagated_pulse.N *
                                self.prop_settings.child("materials", "prop_oversampling").value())
                datasaver.add_data(time_group, DataWithAxes("Temporal intensity", source=DataSource.calculated,
                                                            data=[self.propagated_pulse.intensity],
                                                            axes=[Axis(data=self.propagated_pulse.t, label="time",
                                                                       units="s")]))
                time_phase_group = h5saver.get_set_group(propag_group, "Temporal phase")
                datasaver.add_data(time_phase_group, DataWithAxes("Temporal phase", source=DataSource.calculated,
                                                                  data=[self.propagated_pulse.phase],
                                                                  axes=[Axis(data=self.propagated_pulse.t, label="time",
                                                                             units="s")]))

                settings_str = b'<prop_settings title="Prop. Settings" type="group">'
                settings_str += ioxml.parameter_to_xml_string(self.prop_settings)
                settings_str += ioxml.parameter_to_xml_string(self.pulse_settings)
                settings_str += b"</prop_settings>"

                h5saver.set_attr(propag_group, "settings", settings_str)

            settings_str = b'<All_settings title="All Settings" type="group">'
            settings_str += ioxml.parameter_to_xml_string(self.settings)
            settings_str += ioxml.parameter_to_xml_string(self.pulse_settings)
            settings_str += ioxml.parameter_to_xml_string(self.prop_settings)
            settings_str += b"</All_settings>"

            h5saver.set_attr(h5saver.raw_group, "settings", settings_str)

        except Exception as e:
            logger.exception(str(e))

        h5saver.close_file()

    def save_settings_to_file(self):
        path_to_file = Path(os.path.join(self.resources_dir, 'retriever_settings.h5'))

        msg = QtWidgets.QMessageBox()
        msg.setWindowTitle('Save settings to file')

        if os.path.exists(path_to_file):
            msg.setText('Do you want to overwrite the file with the current settings?')
        else:
            msg.setText('Should I save the current settings to file?')
        msg.setStandardButtons(QtWidgets.QMessageBox.Save | QtWidgets.QMessageBox.Cancel)
        msg.setDefaultButton(QtWidgets.QMessageBox.Save)
        msg.setIcon(QtWidgets.QMessageBox.Question)
        answer = msg.exec_()

        if answer == QtWidgets.QMessageBox.Save:
            self.save_data(path_to_file)

    def recall_settings_from_file(self):
        path_to_file = os.path.join(self.resources_dir, 'retriever_settings.h5')

        if not os.path.exists(path_to_file):
            popup_message("Error", "Did not find a file with saved settings.")
        else:
            self.h5utils.open_file(path_to_file, 'r+')
            dataloader = DataLoader(self.h5utils)
            spectrum = dataloader.load_data('/PyMoDAQFemto/DataIn/FunSpectrum/Data00')
            trace = dataloader.load_data('/PyMoDAQFemto/DataIn/Trace/Data00')
            settings = dataloader.get_node('/PyMoDAQFemto').attrs['settings']
            self.h5utils.close_file()

            saved_settings = ioxml.XML_string_to_parameter(settings)
            saved_param = Parameter.create(name="saved_settings", type="group", children=saved_settings)

            # Algo settings
            for child in putils.iter_children_params(self.settings.child("algo"), childlist=[]):
                path = ["algo"]
                path.extend(self.settings.child("algo").childPath(child))
                self.settings.child(*path).setValue(saved_param.child('dataIN_settings', *path).value())
            # Data Info settings
            for child in putils.iter_children_params(self.settings.child("data_in_info"), childlist=[]):
                path = ["data_in_info"]
                path.extend(self.settings.child("data_in_info").childPath(child))
                self.settings.child(*path).setValue(saved_param.child('dataIN_settings', *path).value())

            if self.data_in is None:
                self.data_in = DataIn(source=DataSource.raw)

            self.data_in.update(
                dict(
                    raw_spectrum=spectrum,
                    spectrum_file_path=path_to_file,
                    spectrum_node_path='/PyMoDAQFemto/DataIn/FunSpectrum/Data00',
                )
            )
            self.update_spectrum_info(self.data_in["raw_spectrum"])
            self.display_spectrum_in()

            # Load trace
            self.settings.child("data_in_info", "loaded_file").setValue(path_to_file)
            self.settings.child("data_in_info", "loaded_node").setValue('/PyMoDAQFemto/DataIn/Trace/Data00')
            self.set_data_in_exp(trace, path_to_file, '/PyMoDAQFemto/DataIn/Trace/Data00')

            # Processing settings
            for child in putils.iter_children_params(self.settings.child("processing"), childlist=[]):
                path = ["processing"]
                path.extend(self.settings.child("processing").childPath(child))
                self.settings.child(*path).setValue(saved_param.child('dataIN_settings', *path).value())
            # Retrieving settings
            for child in putils.iter_children_params(self.settings.child("retrieving"), childlist=[]):
                path = ["retrieving"]
                path.extend(self.settings.child("retrieving").childPath(child))
                self.settings.child(*path).setValue(saved_param.child('dataIN_settings', *path).value())

class RetrieverWorker(QObject):
    result_signal = Signal(SimpleNamespace)
    status_sig = Signal(str)
    callback_sig = Signal(list)

    def __init__(self, data_in, pnps, settings):
        super().__init__()
        self.settings = settings
        self.data_in = data_in
        self.pnps = pnps
        self.retriever = None

    # def send_callback(self, pnps):
    #     self.callback_sig.emit([pnps.Tmn, [pnps.parameter, pnps.process_w], pnps.pulse.field, pnps.field.t])

    @Slot(str)
    def command_retriever(self, command):
        if command == "start":
            self.start_retriever()
        elif command == "stop":
            self.stop_retriever()

    def start_retriever(self):
        retriever_cls = _RETRIEVER_CLASSES[
            self.settings.child("retrieving", "algo_type").value()
        ]
        verbose = self.settings.child("retrieving", "verbose").value()
        max_iter = self.settings.child("retrieving", "max_iter").value()
        fwhm = self.settings.child("retrieving", "pulse_guess", "fwhm").value()
        amplitude = self.settings.child(
            "retrieving", "pulse_guess", "phase_amp"
        ).value()
        fix_spectrum = self.settings.child("retrieving", "fix_spectrum").value()

        uniform_response = self.settings.child("retrieving", "uniform_response").value()
        preprocess2(self.data_in["trace_in"], self.pnps)

        self.retriever = retriever_cls(
            self.pnps,
            logging=True,
            verbose=verbose,
            maxiter=max_iter,
            status_sig=self.status_sig,
            callback=self.callback_sig.emit,
            step_command=QtWidgets.QApplication.processEvents,
        )
        if fix_spectrum and self.retriever.method == "copra":
            self.retriever._retrieve_step = retrieve_step_fix_spectrum.__get__(
                self.retriever
            )
        if not uniform_response:
            self.retriever._error_vector = nonuniform_error_vector.__get__(
                self.retriever
            )

        if (
                self.settings.child("retrieving", "guess_type").value()
                == "Fundamental spectrum"
        ):
            pulse_guess = self.data_in["pulse_in"].copy()
            pulse_guess.spectrum = (1 + 0 * 1j) * np.abs(
                self.data_in["pulse_in"].spectrum
            )
            pulse_guess.spectrum /= (
                    self.data_in["pulse_in"].wl * self.data_in["pulse_in"].wl
            )
            pulse_guess.field /= np.abs(pulse_guess.field).max()

            guess = pulse_guess.spectrum

        elif (
                self.settings.child("retrieving", "guess_type").value() == "Random gaussian"
        ):
            random_gaussian(self.data_in["pulse_in"], fwhm * 1e-15, phase_max=amplitude)
            guess = self.data_in["pulse_in"].spectrum

        self.retriever.retrieve(self.data_in["trace_in"], guess, weights=None)

        self.result_signal.emit(self.retriever.result())

    def stop_retriever(self):
        self.retriever._retrieval_state.running = False


def main():
    # from pymodaq.daq_utils.daq_utils import get_set_preset_path

    app = QtWidgets.QApplication(sys.argv)
    win = QtWidgets.QMainWindow()
    area = DockArea()
    win.setCentralWidget(area)
    win.resize(1000, 500)
    win.setWindowTitle("PyMoDAQ Retriever")

    prog = Retriever(dashboard=None, dockarea=area)
    win.show()
    # try:
    #     prog.load_trace_in(fname='C:\\Data\\2021\\20210315\\Dataset_20210315_001\\Dataset_20210315_001.h5',
    #                         node_path='/Raw_datas/Scan001/Detector000/Data1D/Ch000/Data')
    #     prog.load_spectrum_in(fname='C:\\Users\\weber\\Desktop\\pulse.h5',
    #                         node_path='/Raw_datas/Detector000/Data1D/Ch000/Data')
    #     prog.save_data('C:\\Users\\weber\\Desktop\\pulse_analysis.h5')
    # except:
    #     pass
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
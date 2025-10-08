import numpy as np

# Loading external numpy data
data = np.load('example_trace_numpy.npz')
trace = data["trace"]
delay = data["parameter"]
wavelength = data["wavelength"]

# Load also a fundamental spectrum (not necessary for retrieval but can be used to check)
fundamental = np.load('example_fundamental_numpy.npz')
spectrum = fundamental["spectrum"]
spectrum_wavelength = fundamental["wavelength"]


"""
Conversion Part
================
This will create a file "converted_trace.h5" that can be loaded in PyMoDAQ-Femto
This trace is a "PG-FROG" trace, i.e. the parameter axis is a delay in seconds,
the wavelength axis is in meters.

!! Make sure to select the proper options in the retriever GUI when using the trace! (especially, nlprocess = pg)
"""

from pymodaq_femto.converter import convert_numpy_to_pymodaq_femto

convert_numpy_to_pymodaq_femto("converted_trace.h5", trace, delay, wavelength, parameter_units="s")
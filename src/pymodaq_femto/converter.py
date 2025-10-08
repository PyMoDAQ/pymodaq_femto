
import numpy as np
from pymodaq_gui.h5modules.saving import H5SaverLowLevel
from pymodaq_data.h5modules.data_saving import DataSaverLoader
from pymodaq.utils.data import Axis, DataWithAxes
from pymodaq_data.h5modules.backends import SaveType
from pathlib import Path
from pymodaq_data.data import DataSource


def convert_numpy_to_pymodaq_femto(
    output_file: str,
    trace: np.ndarray,
    parameter_axis: np.ndarray,
    wavelength_axis: np.ndarray,
    parameter_units: str = "s"
) -> None:
    """
    Convert a measured trace (saved as numpy arrays) to a h5 file readable by PyMoDAQ-Femto.

    This function takes numpy arrays representing a measured trace, along with its parameter
    and wavelength axes, and saves the data into an HDF5 file format compatible with PyMoDAQ-Femto.

    All quantities (parameter, wavelength, etc.) must be expressed in standard SI units
    (meters, seconds, etc.), not in femtoseconds or other derived units.

    Parameters
    ----------
    trace : numpy.ndarray
        The 2D array representing the measured trace data.
    parameter_axis : numpy.ndarray
        The 1D array representing the parameter axis (e.g., insertion or delay), in SI units (seconds, meters, etc.).
    wavelength_axis : numpy.ndarray
        The 1D array representing the wavelength axis, in SI units (meters).
    output_file : str
        Path to the output HDF5 file where the converted data will be saved.
    parameter_units : str, optional
        Unit of the parameter axis (default: "s" for seconds). Must be a standard SI unit.

    Notes
    -----
    - The function initializes an HDF5 file using `H5SaverLowLevel` and organizes the data
      into groups for compatibility with PyMoDAQ-Femto.
    - The `DataSaverLoader` is used to manage the saving process.
    """

    # Create saver
    h5saver = H5SaverLowLevel(save_type=SaveType.custom)
    h5saver.init_file(file_name=Path(output_file),
                             raw_group_name='PyMoDAQFemto', new_file=True)
    datasaver = DataSaverLoader(h5saver)
    data_in_group = h5saver.get_set_group(h5saver.raw_group, "DataIn")
    trace_group = h5saver.get_set_group(data_in_group, "Trace")

    # store data into high level object
    measured_data = DataWithAxes("Measured Trace", source=DataSource.raw,
                         data=[trace],
                         axes=[Axis(data=parameter_axis, label="Parameter", units=parameter_units, index=0),
                               Axis(data=wavelength_axis, label="Wavelength", units="m",
                                    index=1)],
                         nav_indexes=(0,))

    # add it to file
    datasaver.add_data(trace_group, measured_data)
    h5saver.close_file()


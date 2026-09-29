PyMoDAQ Femto
#############

PyMoDAQ extension for femtosecond laser pulse characterization (FROG, d-scan, ...).
It provides two applications:

* the **Simulator**, to generate characterization traces from simulated pulses,
* the **Retriever**, to retrieve the pulse (spectral amplitude and phase) from an experimental or simulated trace,
  and propagate it through materials.

Documentation can be found here: https://pymodaq-femto.readthedocs.io/en/latest/index.html

Compatibility
=============

PyMoDAQ Femto 5.2.0 is compatible with **PyMoDAQ 5.2**. It only relies on three PyMoDAQ sub-packages,
pinned to the versions it has been tested with (see Dependencies below): the full ``pymodaq`` package is not needed
to use PyMoDAQ Femto on its own.

PyMoDAQ Femto 5.1.0 is compatible with PyMoDAQ 5.1 only: it does **not** work with PyMoDAQ 5.2 or later.

Versions of PyMoDAQ Femto compatible with PyMoDAQ 3 or 4 are archived in the other branches of this repository
(legacy_v3 and v4).

Installation
============

PyMoDAQ Femto requires Python 3.10 or later (tested with Python 3.14). We advise to install it in a dedicated
conda environment::

    conda create -n pymodaq_femto python=3.14
    conda activate pymodaq_femto
    pip install pymodaq_femto


For development, clone the repository and install it in editable mode::

    git clone https://github.com/PyMoDAQ/pymodaq_femto.git
    cd pymodaq_femto
    pip install -e .

Dependencies
------------

These are installed automatically by pip:

======================  =======  ======================================================================
Package                 Version  Used for
======================  =======  ======================================================================
``pymodaq_utils``       5.2.7    configuration, logging, math and unit utilities
``pymodaq_data``        5.2.10   data objects (``DataWithAxes``) and HDF5 files (through ``h5py``)
``pymodaq_gui``         5.2.9    data viewers, parameter trees and HDF5 browser
``pypret_pymodaq``      any      pulse retrieval algorithms (fork of `pypret`_)
``numpy``, ``scipy``    any      numerical computations
``matplotlib``          any      result plots
``PyQt5``               any      graphical user interface
``pymodaq`` (optional)  5.2.11   running the Retriever from the PyMoDAQ dashboard (``dashboard`` extra)
======================  =======  ======================================================================

The PyMoDAQ packages are pinned to exact versions because new PyMoDAQ releases have previously broken
PyMoDAQ Femto. They will be updated once newer versions have been tested.

.. _pypret: https://github.com/ncgeib/pypret

Launching PyMoDAQ Femto
=======================

Once installed, start either application from your activated environment::

    retriever
    simulator

or equivalently ``python -m pymodaq_femto.retriever`` and ``python -m pymodaq_femto.simulator``.

Usage with data external to PyMoDAQ
===================================

In the /utils/ folder, you will find `an example script <https://github.com/PyMoDAQ/pymodaq_femto/tree/main/src/pymodaq_femto/utils/example_conversion.py>`_
showing how any data (here saved as numpy arrays) can be easily converted into a properly formatted .h5 file,
and used in PyMoDAQ Femto.

License
=======

Published under the MIT license (see the LICENSE file).

GitHub repo: https://github.com/PyMoDAQ/pymodaq_femto

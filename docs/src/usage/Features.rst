Features
========

Overview of PyMoDAQ-Femto
+++++++++++++++++++++++++

The package implements several retrieval algorithms for ultrashort laser pulse measurement methods, such as as frequency-resolved
optical gating (FROG), dispersion scan (d-scan), and more. The application can simulate measurement traces from various pulse shapes,
and apply retrieval algorithms to them. It also works on real experimental measured traces.

PyMoDAQ-Femto is written in `Python`__ and uses Python 3.7+. It is an extension of the `PyMoDAQ`__ package, which
is a Modular Data Acquisition module that can interface any kind of experiment. As such, PyMoDAQ-Femto can natively work on data
measured using PyMoDAQ, although it can also work on any data provided that it is converted to the proper format (see :ref:`convertingdata` for conversion guidelines).

The algorithms implemented in PyMoDAQ-Femto are based on the excellent `pypret`__ package, which provides a common pulse
retrieval algorithm to several pulse measurement methods (see `[Geib2019]`__ for a full description).

    .. _overview:

.. figure:: /image/overview.png
   :alt: overview

   Overview of the two modules of PyMoDAQ-Femto

__ https://docs.python-guide.org/
__ http://pymodaq.cnrs.fr/en/latest/index.html
__ https://pypret.readthedocs.io/en/latest/
__ https://doi.org/10.1364/OPTICA.6.000495

.. _available_methods:

Available methods
+++++++++++++++++

.. list-table:: Available measurement types that can be simulated or retrieved, with their supported non-linear processes.
   :widths: 25 25 50
   :header-rows: 1

   * -  Method
     - Full Name
     - Supported non-linear processes [*]_
   * - frog
     - Frequency-resolved optical gating
     - shg, pg, tg
   * - dscan
     - Dispersion scan
     - shg, thg, sd
   * - ifrog
     - Interferometric frequency-resolved optical gating
     - shg, thg, sd
   * - miips
     - Multiphoton intrapulse interference phase scan
     - shg, thg, sd
   * - tdp
     - Time-domain ptychography
     - shg, thg, sd

.. [*] shg: Second Harmonic Generation, thg: Third Harmonic Generation, sd: Self Diffraction, pg: Polarization Gating, tg: Transient Grating
Welcome to PyMoDAQ-Femto documentation!
=======================================
PyMoDAQ-Femto is an application that allows the temporal characterization of laser pulses (typically of femtosecond or picosecond duration).
The module was initially developed for educational purposes in the framework of the international Femto-UP 2020 School. It has now turned
into a tool used in several research laboratories.

PyMoDAQ-Femto is a 2-in-1 python application dealing with femtosecond laser pulse characterization. It features:

* A user interface called **Simulator** to define short pulses and simulate the non-linear traces they yield. Many
  characterization techniques are available (see :ref:`available_methods`), and several non-linear processes are implemented for each of them.

* A user interface called **Retriever** to run various retrieval algorithms on simulated or experimental traces (acquired using PyMoDAQ or other means)

   .. _simulator_fig:

.. figure:: /image/simulator.png
   :alt: simulator

   PyMoDAQ-Femto's Simulator.



Both modules can be ran as stand-alone applications or plugged as an extension to `PyMoDAQ`__. All together it produces
a framework for complete temporal characterization of shaped ultrashort femtosecond pulses.

__ http://pymodaq.cnrs.fr


Information
***********

GitHub repo: https://github.com/PyMoDAQ/pymodaq_femto

Documentation: http://pymodaq_femto.cnrs.fr/

Based on PyMoDAQ, the `pypret`__ library and the ``pyqtgraph`` library.

PyMoDAQ-Femto is written by Sébastien Weber: sebastien.weber@cemes.fr and Romain Géneaux: romain.geneaux@cea.fr under a
MIT license.

__ https://github.com/ncgeib/pypret


Contribution
************

If you want to contribute see this page: :ref:`contributors`

..
    They use it
    ***********
    See :ref:`feedback`


Citation
********

The module is described at length in the following open-access paper:

    Romain Géneaux and Sébastien Weber, "Femtosecond Pulse Shaping and Characterization: From Simulation to Experimental Pulse Retrieval Using a Python-Based User Friendly Interface". In J. Léonard and C. Hirlimann (Eds.), *Ultrafast Laser Technologies and Applications* (pp. 111-128). EDP Sciences.
    `https://doi.org/10.1051/978-2-7598-2719-0`__

__ https://doi.org/10.1051/978-2-7598-2719-0

If you publish results obtained with the help of the PyMoDAQ-Femto interface, we would appreciate your using this reference.
In that way, you're also helping in its promotion and amelioration.



..
    Changelog
    *********

    Please see :doc:`the changelog </changelog>`.


Index
*****

.. toctree::
   :numbered:
   :maxdepth: 6
   :caption: Contents:

   usage/Features
   usage/Installation
   usage/Simulator
   usage/Retriever
   usage/Contributors


.. toctree::
   :hidden:
   :caption: Related packages

   PyMoDAQ <http://pymodaq.cnrs.fr>
   pypret <https://pypret.readthedocs.io/en/latest/>


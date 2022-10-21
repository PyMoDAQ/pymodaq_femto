Welcome to PyMoDAQ-Femto documentation!
=======================================

PyMoDAQ-Femto is a 2-in-1 python application dealing with femtosecond laser pulse characterization. It features:

* A user interface called **Simulator** to simulate non-linear traces obtained using most known characterization
  techniques (FROG, D-Scan, ... with their various non-linear flavour)
* a user interface called **Retriever** to run various retrieval algorithms on simulated or experimental traces
  (acquired for instance using PyMoDAQ or other means)

   .. _simulator_fig:

.. figure:: /image/simulator.png
   :alt: simulator

   PyMoDAQ-Femto's Simulator.



Both modules can be run as stand-alone application or plugged as an extension to `PyMoDAQ`__. All together it produces
a framework for complete temporal characterization of shaped ultrashort femtosecond pulses.

__ http://pymodaq.cnrs.fr


Information
***********

GitHub repo: https://github.com/PyMoDAQ/pymodaq_femto

Documentation: http://pymodaq_femto.cnrs.fr/

Based on PyMoDAQ and the `pypret`__ library the ``pyqtgraph`` library.


PyMoDAQ-Femto is written by Sébastien Weber: sebastien.weber@cemes.fr and Romain Géneaux: romain.geneaux@cea.fr under a
MIT license.

__ https://github.com/ncgeib/pypret


Contribution
************

If you want to contribute see this page: :ref:`contributors`


They use it
***********
See :ref:`feedback`


Citation
********

By using PyMoDAQ-Femto, you are being asked to cite the article published in Review of Scientific
Instruments `RSI 92, 045104 (2021)`__ when publishing results obtained with the help of its interface.
In that way, you're also helping in its promotion and amelioration.

Changelog
*********

Please see :doc:`the changelog </changelog>`.


Index
*****

.. toctree::
   :maxdepth: 6
   :caption: Contents:

   usage/Installation
   usage/Feedback
   usage/Contributors


.. toctree::
   :hidden:
   :caption: Related packages

   PyMoDAQ <http://pymodaq.cnrs.fr>
   pypret <https://pypret.readthedocs.io/en/latest/>


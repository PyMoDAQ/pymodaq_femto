  .. _section_installation:

Installation
============

.. contents::
   :depth: 1
   :local:
   :backlinks: none

.. highlight:: console

Requirements
------------

PyMoDAQ-Femto runs on **Windows**, **MacOS** and **Linux**, and requires Python 3.10 or later (it is tested with
Python 3.14). We advise to install it with `Miniconda`__ (a light package manager) or `Anaconda`__, in a dedicated
environment: this isolates PyMoDAQ-Femto and its dependencies from your other Python programs.

__ https://docs.conda.io/en/latest/miniconda.html
__ https://www.anaconda.com/download/

PyMoDAQ-Femto 5.2 is compatible with PyMoDAQ 5.2. Versions compatible with PyMoDAQ 3 or 4 are archived in the
legacy_v3 and v4 branches of the `GitHub repository`__.

__ https://github.com/PyMoDAQ/pymodaq_femto

Setting up a new environment
----------------------------

* Download and install Miniconda.
* Open a console (on Windows, the *Anaconda Prompt*).
* Create a new environment called *pymodaq_femto* (any name will do) with a recent Python version::

    conda create -n pymodaq_femto python=3.14

* Activate it, so that only the packages installed within this environment are *seen* by Python::

    conda activate pymodaq_femto

Installing PyMoDAQ-Femto
------------------------

In your activated environment, enter::

    pip install pymodaq_femto

This installs the latest version of PyMoDAQ-Femto and all its dependencies. For a specific version, enter
``pip install pymodaq_femto==x.y.z``.

PyMoDAQ-Femto does not need the full PyMoDAQ package: it only relies on three of its sub-packages
(``pymodaq_utils``, ``pymodaq_data`` and ``pymodaq_gui``), which are pinned to the versions it has been tested with.
The other dependencies are ``pypret_pymodaq`` (the retrieval algorithms), ``numpy``, ``scipy``, ``matplotlib`` and
``PyQt5``.

To also use the Retriever as an extension of the PyMoDAQ dashboard, install the ``dashboard`` option, which adds the
full PyMoDAQ package::

    pip install "pymodaq_femto[dashboard]"

For development, clone the repository and install it in editable mode::

    git clone https://github.com/PyMoDAQ/pymodaq_femto.git
    cd pymodaq_femto
    pip install -e .

  .. _run_module:

Launching PyMoDAQ-Femto
-----------------------

The installation creates two commands in your environment. With the environment activated, enter either:

*  ``simulator``
*  ``retriever``

Alternatively, you can use the full commands:

*  ``python -m pymodaq_femto.simulator``
*  ``python -m pymodaq_femto.retriever``

  .. _shortcut_section:

Creating shortcuts on **Windows**
---------------------------------

Windows users may prefer to start PyMoDAQ-Femto from a shortcut on the desktop
(thanks to Christophe Halgand for the procedure):

* Create a shortcut on your desktop, pointing to any file or program (see :numref:`shortcut_create`).
* Right click on it and open its properties (see :numref:`shortcut_prop`).
* In the *Start in* field, enter the path to the *condabin* folder of your Miniconda or Anaconda installation,
  for instance ``C:\Miniconda3\condabin``.
* In the *Target* field, enter ``C:\Windows\System32\cmd.exe /k conda activate pymodaq_femto & retriever``.
  The shortcut opens a console, activates your environment, then starts the Retriever.
* Repeat with ``simulator`` instead of ``retriever`` to create a shortcut for the Simulator.

   .. _shortcut_create:

.. figure:: /image/installation/shortcut_creation.png
   :alt: shortcut

   Create a shortcut on your desktop

   .. _shortcut_prop:

.. figure:: /image/installation/shortcut_prop.PNG
   :alt: shortcut properties

   Shortcut properties (here with the older *python -m* command)

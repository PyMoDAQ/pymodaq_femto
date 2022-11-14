  .. _section_installation:

Installation
============

.. contents::
   :depth: 1
   :local:
   :backlinks: none

.. highlight:: console

For PyMoDAQ-Femto to run smoothly, you need a Python distribution to be installed. Here are some advices.
On all platforms **Windows**, **MacOS** or **Linux**, `Anaconda`__ or `Miniconda`__ is the advised distribution/package
manager. Environments can be created to deal with different version of packages and isolate the code from other
programs. Anaconda comes with a full set of installed scientific python packages while *Miniconda* is a very
light package manager.

__ https://www.anaconda.com/download/
__ https://docs.conda.io/en/latest/miniconda.html

Setting up a new environment
----------------------------

* Download and install Miniconda3.
* Open a console, and cd to the location of the *condabin* folder, for instance: ``C:\Miniconda3\condabin``
* Create a new environment: ``conda create -n my_env python=3.8``, where my_env is your new environment name. This will create the environment with python version 3.8
  that is currently the recommended one.
* Activate your environment so that only packages installed within this environment will be *seen* by Python:
  ``conda activate my_env``

Installing PyMoDAQ-Femto
------------------------

Easiest part: in your newly created and activated environment enter: ``pip install pymodaq_femto``. This will install the
latest PyMoDAQ-Femto available version and all its dependencies. For a specific version
enter:  ``pip install pymodaq_femto==x.y.z``.

  .. _run_module:

Launching PyMoDAQ-Femto
---------------------------------

During its installation, two scripts have been installed within you environment directory,
this means you can start PyMoDAQ-Femto’s two main functionalities directly writing in your console either:

*  ``simulator``
*  ``retriever``

Alternatively, you can specify the full commands (The *-m* option tells python to look within its *site-packages* folder, where you've just
installed pymodaq_femto):

*  ``python -m pymodaq_femto.simulator``
*  ``python -m pymodaq_femto.retriever``

  .. _shortcut_section:

Creating shortcuts on **Windows**
---------------------------------

Python packages can easily be started from the command line (see :ref:`section_how_to_start`). However, Windows users
will probably prefer using shortcuts on the desktop. Here is how to do it (Thanks to Christophe Halgand for the procedure):

* First create a shortcut (see :numref:`shortcut_create`) on your desktop (pointing to any file or program, it doesn't matter)
* Right click on it and open its properties (see :numref:`shortcut_prop`)
* On the *Start in* field ("Démarrer dans" in french and in the figure), enter the path to the condabin folder of your miniconda or
  anaconda distribution, for instance: ``C:\Miniconda3\condabin``
* On the *Target* field, ("Cible" in french and in the figure), enter this string:
  ``C:\Windows\System32\cmd.exe /k conda activate my_env & python -m pymodaq_femto.retriever``. This means that
  your shortcut will open the windows's command line, then execute your environment activation (*conda activate my_env* bit),
  then finally execute and start **Python**, opening the correct pymodaq_femto file (here *retriever.py*,
  starting the Retriever module, *python -m pymodaq_femto.retriever* bit)
* You're done!
* Do it again for each PyMoDAQ-Femto's module you want (to get the correct python file and it's path, see :ref:`run_module`).



   .. _shortcut_create:

.. figure:: /image/installation/shortcut_creation.png
   :alt: shortcut

   Create a shortcut on your desktop

   .. _shortcut_prop:

.. figure:: /image/installation/shortcut_prop.PNG
   :alt: shortcut properties

   Shortcut properties

  .. _retriever:

Retriever module
================
.. |cemes| image:: /image/logos/logo_cemes.png
   :width: 100
   :alt: cemes

.. |attolab| image:: /image/logos/attolab_logo_carre.jpg
   :width: 100
   :alt: attolab

* PyMoDAQ is used as the core acquisition program of several experiments at CEMES/CNRS and the main
  interface of its HC-IUMI Ultrafast Electron Microscope
* The attolab platform at CEA Saclay started using it in 2019

Institutions using PyMoDAQ
**************************

|cemes| |attolab|


What they think of PyMoDAQ?
***************************

* *"The use of PyMoDAQ has really accelerated our experimental development by allowing to develop a modular acquisition
  system involving very different motorized stages or piezoactuators. It is now running everyday on our experiments,
  100% reliable"*, Dr Arnaud Arbouet, Senior Researcher CEMES/CNRS

* *Pymodaq is a python framework for data acquisition. If your specific device driver is not yet
  implemented, that is the only thing you will have to do. Pymodaq take care of the rest. Graphical
  user interface, synchronization of the instruments and so on, is already implemented. Once you have
  implemented your driver, you can release it for the community. That is how Pymodaq will get more and
  more complete. Of course you need to invest a bit of your time to get used to it, but it is worth it!*, Dr David
  Bresteau, Researcher at CEA Saclay, Attolab platform.

.. note::

  If you are using PyMoDAQ and would like to help to promote the project, please send your feedback to
  `sebastien.weber@cemes.fr <mailto:sebastien.weber@cemes.fr>`_ and we will include your message or logo on this page.
  If you wish to contribute, see :ref:`contributors`.


.. note::

  If you wish to communicate with users of PyMoDAQ, a mailing list exists:
  `pymodaq@services.cnrs.fr <mailto:pymodaq@services.cnrs.fr>`_

.. _convertingdata:

Converting raw data to be used in the retriever
***********************************************

Preamble
++++++++
If the non-linear trace was measured using another acquisition program than PyMoDAQ, it needs to be converted to
the proper format before being loaded into the retriever. PyMoDAQ-Femto uses a binary format known as hdf5__, as described in the PyMoDAQ documentation__.

__ https://www.hdfgroup.org/solutions/hdf5/
__ https://pymodaq.readthedocs.io/en/latest/usage/saving.html

We provide an example script to convert raw data, located inside the PyMoDAQ-Femto module, under :file:`pymodaq_femto/utils/convert_to_pymodaq_compatible.py`.

| If you are using **conda** with a dedicated environment as suggested in the :ref:`section_installation` section, this folder will be located inside the :file:`site-packages` folder, for instance in Windows something like:
| :file:`C:/Miniconda/envs/your_environment_name/Lib/site-packages/pymodaq_femto`
| This :file:`/utils` folder is also easily accessible on GitHub__.

If it is not already the case, raw data should be converted to numpy arrays. 5 arrays are needed:

* The 2D trace     [N x M numpy array]
* An array corresponding to the parameter axis (delay in Frog, glass insertion in Dscan, etc.) in physical units [N x 1 numpy array]
* An array with the wavelength axis of the trace [M x 1 numpy array]
* The fundamental spectrum (spectrum of light before non-linear conversion) [P x 1 numpy array]
* The wavelength axis of the fundamental spectrum  [P x 1 numpy array]

.. note::
    The retriever has an option to rescale any input array, so parameter or wavelength axes can be saved in any units. That being said, it is usually easier to save all data in standard units everytime (meters for wavelengths and insertions, seconds for delays, etc.).

The fundamental spectrum doesn't need to be on the same wavelength axis as the 2D trace, they will get interpolated on a common axis during retrieval.
The role of the fundamental is to compare retrieved spectrum with measured one (a good measure of the quality of retriever), and can also be used as an initial guess for the algorithm. But if you don't have one for every trace, just use any spectrum you have and the algorithm will still work.

__ https://github.com/CEMES-CNRS/pymodaq_femto/tree/main/src/pymodaq_femto/utils

Conversion
++++++++++
Once the 5 numpy arrays are loaded, you can use the utility functions of :file:`pymodaq_femto/utils/convert_to_pymodaq_compatible.py` to create a new .h5 file
and add all data to it, with the proper structure. This is copied here for convenience::

    # Open file and create scan node
    saver.init_file(addhoc_file_path=str(pathToSave.joinpath(fileName)), update_h5=True)
    scannode = saver.add_scan_group()
    scannode.set_attr('scan_type', "Scan1D")

    # Add all data
    saver.add_exp_insertion(scannode, parameter_axis)
    saver.add_exp_trace(scannode, trace_data, spectrum_trace_axis)
    saver.add_exp_fundamental(scannode, spectrum_fundamental_intensity, spectrum_fundamental_axis_wavelength)
    saver.close_file()

In this example, the 5 arrays are ``trace_data``, ``parameter_axis``, ``spectrum_trace_axis``, ``spectrum_fundamental_intensity`` and ``spectrum_fundamental_axis_wavelength``.

One raw DScan measurement is provided in :file:`pymodaq_femto/utils/raw_scans/example_measured_dscan_to_convert.h5`, it contains the 5 numpy arrays casted in the proper format.
It allows to run :file:`convert_to_pymodaq_compatible.py` directly. If all proper dependencies are installed, the script should save a converted file into :file:`pymodaq_femto/utils/converted_scans/`
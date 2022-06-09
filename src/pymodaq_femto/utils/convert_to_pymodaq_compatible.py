# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 18:30:59 2021

@author: Romain Géneaux
"""

import pypret
from pymodaq.daq_utils.h5modules import H5SaverBase
import os


class DScanCustomSaver(H5SaverBase):
    def add_exp_trace(self, node, trace, wl, label="Wavelength", units="m"):
        traceToSave = dict(data=trace, x_axis=dict(data=wl, label=label, units=units))
        det = self.get_set_group(node, "DSCAN_Measurement", title='DScan')
        self.set_attr(det, 'type', 'detector')
        # det = self.add_det_group(node, title='DScan')
        self.add_data(det, traceToSave, title='DScan trace')

    def add_exp_insertion(self, node, insertion, label="Insertion", units="m"):
        ax = self.add_navigation_axis(insertion, node, axis='x_axis', title=label)
        self.set_attr(ax, "label", label)
        self.set_attr(ax, "units", units)
        self.set_attr(ax, "nav_index", 0)

    def add_exp_fundamental(self, node, spectrum, wl, label="Wavelength", units="m"):
        spectrumToSave = dict(data=spectrum, x_axis=dict(data=wl, label=label, units=units))
        # det = self.add_det_group(node, title='Fundamental')
        det = self.get_set_group(node, "Fundamental_spectrum", title='Fundamental')
        self.add_data(det, spectrumToSave, title='Fundamental spectrum')


if __name__ == '__main__':
    from pathlib import Path
    saver = DScanCustomSaver()

    # # One particular implementation of saved data
    path_root = Path(__file__).parent.parent.parent.parent
    pathToLoad = path_root.joinpath("data/dscan_bank/")
    pathToSave = path_root.joinpath("data/dscan_bank_pymodaq/")
    if not pathToSave.is_dir():
        pathToSave.mkdir()
    fileIndex = 6
    # Load data
    fileName = "measured_dscan_" + str(fileIndex) + ".h5"
    pulse_from_spectrum, retrieved_pulse, measured_dscan, retrieved_dscan, raw_spectrum, folder = pypret.load(
        str(pathToLoad.joinpath(fileName)))
    insertion = measured_dscan.axes[0]
    spectrum_trace_axis = measured_dscan.axes[1]
    spectrum_fundamental_intensity = raw_spectrum.intensity
    spectrum_fundamental_axis_wavelength = raw_spectrum.wl
    trace_data = measured_dscan.data

    # # template to use:
    # trace_data = my_experimental_data_trace_as_2D_numpy_array
    # insertion = my_experimental_insertion_as_numpy_array
    # spectrum_trace_axis = my_experimental_trace_wavelength_axis_as_numpy_array
    # spectrum_fundamental_intensity = my_experimental_fundamental_spectrum_as_numpy_array
    # spectrum_fundamental_axis_wavelength = my_experimental_fundamental_spectrum_wavelength_axis_as_numpy_array

    # Open file and create scan node
    saver.init_file(addhoc_file_path=str(pathToSave.joinpath(fileName)), update_h5=True)
    scannode = saver.add_scan_group()
    scannode.set_attr('scan_type', "Scan1D")

    # Add all data
    saver.add_exp_insertion(scannode, insertion)
    saver.add_exp_trace(scannode, trace_data, spectrum_trace_axis)
    saver.add_exp_fundamental(scannode, raw_spectrum.intensity, spectrum_fundamental_axis_wavelength)
    saver.close_file()
    pass


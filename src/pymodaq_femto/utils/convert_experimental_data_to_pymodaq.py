# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 18:30:59 2021

@author: Romain Géneaux

Warning: ONLY WORKS IN DEBUGGER MODE FOR NOW
SET A BREAKPOINT AT LINE 53 (saver = DScanCustomSaver())
AND CONTINUE EXECUTION
"""

import pypret
from pymodaq.daq_utils.h5modules import H5Saver
import os


class DScanCustomSaver(H5Saver):
    def add_exp_trace(self, node, trace, wl, label="Wavelength", units="m"):
        traceToSave = dict(data=trace, x_axis=dict(data=wl, label=label, units=units))
        det = self.get_set_group(node, "Trace_Measurement", title='Trace')
        self.set_attr(det, 'type', 'detector')
        self.add_data(det, traceToSave, title='2D Trace')

    def add_exp_parameter(self, node, insertion, label="Insertion", units="m"):
        ax = self.add_navigation_axis(insertion, node, axis='x_axis', title=label)
        self.set_attr(ax, "label", label)
        self.set_attr(ax, "units", units)
        self.set_attr(ax, "nav_index", 0)

    def add_exp_fundamental(self, node, spectrum, wl, label="Wavelength", units="m"):
        spectrumToSave = dict(data=spectrum, x_axis=dict(data=wl, label=label, units=units))
        det = self.get_set_group(node, "Fundamental_spectrum", title='Fundamental')
        self.add_data(det, spectrumToSave, title='Fundamental spectrum')

    def convert_file(self, new_filepath, trace, parameter, trace_wl, fundamental, fundamental_wl,
                     parameter_label="Insertion", parameter_unit ="m",
                     fundamental_wl_label="Wavelength", fundamental_wl_unit="m"):

        self.init_file(addhoc_file_path=new_filepath, update_h5=True)

        scannode = self.add_scan_group()
        scannode.set_attr('scan_type', "Scan1D")

        self.add_exp_parameter(scannode, parameter, label=parameter_label, units=parameter_unit)
        self.add_exp_trace(scannode, trace, trace_wl)
        self.add_exp_fundamental(scannode, fundamental, fundamental_wl, label=fundamental_wl_label, units=fundamental_wl_unit)

        self.h5file.close()

pathToLoad = "./../../../tests/dscan_bank/"
pathToSave = "./../../../tests/dscan_convert/"

saver = DScanCustomSaver()

for fileIndex in range(18):
    # Original file
    fileName = "measured_dscan_" + str(fileIndex) + ".h5"
    new_filepath = pathToSave + fileName

    pulse_from_spectrum, retrieved_pulse, measured_dscan, retrieved_dscan, raw_spectrum, folder = pypret.load(
        pathToLoad + fileName)

    # Data to save
    trace = measured_dscan.data
    parameter = measured_dscan.axes[0]
    trace_wl = measured_dscan.axes[1]

    fundamental = raw_spectrum.intensity
    fundamental_wl = raw_spectrum.wl

    #Save
    saver.convert_file(new_filepath, trace, parameter, trace_wl, fundamental, fundamental_wl)
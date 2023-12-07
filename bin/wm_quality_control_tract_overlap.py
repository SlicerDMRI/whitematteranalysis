#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import os
import time
import warnings

import numpy as np
import vtk

import whitematteranalysis as wma
from whitematteranalysis.utils.opt_pckg import optional_package

matplotlib, have_mpl, _ = optional_package("matplotlib")
plt, _, _ = optional_package("matplotlib.pyplot")

if not have_mpl:
    warnings.warn(matplotlib._msg)
    warnings.warn("Cannot plot quality control data.")


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Perform quality control to render the overlap of two tracts.",
        epilog="Written by Fan Zhang (fzhang@bwh.harvard.edu).")
    parser.add_argument(
        'inputTract1',
        help='Input tract 1 as vtkPolyData (.vtk or .vtp).')
    parser.add_argument(
        'inputTract2',
        help='Input tract 2 as vtkPolyData (.vtk or .vtp).')
    parser.add_argument(
        'outputDirectory',
        help='Quality control information will be stored in the output directory, which will be created if it does not exist.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    print(f"<{os.path.basename(__file__)}> Starting...")
    
    if not os.path.exists(args.inputTract1):
        print(f"<{os.path.basename(__file__)}> Error: Input tract 1 {args.inputTract1} does not exist.")
        exit()
    
    if not os.path.exists(args.inputTract2):
        print(f"<{os.path.basename(__file__)}> Error: Input tract 2 {args.inputTract2} does not exist.")
        exit()
    
    output_dir = args.outputDirectory
    if not os.path.exists(output_dir):
        print(f"<{os.path.basename(__file__)}> Output directory {output_dir} does not exist, creating it.")
        os.makedirs(output_dir)
    
    input_polydatas = []
    input_polydatas.append(args.inputTract1)
    input_polydatas.append(args.inputTract2)
    
    number_of_subjects = len(input_polydatas)
    
    if have_mpl:
        plt.figure(1)
    
    # Loop over subjects and check each
    subject_idx = 1
    appender = vtk.vtkAppendPolyData()
    for fname in input_polydatas:
        subject_id = os.path.splitext(os.path.basename(fname))[0]
        print(f"Subject {subject_idx} / {number_of_subjects} ID: {subject_id}")
    
        # Read data
        pd = wma.io.read_polydata(fname)
    
        pd2, lengths, step_size = wma.filter.preprocess(pd, 20, return_lengths=True, verbose=False)
    
        number_rendered_fibers = 800
        pd3 = wma.filter.downsample(pd2, number_rendered_fibers, verbose=False)
        mask = np.ones(number_rendered_fibers)
        colors = np.multiply(mask, subject_idx)
        pd3 = wma.filter.mask(pd3, mask, colors, verbose=False)
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            appender.AddInputData(pd3)
        else:
            appender.AddInput(pd3)
    
        # release data
        del pd
        del pd2
        del pd3
        subject_idx += 1
    
    print(f"<{os.path.basename(__file__)}> Final step: rendering two vtk files together.")
    appender.Update()
    pd_all = appender.GetOutput()
    ren = wma.render.render(pd_all, colormap='hot', verbose=False)
    ren.save_views(output_dir, "tract_overlap")
    del ren
    del appender

if __name__ == '__main__':
    main()

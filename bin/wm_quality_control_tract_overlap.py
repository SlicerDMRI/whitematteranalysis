#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Run registration on the test dataset.

import argparse
import os
import numpy
import vtk
import time

try:
    import whitematteranalysis as wma
except:
    print(f"<{os.path.basename(__file__)}> Error importing white matter analysis package\n")
    raise

HAVE_PLT = 1

try:
    import matplotlib.pyplot as plt
except:
    print(f"<{os.path.basename(__file__)}> Error importing matplotlib.pyplot package, can't plot quality control data.\n")
    HAVE_PLT = 0    

def main():
    #-----------------
    # Parse arguments
    #-----------------
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
     
    args = parser.parse_args()
    
    print(f"<{os.path.basename(__file__)}> Starting...")
    
    if not os.path.exists(args.inputTract1):
        print(f"<{os.path.basename(__file__)}> Error: Input tract 1", args.inputTract1, "does not exist.")
        exit()
    
    if not os.path.exists(args.inputTract2):
        print(f"<{os.path.basename(__file__)}> Error: Input tract 2", args.inputTract2, "does not exist.")
        exit()
    
    output_dir = args.outputDirectory
    if not os.path.exists(output_dir):
        print(f"<{os.path.basename(__file__)}> Output directory", output_dir, "does not exist, creating it.")
        os.makedirs(output_dir)
    
    input_polydatas = []
    input_polydatas.append(args.inputTract1)
    input_polydatas.append(args.inputTract2)
    
    number_of_subjects = len(input_polydatas)
    
    if HAVE_PLT:
        plt.figure(1)
    
    # Loop over subjects and check each
    subject_idx = 1
    appender = vtk.vtkAppendPolyData()
    for fname in input_polydatas:
        subject_id = os.path.splitext(os.path.basename(fname))[0]
        print("Subject ", subject_idx, "/", number_of_subjects, "ID:", subject_id)
    
        # Read data
        pd = wma.io.read_polydata(fname)
    
        pd2, lengths, step_size = wma.filter.preprocess(pd, 20, return_lengths=True, verbose=False)
    
        number_rendered_fibers = 800
        pd3 = wma.filter.downsample(pd2, number_rendered_fibers, verbose=False)
        mask = numpy.ones(number_rendered_fibers)
        colors = numpy.multiply(mask, subject_idx)
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

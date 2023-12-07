#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import time

import numpy as np
import vtk

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Compare two tractography files to see how much registration has changed the points. The files must have the exact same size and lines and points, e.g. the original file and the file that results from applying a transform to a file.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.1 (2007): 1562-1575.\"")
    parser.add_argument(
        'inputFile1',
        help='A file of whole-brain tractography as vtkPolyData (.vtk or .vtp).')
    parser.add_argument(
        'inputFile2',
        help='A file of whole-brain tractography as vtkPolyData (.vtk or .vtp).')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    pd1 = wma.io.read_polydata(args.inputFile1)
    pd2 = wma.io.read_polydata(args.inputFile2)
    
    # Go through the whole polydata and compute the differences
    # Estimate this by converting to arrays and computing distances
    # figure out hemisphere labels for this cluster
    farray1 = wma.fibers.FiberArray()
    farray1.convert_from_polydata(pd1, points_per_fiber=30)
    farray2 = wma.fibers.FiberArray()
    farray2.convert_from_polydata(pd2, points_per_fiber=30)
    
    # compute the distances between the matching points
    ddx = farray1.fiber_array_r - farray2.fiber_array_r
    ddy = farray1.fiber_array_a - farray2.fiber_array_a
    ddz = farray1.fiber_array_s - farray2.fiber_array_s
    
    dx = np.square(ddx)
    dy = np.square(ddy)
    dz = np.square(ddz)
    
    distance = np.sqrt(dx + dy + dz)
    
    # print information
    print("MIN\tMAX\tMEAN\tSTD\tMEDIAN")
    print(np.min(distance), '\t', np.max(distance), '\t', np.mean(distance), '\t', np.std(distance), '\t', np.median(distance))
     
            
    # figure out what percentage is less than 1mm, 2mm, 3mm, 4mm, 5mm etc.
    sz = farray1.fiber_array_r.shape
    number_of_points = float(sz[0]*sz[1])

    print(f'{100 * np.sum(distance < 0.001) / number_of_points} percent < 0.001 mm')
    print(f'{100 * np.sum(distance < 0.01) / number_of_points} percent < 0.01 mm')
    print(f'{100 * np.sum(distance < 0.1) / number_of_points} percent < 0.1 mm')
    print(f'{100 * np.sum(distance < 0.5) / number_of_points} percent < 0.5 mm')
    print(f'{100 * np.sum(distance < 1) / number_of_points} percent < 1 mm')
    print(f'{100 * np.sum(distance < 1.2) / number_of_points} percent < 1.2 mm')
    print(f'{100 * np.sum(distance < 1.3) / number_of_points} percent < 1.3 mm')
    print(f'{100 * np.sum(distance < 1.4) / number_of_points} percent < 1.4 mm')
    print(f'{100 * np.sum(distance < 1.5) / number_of_points} percent < 1.5 mm')
    print(f'{100 * np.sum(distance < 2) / number_of_points} percent < 2 mm')
    print(f'{100 * np.sum(distance < 3) / number_of_points} percent < 3 mm')
    print(f'{100 * np.sum(distance < 4) / number_of_points} percent < 4 mm')
    print(f'{100 * np.sum(distance < 5) / number_of_points} percent < 5 mm')
    
    # output a new polydata showing what the distances are in case there is a pattern

if __name__ == '__main__':
    main()

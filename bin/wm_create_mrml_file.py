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
        description="Create a MRML file that includes all vtk and vtp files in the input directory. Color tracts with different colors in the scene.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
    parser.add_argument(
        'inputDirectory',
        help='A directory of whole-brain tractography as vtkPolyData (.vtk or .vtp). Output MRML scene file scene.mrml will be stored here and will point to all tractography files.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputDirectory):
        print(f"<{os.path.basename(__file__)}> Error: Input directory {args.inputDirectory} does not exist.")
        exit()
    
    mrml_filename = "scene.mrml"
    
    input_polydatas = wma.io.list_vtk_files(args.inputDirectory)
    number_of_files = len(input_polydatas)
    print(f"<{os.path.basename(__file__)}> Found {number_of_files} vtk files in input directory {args.inputDirectory}")
    
    # define R, G, B colors
    # hack a colormap. 0..255 values for each
    step = int(100*255.0 / (number_of_files-1))
    print(step, number_of_files)
    R = np.array(list(range(0,100*255+1, step))) / 100.0
    G = np.abs(list(range(100*-127,100*128+1, step)))* 2.0 / 100.0
    B = np.array(list(range(100*255+1,0, -step))) / 100.0
    
    #print len(R), len (G), len(B)
    #print R
    #print G
    #print B
    
    colors = list()
    idx = 0
    for pd in input_polydatas:
        colors.append([R[idx], G[idx],B[idx]])
        idx += 1
    colors = np.array(colors)
    print(colors)
    wma.mrml.write(input_polydatas, colors, os.path.join(args.inputDirectory, mrml_filename), ratio=1.0)

if __name__ == '__main__':
    main()

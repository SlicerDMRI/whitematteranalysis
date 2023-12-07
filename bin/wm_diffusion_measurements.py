#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
from warnings import warn

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Compute diffusion scalar measurements (such as FA, MD, etc). This script reports the mean statistics of each fiber cluster (or fiber tract) within the input folder.",
        epilog="Written by Fan Zhang (fzhang@bwh.harvard.edu)")
    parser.add_argument(
        'inputDirectory',
        help='Directory of fiber clustering results obtained by <wm_cluster_from_atlas.py> of multiple subjects. Make sure only the fiber clustering results are stored in this folder, making one subdirectory corresponding to one subject.')
    parser.add_argument(
        'outputCSV',
        help='Directory of output CSV files of fiber scalar measurement (computed using Slicer FiberTractMeasurements module).')
    parser.add_argument(
        'Slicer',
        help='Path of 3D Slicer.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputDirectory):
        print(f"Error: Input directory {args.inputDirectory} does not exist.")
        exit()
    
    cli= args.Slicer.split()[-1]
    if not os.path.exists(cli):
        warn(f"{cli} module could not be found, program may fail.")
    
    
    module_FTSM = args.Slicer + ' '
    
    outdir = os.path.split(args.outputCSV)[0]
    if not os.path.exists(outdir):
        print(f"Output directory {outdir} does not exist, creating it.")
        os.makedirs(outdir)
    
    print(f"<{os.path.basename(__file__)}>. Starting scalar measurement extraction.")
    print("")
    print(f"=====input directory======\n {args.inputDirectory}")
    print(f"=====output directory=====\n {outdir}")
    print(f"=====3D Slicer====\n {args.Slicer}")
    print("==========================")
    
    os.system(
        f"{module_FTSM} --inputtype Fibers_File_Folder --format Column_Hierarchy --separator Comma  --inputdirectory {args.inputDirectory} --outputfile {args.outputCSV}")
    
    print(f"<{os.path.basename(__file__)}> Measurements done at: {args.outputCSV}")

if __name__ == '__main__':
    main()

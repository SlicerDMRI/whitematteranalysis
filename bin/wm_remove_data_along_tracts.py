#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import multiprocessing
import os

import numpy as np
from joblib import Parallel, delayed

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Applies preprocessing to input directory. Downsamples, removes short fibers. Preserves tensors and scalar point data along retained fibers.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu")
    parser.add_argument(
        'inputDirectory',
        help='Contains whole-brain tractography as vtkPolyData file(s).')
    parser.add_argument(
        'outputDirectory',
        help='The output directory should be a new empty directory. It will be created if needed.')
    ## parser.add_argument(
    ##     '-f', action="store", dest="numberOfFibers", type=int,
    ##     help='Number of fibers to keep from each dataset.')
    ## parser.add_argument(
    ##     '-l', action="store", dest="fiberLength", type=int,
    ##     help='Minimum length (in mm) of fibers to keep.')
    parser.add_argument(
        '-j', action="store", dest="numberOfJobs", type=int,
        help='Number of processors to use.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputDirectory):
        print(f"Error: Input directory {args.inputDirectory} does not exist.")
        exit()
    
    outdir = args.outputDirectory
    if not os.path.exists(outdir):
        print(f"Output directory {outdir} does not exist, creating it.")
        os.makedirs(outdir)
    
    print("")
    print(f"=====input directory======\n {args.inputDirectory}")
    print(f"=====output directory=====\n {args.outputDirectory}")
    print("==========================")
    
    print(f'CPUs detected: {multiprocessing.cpu_count()}')
    if args.numberOfJobs is not None:
        parallel_jobs = args.numberOfJobs
    else:
        parallel_jobs = multiprocessing.cpu_count()
    print(f'Using N jobs: {parallel_jobs}')
    
    
    print("==========================")
    
    # Loop over input DWIs
    inputPolyDatas = wma.io.list_vtk_files(args.inputDirectory)
    
    print(f"<{os.path.basename(__file__)}> Input number of files: f{len(inputPolyDatas)}")
    
    # for testing
    #inputPolyDatas = inputPolyDatas[0:2]
    
    def pipeline(inputPolyDatas, sidx, args):
        # get subject identifier from unique input filename
        # -------------------
        #subjectID = os.path.splitext(os.path.basename(inputPolyDatas[sidx]))[0]
        fname = os.path.basename(inputPolyDatas[sidx])
        print(f"<{os.path.basename(__file__)}> {sidx + 1} / {len(inputPolyDatas)}")
    
        # read input vtk data
        # -------------------
        inpd = wma.io.read_polydata(inputPolyDatas[sidx])
        num_lines = inpd.GetNumberOfLines()
        fiber_mask = np.ones(num_lines)
        outpd = wma.filter.mask(inpd, fiber_mask, preserve_point_data=False, preserve_cell_data=False, verbose=False)
        print(f"Number of fibers retained: {outpd.GetNumberOfLines()} / {num_lines}")
            
        # outputs
        # -------------------
        fname = os.path.join(args.outputDirectory, fname)
        try:
            print(f"Writing output polydata {fname}...")
            wma.io.write_polydata(outpd, fname)
        except:
            print("Unknown exception in IO")
            raise
        del outpd
        del inpd
    
    # loop over all inputs
    Parallel(n_jobs=parallel_jobs, verbose=0)(
            delayed(pipeline)(inputPolyDatas, sidx, args)
            for sidx in range(0, len(inputPolyDatas)))
    
    print("Launched all jobs")
    exit()

if __name__ == '__main__':
    main()

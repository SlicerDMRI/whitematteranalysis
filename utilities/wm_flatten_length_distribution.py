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
        description="Samples fibers such that the output distribution of fiber lengths is approximately flat, within the input length range.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu")
    parser.add_argument(
        'inputDirectory',
        help='Contains whole-brain tractography as vtkPolyData file(s).')
    parser.add_argument(
        'outputDirectory',
        help='The output directory should be a new empty directory. It will be created if needed.')
    parser.add_argument(
        '-f', action="store", dest="numberOfFibers", type=int,
        help='Number of fibers to keep from each dataset.')
    parser.add_argument(
        '-lmin', action="store", dest="fiberLengthMin", type=int,
        help='Minimum length (in mm) of fibers to keep.')
    parser.add_argument(
        '-lmax', action="store", dest="fiberLengthMax", type=int,
        help='Maximum length (in mm) of fibers to keep.')
    parser.add_argument(
        '-nbins', action="store", dest="numberOfBins", type=int, default=10,
        help='Number of bins for sampling the fiber length distribution.')
    parser.add_argument(
        '-nf', action="store", dest="fibersPerBin", type=int, default=1000,
        help='Number of fibers sampled per bin. If you have specified the total number of fibers, that will override this parameter.')
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

    print(f"{os.path.basename(__file__)}. Starting streamline length distribution flattening.")
    print("")
    print(f"=====input directory======\n {args.inputDirectory}")
    print(f"=====output directory=====\n {args.outputDirectory}")
    print("==========================")

    if args.numberOfFibers is not None:
        print(f"fibers to retain per subject: {args.numberOfFibers}")
        args.fibersPerBin = np.divide(args.numberOfFibers,args.numberOfBins)
    else:
        print("fibers to retain per subject: ALL")

    if args.fiberLengthMin is not None:
        print(f"minimum length of fibers to retain (in mm): {args.fiberLengthMin}")
    else:
        print(f"minimum length of fibers to retain (in mm): 0")

    if args.fiberLengthMax is not None:
        print("maximum length of fibers to retain (in mm): {args.fiberLengthMax}")

    print(f"Bins: {args.numberOfBins}")
    print(f"Fibers per bin: {args.fibersPerBin}")

    if args.numberOfJobs is not None:
        parallel_jobs = args.numberOfJobs
    else:
        parallel_jobs = 1
    print(f'Using N jobs: {parallel_jobs}')


    print("==========================")

    # Loop over input DWIs
    inputMask1 = f"{args.inputDirectory}/*.vtk"
    inputMask2 = f"{args.inputDirectory}/*.vtp"

    inputPolyDatas = glob.glob(inputMask1) + glob.glob(inputMask2)

    print(f"<{os.path.basename(__file__)}> Input number of files: {len(inputPolyDatas)}")

    # for testing
    #inputPolyDatas = inputPolyDatas[0:2]

    def pipeline(inputPolyDatas, sidx, args):
        # get subject identifier from unique input filename
        # -------------------
        subjectID = os.path.splitext(os.path.basename(inputPolyDatas[sidx]))[0]
        id_msg = f"<{os.path.basename(__file__)}> {sidx + 1} / {len(inputPolyDatas)}"
        msg = f"**Starting subject: {subjectID}"
        print(id_msg + msg)

        # read input vtk data
        # -------------------
        msg = f"**Reading input: {subjectID}"
        print(id_msg + msg)

        wm = wma.io.read_polydata(inputPolyDatas[sidx])

        num_lines = wm.GetNumberOfLines()
        #print "Input number of fibers", num_lines

        wm2 = wma.filter.flatten_length_distribution(wm, args.fiberLengthMin, args.fiberLengthMax, args.numberOfBins, args.fibersPerBin, verbose=True)

        # outputs
        # -------------------
        msg = f"**Writing output data for subject: {subjectID}"
        print(id_msg, msg)

        fname = os.path.join(args.outputDirectory, f'{subjectID}_flat.vtp')
        try:
            print(f"Writing output polydata {fname}...")
            wma.io.write_polydata(wm2, fname)
            print(f"Wrote output {fname}.")
        except:
            print("Unknown exception in IO")
            raise
        del wm2

    # loop over all inputs
    Parallel(n_jobs=parallel_jobs, verbose=0)(
            delayed(pipeline)(inputPolyDatas, sidx, args)
            for sidx in range(0, len(inputPolyDatas)))

    exit()


if __name__ == "__main__":
    main()

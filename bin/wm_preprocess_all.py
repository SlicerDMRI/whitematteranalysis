#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import multiprocessing
import os
from pathlib import Path

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
    parser.add_argument(
        '-f', action="store", dest="numberOfFibers", type=int,
        help='Number of fibers to keep from each dataset.')
    parser.add_argument(
        '-l', action="store", dest="fiberLength", type=int,
        help='Minimum length (in mm) of fibers to keep.')
    parser.add_argument(
        '-lmax', action="store", dest="maxFiberLength", type=int,
        help='Maximum length (in mm) of fibers to keep.')
    parser.add_argument(
        '-j', action="store", dest="numberOfJobs", type=int,
        help='Number of processors to use.')
    parser.add_argument(
        '-retaindata', action='store_true', dest="flag_retaindata",
        help='If given, all point and cell data stored along the tractography will be retained.')
    parser.add_argument(
        '--nonidentical', action='store_true',
        help='Obtain nonidentical results across runs for downsampling.')

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

    print(f"{os.path.basename(__file__)}. Starting white matter preprocessing.")
    print("")
    print(f"=====input directory======\n {args.inputDirectory}")
    print(f"=====output directory=====\n {args.outputDirectory}")
    print("==========================")

    if args.numberOfFibers is not None:
        print(f"fibers to retain per subject: {args.numberOfFibers}")
    else:
        print("fibers to retain per subject: ALL")

    if args.fiberLength is not None:
        print(f"minimum length of fibers to retain (in mm): {args.fiberLength}")
    else:
        print("minimum length of fibers to retain (in mm): 0")

    print(f'CPUs detected: {multiprocessing.cpu_count()}')
    if args.numberOfJobs is not None:
        parallel_jobs = args.numberOfJobs
    else:
        parallel_jobs = multiprocessing.cpu_count()
    print(f'Using N jobs: {parallel_jobs}')
    
    if args.flag_retaindata:
        print("Retain all data stored along the tractography.")
    else:
        print("Remove all data stored along the tractography and only keep fiber streamlines.")
    retaindata = args.flag_retaindata

    random_seed = 1234
    if args.nonidentical:
        random_seed = None

    print("==========================")

    # Loop over input DWIs
    inputPolyDatas = wma.io.list_vtk_files(args.inputDirectory)

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

        # remove short fibers
        # -------------------
        wm2 = None
        if args.fiberLength is not None or args.maxFiberLength is not None:
            msg = f"**Preprocessing: {subjectID}"
            print(id_msg + msg)

            if args.fiberLength is not None:
                minlen = args.fiberLength
            else:
                minlen = 0

            if args.maxFiberLength is not None:
                maxlen = args.maxFiberLength
            else:
                maxlen = None

            wm2 = wma.filter.preprocess(wm, minlen, preserve_point_data=retaindata, preserve_cell_data=retaindata, verbose=False, max_length_mm=maxlen)
            print(f"Number of fibers retained (length threshold {args.fiberLength}): {wm2.GetNumberOfLines()} / {num_lines}")

        if wm2 is None:
            wm2 = wm
        else:
            del wm

        # downsample
        # -------------------
        wm3 = None
        if args.numberOfFibers is not None:
            msg = f"**Downsampling input: {subjectID} number of fibers: {args.numberOfFibers}"
            print(id_msg + msg)

            # , preserve_point_data=True needs editing of preprocess function to use mask function
            wm3 = wma.filter.downsample(wm2, args.numberOfFibers, preserve_point_data=retaindata, preserve_cell_data=retaindata, verbose=False, random_seed=random_seed)
            print(f"Number of fibers retained: {wm3.GetNumberOfLines()} / {num_lines}")

        if wm3 is None:
            wm3 = wm2
        else:
            del wm2

        # outputs
        # -------------------
        msg = f"**Writing output data for subject: {subjectID}"
        print(id_msg, msg)

        ext = Path(inputPolyDatas[sidx]).suffixes[0][1:]
        fname = os.path.join(args.outputDirectory, f"{subjectID}_pp.{ext}")

        try:
            print(f"Writing output polydata {fname}...")
            wma.io.write_polydata(wm3, fname)
            print(f"Wrote output {fname}.")
        except:
            print("Unknown exception in IO")
            raise
        del wm3

    # loop over all inputs
    Parallel(n_jobs=parallel_jobs, verbose=0)(
            delayed(pipeline)(inputPolyDatas, sidx, args)
            for sidx in range(0, len(inputPolyDatas)))

    exit()

if __name__ == '__main__':
    main()

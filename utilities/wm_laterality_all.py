#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import os

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Applies white matter laterality pipeline to input directory.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu")
    parser.add_argument(
        'inputDirectory',
        help='Contains whole-brain tractography as vtkPolyData file(s).')
    parser.add_argument(
        'outputDirectory',
        help='The output directory should be a new empty directory.')
    parser.add_argument(
        '-f', action="store", dest="numberOfFibers", type=int,
        help='Number of fibers to analyze from each dataset.')
    parser.add_argument(
        '-l', action="store", dest="fiberLength", type=int,
        help='Minimum length (in mm) of fibers to analyze.')
    parser.add_argument(
        '-s', action="store", dest="sigma", type=float,
        help='Sigma for laterality computation. Useful values are 10-50 (mm).')
    parser.add_argument(
        '-t', action="store", dest="threshold", type=float,
        help='Threshold lower fiber distances to 0. Useful range 0-5mm.')
    parser.add_argument(
        '-rm_outlier', action='store_true', dest="flag_removeOutliers")
    parser.add_argument(
        '-equal_fiber_num', action='store_true', dest="flag_equalFibers",
        help='To analyze an equal number of fibers per hemisphere')
    parser.add_argument(
        '-fibers_per_hem', action="store", dest="numberOfFibersPerHem", type=int,
        help='Number of fibers to analyze from each hemisphere.')

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
        print(f"<{os.path.basename(__file__)}> Output directory {outdir} does not exist, creating it.")
        os.makedirs(outdir)

    print("f<{os.path.basename(__file__)}> Starting white matter laterality computation.")
    print("")
    print(f"=====input directory======\n {args.inputDirectory}")
    print(f"=====output directory=====\n {args.outputDirectory}")
    print("==========================")

    if args.numberOfFibers is not None:
        print("fibers to analyze per subject: ", args.numberOfFibers)
    else:
        print("fibers to analyze per subject: ALL")

    if args.fiberLength is not None:
        print("minimum length of fibers to analyze (in mm): ", args.fiberLength)

    if args.sigma is not None:
        print("sigma for laterality index computation: ", args.sigma)

    if args.threshold is not None:
        print("intra-fiber distance threshold (in mm): ", args.threshold)

    if args.flag_removeOutliers:
        print("Automatic outlier removal is ON.")
    else:
        print("Automatic outlier removal is OFF.")

    if args.flag_equalFibers:
        print("Use equal fiber number from each hemisphere is ON.")
    else:
        print("Use equal fiber number from each hemisphere is OFF. Using input fiber number.")

    if args.numberOfFibersPerHem is not None:
        print(f"fibers to analyze per hemisphere: {args.numberOfFibersPerHem}")
    else:
        print("fibers to analyze per hemisphere: all or equal")

    print("==========================")

    inputPolyDatas = wma.io.list_vtk_files(args.inputDirectory)

    print(f"<{os.path.basename(__file__)}> Input number of files: ", len(inputPolyDatas))

    # loop over all inputs
    for sidx in range(0, len(inputPolyDatas)):

        # get subject identifier from unique input filename
        # -------------------
        subjectID = os.path.splitext(os.path.basename(inputPolyDatas[sidx]))[0]
        id_msg = f"<{os.path.basename(__file__)}> {sidx + 1} / {len(inputPolyDatas)}"
        msg = f"**Starting subject: {subjectID}"
        print(id_msg, msg)

        # read input vtk data
        # -------------------
        msg = f"**Reading input: {subjectID}"
        print(id_msg, msg)

        wm = wma.io.read_polydata(inputPolyDatas[sidx])

        # remove short fibers
        # -------------------
        if args.fiberLength is not None:
            msg = f"**Preprocessing: {subjectID}"
            print(id_msg, msg)

            wm = wma.filter.preprocess(wm, args.fiberLength, remove_u=True, remove_u_endpoint_dist=50, remove_brainstem=True)
            print(f"Number of fibers retained: {wm.GetNumberOfLines()}")

        # remove outlier fibers
        # -------------------
        if args.flag_removeOutliers:
            msg = f"**Removing outliers: {subjectID}"
            print(id_msg, msg)

            # if it's huge downsample to twice requested size first
            if args.numberOfFibers is not None:
                if (wm.GetNumberOfLines() > args.numberOfFibers * 2):
                    wm = wma.filter.downsample(
                        wm, args.numberOfFibers * 2)
                    print(wm.GetNumberOfLines())

            outlierThreshold = 10
            wm = wma.filter.remove_outliers(wm, outlierThreshold)

        # downsample if requested
        # -------------------
        if args.numberOfFibers is not None:
            msg = f"**Downsampling input: {subjectID}"
            print(id_msg, msg)

            wm = wma.filter.downsample(wm, args.numberOfFibers)

        # midsagittal alignment is already done
        wm_align = wm

        # compute laterality on each dataset
        # -------------------
        msg = f"**Computing laterality: {subjectID}"
        print(id_msg, msg)

        laterality = wma.laterality.WhiteMatterLaterality()
        if args.sigma is not None:
            laterality.sigma = args.sigma
        if args.threshold is not None:
            laterality.threshold = float(args.threshold)
        else:
            laterality.threshold = 0.0
        if args.flag_equalFibers:
            laterality.use_equal_fibers = True
        else:
            laterality.use_equal_fibers = False

        if args.numberOfFibersPerHem is not None:
            laterality.fibers_per_hemisphere = args.numberOfFibersPerHem

        laterality_results = laterality.compute(wm_align)

        # outputs
        # -------------------
        msg = f"**Writing output data for subject: {subjectID}"
        print(id_msg, msg)

        outdir = os.path.join(args.outputDirectory, subjectID)
        try:
            #print "Writing output files..."
            laterality_results.write(outdir)
            print("wrote output")
        except:
            print("Unknown exception in IO")
            raise

    exit()


if __name__ == "__main__":
    main()

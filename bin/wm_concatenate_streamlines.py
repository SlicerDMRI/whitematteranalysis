#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Concatenate streamlines contained in different tractograms. Point and cell
scalar and tensor data are not preserved.
"""

import os
import argparse

import vtk

import whitematteranalysis as wma


def pipeline(fnames, out_fname, verbose=False):

    out_polydata = vtk.vtkPolyData()

    # Get subject identifier from unique input filename
    for index, fname in enumerate(fnames):
        sub_id = os.path.splitext(os.path.basename(fname))[0]
        id_msg = f"{os.path.basename(__file__)} {index + 1} / {len(fnames)}"
        msg = f"**Starting subject: {sub_id}"
        print(id_msg + msg)

        # Read tractogram
        msg = f"**Reading input: {sub_id}"
        print(id_msg + msg)

        polydata = wma.io.read_polydata(fname)
        print(f"Number of streamlines: {polydata.GetNumberOfLines()}")

        # Concatenate
        out_polydata = wma.filter.concatenate_streamlines(
            [out_polydata, polydata],
            _verbose=verbose
        )

    print(f"Number of streamlines concatenated: {out_polydata.GetNumberOfLines()}")

    # Output
    try:
        print(f"Writing output polydata {out_fname}...")
        wma.io.write_polydata(out_polydata, out_fname)
        print(f"Wrote output {out_fname}.")
    except:
        print("Unknown exception in IO")
        raise


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'inputDirectory',
        help='Directory containing tractography files (.*vtk|.*vtp).')
    parser.add_argument(
        'outputFilename',
        help='Output filename (.*vtk|.*vtp).')
    parser.add_argument(
        '-verbose', action='store_true', dest="flag_verbose",
        help='Verbose. Run with -verbose to print operation information.')

    return parser


def _parse_args(parser):

    args = parser.parse_args()
    return args


def main():
    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputDirectory):
        print(f"Error: Input directory {args.inputDirectory} does not exist.")
        exit()

    out_fname = args.outputFilename
    if os.path.exists(out_fname):
        msg = f"Output file {out_fname} exists. Remove or rename the output file."
        parser.error(msg)

    print(f"{os.path.basename(__file__)}. Starting streamline concatenation.")
    print("")
    print("=====input directory======\n", args.inputDirectory)
    print("=====output filename=====\n", args.outputFilename)
    print("==========================")

    verbose = args.flag_verbose

    fnames = wma.io.list_vtk_files(args.inputDirectory)

    print(f"<{os.path.basename(__file__)}> Input number of files: ", len(fnames))

    pipeline(fnames, out_fname, verbose=verbose)

    exit()


if __name__ == "__main__":
    main()

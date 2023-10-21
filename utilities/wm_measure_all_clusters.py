#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import numpy as np

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Measures mean FA, etc. in tractography clusters in a directory.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
    parser.add_argument(
        'inputDirectory',
        help='A directory of tractography as vtkPolyData (.vtk or .vtp).')
    parser.add_argument(
        'outputFile',
        help='Output measurements will be recorded here.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputDirectory):
        print(f"Error: Input directory {args.inputDirectory} does not exist.")
        exit()

    print(f"<{os.path.basename(__file__)}> Starting computation.")
    print("")
    print(f"=====input directory ======\n {args.inputDirectory}")
    print(f"=====output file =====\n {args.outputFile}")
    print("==========================")
    print("")

    def compute_point_data_stats(pd, array_name):
        point_array = pd.GetPointData().GetArray(array_name)
        if point_array is None:
            return None
        # make sure this is a one-component scalar
        if point_array.GetNumberOfComponents() > 1:
            print(f"Error in compute_point_data_stats: Array {array_name} has more than one component {point_array.GetNumberOfComponents()}.")
            return None
        print(point_array)
        num_points = pd.GetNumberOfPoints()
        points_copy = np.zeros(num_points)

        for pidx in range(0, num_points):
            if (pidx % 1000) == 0:
                print(f"Point {pidx} / {num_points}")
            # this assumes we have scalars here
            points_copy[pidx] = point_array.GetTuple(pidx)[0]

        points_mean = np.mean(points_copy)
        #points_std = np.std(points_copy)
        #points_median = np.median(points_copy)

        print(f"Mean {array_name} : {points_mean}")

        return points_mean



    input_polydatas = wma.io.list_vtk_files(args.inputDirectory)
    number_of_clusters = len(input_polydatas)

    input_polydatas = input_polydatas[0:10]

    print(f"<{os.path.basename(__file__)}> Input number of vtk/vtp files: {number_of_clusters}")

    scalars = ['FA', 'Trace', 'FA1', 'FA2', 'Trace1', 'Trace2']

    output_rows = list()

    # read in data
    input_pds = list()
    for fname in input_polydatas:
        print(fname)
        # read data
        print(f"<{os.path.basename(__file__)}> Reading input file: {fname}")
        pd = wma.io.read_polydata(fname)
        # preprocessing step: minimum length
        print(f"<{os.path.basename(__file__)}> Computing stats for input file: {fname}")
        output_row = list()
        output_row.append(fname)
        for sc in scalars:
            stat = compute_point_data_stats(pd, sc)
            if stat is not None:
                output_row.append(sc)
                output_row.append(stat)
        output_rows.append(output_row)

    # output a csv file
    f = open(args.outputFile, 'w')

    for row in output_rows:
        outstr = ''
        for item in row:
            if outstr != '':
                outstr = f'{outstr}\t'
            outstr = f'{outstr}{str(item)}'
        outstr = f'{outstr}\n'
        f.write(outstr)

    f.close()


if __name__ == "__main__":
    main()

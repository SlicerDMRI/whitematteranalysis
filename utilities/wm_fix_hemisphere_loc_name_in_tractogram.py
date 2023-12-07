#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Fix hemisphere location property name in tractogram: older versions of WMA
contain a typo in the hemisphere location VTK array scalar property. Change
`HemisphereLocataion` to `HemisphereLocation`.
"""

import argparse
import os

import vtk

import whitematteranalysis as wma


def fix_hemisphere_location_property_name(polydata):

    out_polydata = vtk.vtkPolyData()
    point_data = polydata.GetPointData()

    if point_data.GetNumberOfArrays():
        point_data_array_indices = list(range(point_data.GetNumberOfArrays()))

        for idx in point_data_array_indices:
            array = point_data.GetArray(idx)

            if array.GetName() == "HemisphereLocataion":
                print("HemisphereLocataion is in the input data: re-writing to HemisphereLocation.")
            elif array.GetName() == "HemisphereLocation":
                print("HemisphereLocation is in the input data: no re-write needed.")
            else:
                print("Hemisphere location property is not in the input data: no re-write needed.")

    vtk_array = vtk.vtkDoubleArray()
    vtk_array.SetName("HemisphereLocation")

    polydata.GetLines().InitTraversal()

    for l_idx in range(polydata.GetNumberOfLines()):
        point_ids = vtk.vtkIdList()
        polydata.GetLines().GetNextCell(point_ids)

        for p_idx in range(point_ids.GetNumberOfIds()):
            vtk_array.InsertNextTuple1(l_idx)

    out_polydata.GetPointData().AddArray(vtk_array)

    return out_polydata


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('in_fname', help='Input filename (.*vtk|.*vtp).')
    parser.add_argument('out_fname', help='Output filename (.*vtk|.*vtp).')

    return parser


def _parse_args(parser):

    args = parser.parse_args()
    return args


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    out_fname = args.outputFilename
    if os.path.exists(out_fname):
        msg = f"Output file {out_fname} exists. Remove or rename the output file."
        parser.error(msg)

    print(f"{os.path.basename(__file__)}. Starting.")
    print("")
    print(f"=====input filename======\n {args.in_fname}")
    print(f"=====output filename=====\n {args.out_fname}")
    print("==========================")

    polydata = wma.io.read_polydata(args.in_fname)
    out_polydata = fix_hemisphere_location_property_name(polydata)
    wma.io.write_polydata(out_polydata, args.out_fname)


if __name__ == "__main__":
    main()

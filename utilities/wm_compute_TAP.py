#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import argparse
import os
import numpy
import vtk
import whitematteranalysis as wma


def main():
    parser = argparse.ArgumentParser(
        description="Compute Tract Anatomical Profile",
        epilog="Written by Fan Zhang, zhangfanmark@gmail.com.")
    parser.add_argument("-v", "--version",
        action="version", default=argparse.SUPPRESS,
        version='1.0',
        help="Show program's version number and exit")
    parser.add_argument(
        'inputDirectory',
        help='A directory containing subdirectories for all clustered subjects.')

    args = parser.parse_args()

    if not os.path.isdir(args.inputDirectory):
        print("Error: Input directory", args.inputDirectory, "does not exist or is not a directory.")
        exit()

    def list_vtk_files(input_dir):
        input_mask = f"{input_dir}/*.vtk"
        input_mask2 = f"{input_dir}/*.vtp"
        input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
        input_pd_fnames = sorted(input_pd_fnames)
        return(input_pd_fnames)

    def compute_TAP(inpd):

        inpointdata = inpd.GetPointData()
        NoS = inpd.GetNumberOfLines()

        TAP_str = ""
        if NoS != 0:
            point_data_array_indices = list(range(inpointdata.GetNumberOfArrays()))
            for idx in point_data_array_indices:
                array = inpointdata.GetArray(idx)

                if array.GetName() == 'region_label':
                    inpd.GetLines().InitTraversal()
                    fiber_regions = []
                    for lidx in range(0, inpd.GetNumberOfLines()):
                        ptids = vtk.vtkIdList()
                        inpd.GetLines().GetNextCell(ptids)

                        regions = numpy.zeros(ptids.GetNumberOfIds())
                        for pidx in range(0, ptids.GetNumberOfIds()):
                            regions[pidx] = array.GetTuple(ptids.GetId(pidx))[0]

                        regions = numpy.unique(regions)
                        fiber_regions.append(regions)

                    all_regions = numpy.unique(numpy.concatenate(fiber_regions))

                    r_prob = []
                    r_list = []
                    for r in all_regions:
                        count = 0
                        for fiber_region in fiber_regions:
                            if r in fiber_region:
                                count = count + 1

                        prob = count / NoS
                        r_prob.append(prob)
                        r_list.append(r)

                    r_prob = numpy.array(r_prob)
                    r_list = numpy.array(r_list)
                    arg = numpy.argsort(-r_prob)
                    r_prob = r_prob[arg]
                    r_list = r_list[arg]

                    for r, p in zip(r_list, r_prob):
                        TAP_str += "'%d:%0.6f'," % (r, p)

        return TAP_str


    vtk_files = list_vtk_files(args.inputDirectory)

    all_TAP_str = "Cluster Index,Region:Probability,\n"
    for vtk_file in vtk_files:
        print(vtk_file)
        tract_name = os.path.basename(vtk_file).replace(".vtp", "").replace(".vtk", "")
        pd = wma.io.read_polydata(vtk_file)
        TAP_str = compute_TAP(pd)
        all_TAP_str += tract_name + "," + TAP_str + "\n"

    print(all_TAP_str)

    f = open(os.path.join(args.inputDirectory, 'TAP.csv'), 'a')
    f.write(all_TAP_str)
    f.close()


if __name__ == "__main__":
    main()

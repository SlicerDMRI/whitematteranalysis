#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import os

import numpy as np

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Decide a cluster belonging to commissural or hemispheric in the atlas.",
        epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")
    parser.add_argument(
        'inputAtlasDirectory',
        help='Directory of the fiber clustering atlas. An output CSV file will be generated, showing the location of each fiber cluster: left-hemisphere, right-hemisphere, commissural, or Not Given if the fibers of a cluster are relatively even distributed in these three regions. ')
    parser.add_argument(
        '-pthresh', action="store", dest="hemispherePercentThreshold", type=float, default=0.6,
        help='(A same pthresh value should be used in <wm_separate_clusters_by_hemisphere.py>) The percent of a fiber that has to be in one hemisphere to consider the fiber as part of that hemisphere (rather than a commissural fiber). Default number is 0.6, where a higher number will tend to label fewer fibers as hemispheric and more fibers as commissural (not strictly in one hemisphere or the other), while a lower number will be stricter about what is classified as commissural.')
    parser.add_argument(
        '-advanced_times_threshold', action="store", dest="advanced_times_threshold", type=float, default=3,
        help='(Advanced parameter should be used according to applications.) For one cluster, if it has the hemispheric (left or right) fibers three times more than the commissural part, this cluster will be considered as a hemispheric cluster. In a similar way, one cluster needs to have the commissural fibers three times more than the hemispheric (left and right) fibers to be considered as a commissural cluster. Users can change -advanced_num_threshold to have more strict classification. For example, if advanced_num_threshold = 5, the cluster needs to have the commissural fibers five times more than the hemispheric (left and right) fibers to be considered as a commissural cluster.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputAtlasDirectory):
        print("Error: Input directory", args.inputTractDirectory, "does not exist.")
        exit()

    print(f"<{os.path.basename(__file__)}>. Starting processing.")
    print("")
    print(f"=====input atlas directory======\n {args.inputAtlasDirectory}")
    print(f"=====pthresh====\n {args.hemispherePercentThreshold}")
    print(f"=====advanced_times_threshold====\n {args.advanced_times_threshold}")

    def list_cluster_files(input_dir):
        # Find input files
        input_mask = f"{input_dir}/cluster_*.vtk"
        input_mask2 = f"{input_dir}/cluster_*.vtp"
        input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
        input_pd_fnames = sorted(input_pd_fnames)
        return(input_pd_fnames)

    input_polydatas = list_cluster_files(args.inputAtlasDirectory)
    number_of_clusters = len(input_polydatas)

    points_per_fiber = 40
    print("%20s:  %10s  %10s  %10s  -  %10s" \
              % ('cluster', 'left', 'right', 'comm', 'locations'))

    output_file = open(os.path.join(args.inputAtlasDirectory, 'clusters_location.txt'), 'w')
    outstr = 'cluster' + '\t'+ 'left hemispheric fibers' + '\t'+ 'right hemispheric fibers' + '\t'+ 'commissural fibers' + '\t'+ 'location' + '\n'

    hemi_list = []
    comm_list = []
    ng_list = []
    for fname in input_polydatas:

        fname_base = os.path.basename(fname)
        pd = wma.io.read_polydata(fname)

        fibers = wma.fibers.FiberArray()
        fibers.points_per_fiber = points_per_fiber
        fibers.hemisphere_percent_threshold = args.hemispherePercentThreshold
        # must request hemisphere computation from object
        fibers.hemispheres = True
        # Now convert to array with points and hemispheres as above
        fibers.convert_from_polydata(pd)

        num_left = len(fibers.index_left_hem)
        num_right = len(fibers.index_right_hem)
        num_comm = len(fibers.index_commissure)

        location = ''
        if num_comm == 0:
            location = 'hemispheric'
            hemi_list.append(fname_base)
        elif num_left > args.advanced_times_threshold * num_comm or num_right > args.advanced_times_threshold * num_comm:
            location = 'hemispheric'
            hemi_list.append(fname_base)
        elif num_comm > args.advanced_times_threshold * num_left and num_comm > args.advanced_times_threshold * num_right:
            location = 'commissural'
            comm_list.append(fname_base)
        else:
            location = 'Not Given'
            ng_list.append(fname_base)

        print(f"{fname_base:20s}:  {num_left:10s}  {num_right:10s}  {num_comm:10s}  -  {location:10s}")

        outstr = f'{outstr}{fname_base}\t{str(num_left)}\t{str(num_right)}\t{str(num_comm)}\t{location}\n'

    output_file.write(outstr)
    output_file.close()

    # hemisphere
    number_of_files = len(hemi_list)
    step = int(100 * 255.0 / (number_of_files - 1))
    R = np.array(list(range(0, 100 * 255 + 1, step))) / 100.0
    G = np.abs(list(range(100 * -127, 100 * 128 + 1, step))) * 2.0 / 100.0
    B = np.array(list(range(100 * 255 + 1, 0, -step))) / 100.0

    colors = list()
    idx = 0
    for pd in hemi_list:
        colors.append([R[idx], G[idx], B[idx]])
        idx += 1
    colors = np.array(colors)
    hemi_mrml_filename = os.path.join(args.inputAtlasDirectory, "clusters_location_hemispheric.mrml")
    wma.mrml.write(hemi_list, colors, hemi_mrml_filename, ratio=1.0)

    #commissural
    number_of_files = len(comm_list)
    step = int(100 * 255.0 / (number_of_files - 1))
    R = np.array(list(range(0, 100 * 255 + 1, step))) / 100.0
    G = np.abs(list(range(100 * -127, 100 * 128 + 1, step))) * 2.0 / 100.0
    B = np.array(list(range(100 * 255 + 1, 0, -step))) / 100.0

    colors = list()
    idx = 0
    for pd in comm_list:
        colors.append([R[idx], G[idx], B[idx]])
        idx += 1
    colors = np.array(colors)
    comm_mrml_filename = os.path.join(args.inputAtlasDirectory, "clusters_location_commissural.mrml")
    wma.mrml.write(comm_list, colors, comm_mrml_filename, ratio=1.0)

    #Not Given
    number_of_files = len(ng_list)
    step = int(100 * 255.0 / (number_of_files - 1))
    R = np.array(list(range(0, 100 * 255 + 1, step))) / 100.0
    G = np.abs(list(range(100 * -127, 100 * 128 + 1, step))) * 2.0 / 100.0
    B = np.array(list(range(100 * 255 + 1, 0, -step))) / 100.0

    colors = list()
    idx = 0
    for pd in ng_list:
        colors.append([R[idx], G[idx], B[idx]])
        idx += 1
    colors = np.array(colors)
    ng_mrml_filename = os.path.join(args.inputAtlasDirectory, "clusters_location_not_given.mrml")
    wma.mrml.write(ng_list, colors, ng_mrml_filename, ratio=1.0)


if __name__ == "__main__":
    main()

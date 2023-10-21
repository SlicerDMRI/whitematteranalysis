#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import os

from joblib import Parallel, delayed

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Measure overlaps of fiber clusters with cortical parcellation or fMRI functional areas. This is based on the 3D Slicer module FiberEndPointFromLabelMap.",
        epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")
    parser.add_argument(
        'inputTractDirectory',
        help='Directory of fiber clustering results obtained by <wm_cluster_from_atlas.py> of multiple subjects. Make sure only the fiber clustering results are stored in this folder, making one subdirectory corresponding to one subject.')
    parser.add_argument(
        'inputLabelMapDirectory',
        help='Contains the parcellation or functional areas as label map files. Make sure that the input tract files and the label map files match each other in alphabetical order.')
    parser.add_argument(
        'outputDirectory',
        help='Directory of output CSV files that shows the percentage of fiber clusters\' endpoint connecting to certain regions in the label map.')
    parser.add_argument(
        'modulePath',
        help='Path of the 3D Slicer FiberEndPointFromLabelMap module.')
    parser.add_argument(
        '-j', action="store", dest="numberOfJobs", type=int,
        help='Number of processors to use.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputTractDirectory):
        print(f"Error: Input directory {args.inputTractDirectory} does not exist.")
        exit()

    if not os.path.isdir(args.inputLabelMapDirectory):
        print(f"Error: Input label map directory {args.inputLabelMapDirectory} does not exist.")
        exit()

    if not os.path.exists(args.modulePath):
        print(f"Error: FiberEndPointFromLabelMap {args.modulePath} does not exist.")
        exit()

    if args.numberOfJobs is not None:
        number_of_jobs = args.numberOfJobs
    else:
        number_of_jobs = 1

    if not os.path.exists(args.outputDirectory):
        print(f"Output directory {args.outputDirectory} does not exist, creating it.")
        os.makedirs(args.outputDirectory)

    print(f"<{os.path.basename(__file__)}>. Starting processing.")
    print("")
    print(f"=====input fiber cluster directory======\n {args.inputTractDirectory}")
    print(f"=====input label map directory======\n {args.inputLabelMapDirectory}")
    print(f"=====output directory=====\n {args.outputDirectory}")
    print(f"=====module path====\n {args.modulePath}")
    print(f"=====using N jobs: {number_of_jobs} ====\n")

    tract_dir_list = os.listdir(args.inputTractDirectory)
    tract_dir_list = sorted(tract_dir_list)

    print(f"<{os.path.basename(__file__)}> found {len(tract_dir_list)} subjects.")

    def list_label_map_files(input_dir):
        # Find input files
        input_mask = f"{input_dir}/*.nrrd"
        input_mask2 = f"{input_dir}/*.nhdr"
        input_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
        input_fnames = sorted(input_fnames)
        return(input_fnames)

    label_map_file_list = list_label_map_files(args.inputLabelMapDirectory)

    print(f"<{os.path.basename(__file__)}> found {len(label_map_file_list)} label maps.\n")

    if len(tract_dir_list) != len(label_map_file_list):
        print(f"Error: The number of subjects {len(tract_dir_list)} should be equal to the number of label maps {len(label_map_file_list)}")
        exit()

    def extract_endpoint(tract_dir, lalel_map_file, args):

        pds = wma.io.list_vtk_files(os.path.join(args.inputTractDirectory, tract_dir))
        print(f"<{os.path.basename(__file__)}> Computing: {os.path.join(args.inputTractDirectory, tract_dir)}")
        print(f"                            using {lalel_map_file}")
        print(f"                            with: {len(pds)} vtk/vtp files.")

        sub_name = os.path.split(tract_dir)[1]
        os.system(f"{args.modulePath} {lalel_map_file} {os.path.join(args.inputTractDirectory, tract_dir)} {os.path.join(args.outputDirectory, sub_name + '_endpoint.txt')} > {os.path.join(args.outputDirectory, 'log' + sub_name)}")

    Parallel(n_jobs=number_of_jobs, verbose=1)(
        delayed(extract_endpoint)(tract_dir, label_map_file, args)
        for tract_dir, label_map_file in zip(tract_dir_list, label_map_file_list))

    def list_txt_files(input_dir):
        # Find input files
        input_mask = f"{input_dir}/*.txt"
        input_fnames = glob.glob(input_mask)
        input_fnames = sorted(input_fnames)

        return input_fnames

    endpoint_txt_list = list_txt_files(args.outputDirectory)
    print(f"<{os.path.basename(__file__)}> Endpoint analysis were measured for {len(endpoint_txt_list)} subjects.")

    if len(tract_dir_list) != len(endpoint_txt_list):
        print("Error: The numbers of inputs and outputs are different. Check the log file of each subject.")
    else:
        os.system(f"rm -rf {os.path.join(args.outputDirectory, 'log*')}")


if __name__ == "__main__":
    main()

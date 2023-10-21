#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import os

import numpy as np

import whitematteranalysis as wma


def list_nhdr_files(input_dir):
    # Find input files
    input_mask = f"{input_dir}/*.nhdr"
    input_mask2 = f"{input_dir}/*.nrrd"
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Measures mean FA, etc. in tractography clusters in a directory.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
    parser.add_argument(
        'inputDirectoryDWI',
        help='A directory of DWIs (nhdr or nrrd).')
    parser.add_argument(
        'inputDirectoryMask',
        help='A directory of masks (nhdr or nrrd).')
    parser.add_argument(
        'outputDirectory',
        help='Output data will be saved here.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputDirectoryDWI):
        print(f"Error: Input directory {args.inputDirectory} does not exist.")
        exit()

    if not os.path.isdir(args.inputDirectoryMask):
        print(f"Error: Input directory {args.inputDirectory} does not exist.")
        exit()

    outdir = args.outputDirectory
    if not os.path.exists(outdir):
        print(f"<{os.path.basename(__file__)}> Output directory {outdir} does not exist, creating it.")
        os.makedirs(outdir)

    # get inputs
    dwi_list = list_nhdr_files(args.inputDirectoryDWI)
    mask_list = list_nhdr_files(args.inputDirectoryMask)

    # create FA images
    for (dwi, mask) in zip(dwi_list, mask_list):
        subject_id = os.path.splitext(os.path.basename(dwi))[0]
        fname_out_dti = os.path.join(args.outputDirectory, f'{subject_id}_DTI.nhdr')
        fname_out_b0 = os.path.join(args.outputDirectory, f'{subject_id}_B0.nhdr')
        fname_out_fa = os.path.join(args.outputDirectory, f'{subject_id}_FA.nhdr')
        print(f"/Applications/Slicer.app/Contents/lib/Slicer-4.4/cli-modules/DWIToDTIEstimation -m {mask} {dwi} {fname_out_dti} {fname_out_b0}")
        print(f"/Applications/Slicer.app/Contents/lib/Slicer-4.4/cli-modules/DiffusionTensorScalarMeasurements {fname_out_dti} {fname_out_fa} -e FractionalAnisotropy &")


if __name__ == "__main__":
    main()

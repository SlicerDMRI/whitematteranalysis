#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import os

import pandas as pd
import vtk

import whitematteranalysis as wma
from whitematteranalysis.anatomy.org_atlas_utils import (
    ORGAtlasBundleFileHeading, add_org_atlas_prefix, get_bundle_short_name,
    get_commissural_augmented_bundles, get_hemispheric_mono_bundles)
from whitematteranalysis.data.atlas.utils import (ORGAtlasVersion,
                                                  get_local_atlas_bundle_fname)


def get_hemispheric_mono_bundle_names(_df):

    hemis_mono_df = get_hemispheric_mono_bundles(_df)
    hemis_names = hemis_mono_df[ORGAtlasBundleFileHeading.SHORT_NAME.value].values.tolist()
    return add_org_atlas_prefix(hemis_names)


def get_commissural_bundle_names(_df):

    com_df = get_commissural_augmented_bundles(_df)
    com_names = get_bundle_short_name(com_df, com_df[
        ORGAtlasBundleFileHeading.ID.value].values)
    return add_org_atlas_prefix(com_names)


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Append multiple fiber clusters into anatomical tracts based on the ORG atlas.",
        epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")
    parser.add_argument(
        'inputDirectory',
        help='Contains fiber clusters as vtkPolyData file(s).')
    parser.add_argument(
        'atlasDirectory',
        help='The ORG atlas folder that contains the anatomical tract MRML files.')
    parser.add_argument(
        'outputDirectory',
        help='The output directory should be a new empty directory. It will be created if needed.')
    parser.add_argument(
        '--version',
        choices=ORGAtlasVersion._member_names_,
        default=ORGAtlasVersion.V1_2.name,
        help="Atlas version.",
    )

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    inputdir = os.path.abspath(args.inputDirectory)
    if not os.path.isdir(args.inputDirectory):
        print(f"<{os.path.basename(__file__)}> Error: Input directory {args.inputDirectory} does not exist.")
        exit()

    inputdir_left = os.path.join(inputdir, 'tracts_left_hemisphere')
    inputdir_right = os.path.join(inputdir, 'tracts_right_hemisphere')
    inputdir_comm = os.path.join(inputdir, 'tracts_commissural')

    if not os.path.isdir(inputdir_left):
        print(f"<{os.path.basename(__file__)}> Error: Input directory {inputdir_left} does not exist.")
        exit()
    if not os.path.isdir(inputdir_right):
        print(f"<{os.path.basename(__file__)}> Error: Input directory {inputdir_right}does not exist.")
        exit()
    if not os.path.isdir(inputdir_comm):
        print(f"<{os.path.basename(__file__)}> Error: Input directory {inputdir_comm} does not exist.")
        exit()

    atlasdir = os.path.abspath(args.atlasDirectory)
    if not os.path.isdir(args.atlasDirectory):
        print(f"<{os.path.basename(__file__)}> Error: Atlas directory {args.atlasDirectory} does not exist.")
        exit()

    def list_mrml_files(input_dir):
        # Find input files
        input_mask = f"{input_dir}/T*.mrml"
        input_mrml_fnames = glob.glob(input_mask)
        return (input_mrml_fnames)

    mrml_files = list_mrml_files(atlasdir)

    if len(mrml_files) == 0:
        print(f"<{os.path.basename(__file__)}> Error: There is no mrml files in the input atlas folder.")
    else:
        print(f"<{os.path.basename(__file__)}> {len(mrml_files)-1} mrml files are detected.")

    outdir = os.path.abspath(args.outputDirectory)
    if not os.path.exists(args.outputDirectory):
        print(f"<{os.path.basename(__file__)}> Output directory {args.outputDirectory} does not exist, creating it.")
        os.makedirs(outdir)

    def output_appended_tract(cluster_vtp_list, outputfile):
        appender = vtk.vtkAppendPolyData()
        for c_idx in range(len(cluster_vtp_list)):
            cluster_vtp = cluster_vtp_list[c_idx]
            pd_cluster = wma.io.read_polydata(cluster_vtp)

            vtk_array = vtk.vtkIntArray()
            vtk_array.SetName('cluster_idx')
            for p_idx in range(0, pd_cluster.GetNumberOfPoints()):
                vtk_array.InsertNextTuple1(int(c_idx))

            pd_cluster.GetPointData().AddArray(vtk_array)

            #print '<wm_append_clusters_to_anatomical_tracts>', cluster_vtp, ', number of fibers', pd_cluster.GetNumberOfLines()

            if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                appender.AddInputData(pd_cluster)
            else:
                appender.AddInput(pd_cluster)

        appender.Update()
        pd_appended_cluster = appender.GetOutput()

        wma.io.write_polydata(pd_appended_cluster, outputfile)

    fname = get_local_atlas_bundle_fname(ORGAtlasVersion(args.version))
    df = pd.read_csv(fname, sep=",")

    hemispheric_tracts = get_hemispheric_mono_bundle_names(df)
    commissural_tracts = get_commissural_bundle_names(df)

    print(f"<{os.path.basename(__file__)}> hemispheric tracts (left and right): ")
    tract_idx = 1
    for tract in hemispheric_tracts:

        print(f" * {tract_idx} - {tract}")
        tract_idx = tract_idx + 1
        mrml = os.path.join(atlasdir, tract+".mrml")

        if not os.path.exists(mrml):
            print(f"<{os.path.basename(__file__)}> Error: Cannot locate {mrml}")
            exit()

        cluster_vtp_list_left = list()
        cluster_vtp_list_right = list()
        f = open(mrml)
        for line in f:
            idx = line.find('.vtp')
            if idx > 0:
                cluster_vtp_filename = line[idx-13:idx+4]

                cluster_vtp_filename_left = os.path.join(inputdir_left, cluster_vtp_filename)
                if not os.path.exists(cluster_vtp_filename_left):
                    print(f"<{os.path.basename(__file__)}> Error: {cluster_vtp_filename_left} does not exist.")
                    exit()
                cluster_vtp_list_left.append(cluster_vtp_filename_left)

                cluster_vtp_filename_right = os.path.join(inputdir_right, cluster_vtp_filename)
                if not os.path.exists(cluster_vtp_filename_right):
                    print(f"<{os.path.basename(__file__)}> Error: {cluster_vtp_filename_right} does not exist.")
                    exit()
                cluster_vtp_list_right.append(cluster_vtp_filename_right)

        output_tract_left = os.path.join(outdir, tract + '_left.vtp')
        output_tract_right = os.path.join(outdir, tract + '_right.vtp')

        output_appended_tract(cluster_vtp_list_left, output_tract_left)
        output_appended_tract(cluster_vtp_list_right, output_tract_right)

    print(f"<{os.path.basename(__file__)}> commissural tracts: ")
    for tract in commissural_tracts:

        print(f" * {tract_idx} - {tract}")
        tract_idx = tract_idx + 1

        mrml = os.path.join(atlasdir, tract+".mrml")

        if not os.path.exists(mrml):
            print(f"<{os.path.basename(__file__)}> Error: Cannot locate {mrml}")
            exit()

        cluster_vtp_list_comm = list()
        f = open(mrml)
        for line in f:
            idx = line.find('.vtp')
            if idx > 0:
                cluster_vtp_filename = line[idx-13:idx+4]

                cluster_vtp_filename_comm = os.path.join(inputdir_comm, cluster_vtp_filename)
                if not os.path.exists(cluster_vtp_filename_comm):
                    print(f"<{os.path.basename(__file__)}> Error {cluster_vtp_filename_comm} does not exist.")
                    exit()
                cluster_vtp_list_comm.append(cluster_vtp_filename_comm)

        output_tract_comm = os.path.join(outdir, tract + '.vtp')

        output_appended_tract(cluster_vtp_list_comm, output_tract_comm)

    def list_cluster_files(input_dir):
        # Find input files
        input_mask = f"{input_dir}/T_*.vtk"
        input_mask2 = f"{input_dir}/T_*.vtp"
        input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
        input_pd_fnames = sorted(input_pd_fnames)
        return(input_pd_fnames)

    list_tracts= list_cluster_files(outdir)

    print('')
    print(f'<{os.path.basename(__file__)}> Appended tracts can be found at {outdir}\n')
    print(f'<{os.path.basename(__file__)}> A total of {len(list_tracts)} tracts\n')


if __name__ == "__main__":
    main()

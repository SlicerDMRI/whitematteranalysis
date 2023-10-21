#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import os
import warnings

import numpy as np

import whitematteranalysis as wma
from whitematteranalysis.utils.opt_pckg import optional_package

matplotlib, have_mpl, _ = optional_package("matplotlib")
plt, _, _ = optional_package("matplotlib.pyplot")

if have_mpl:
    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use("Agg")
else:
    warnings.warn(matplotlib._msg)
    warnings.warn("Cannot plot quality control data.")


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Quality control of fiber clustering results across multiple subjects.",
        epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")
    parser.add_argument(
        'inputDirectory',
        help='Directory of fiber clustering results obtained by <wm_cluster_from_atlas.py> of multiple subjects. Make sure only the fiber clustering results are stored in this folder, making one subdirectory corresponding to one subject.')
    parser.add_argument(
        'outputDirectory',
        help='Quality control information will be stored in the output directory, which will be created if it does not exist.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputDirectory):
        print(f"Error: Input directory {args.inputDirectory} does not exist.")
        exit()
    
    output_dir = args.outputDirectory
    if not os.path.exists(output_dir):
        print(f"<{os.path.basename(__file__)}> Output directory {output_dir} does not exist, creating it.")
        os.makedirs(output_dir)
    
    subject_list = os.listdir(args.inputDirectory)
    print(f"<{os.path.basename(__file__)}> found {len(subject_list)} subjects.")
    
    # Check if all subjects have the same number of clusters.
    # This can also help find the sub folders that are not the fiber clustering results.
    num_of_subjects = len(subject_list)
    flag_same_num_clusters = 1
    for sidx in range(0, num_of_subjects):
        sub = subject_list[sidx]
        sub_dir = os.path.join(args.inputDirectory, sub)
        input_mask = f"{sub_dir}/cluster_*.vtk"
        input_mask2 = f"{sub_dir}/cluster_*.vtp"
        cluster_polydatas = glob.glob(input_mask) + glob.glob(input_mask2)
        cluster_polydatas = sorted(cluster_polydatas)
        print(f"   {sub} has {len(cluster_polydatas)} clusters.")
        if sidx == 0:
            num_of_clusters = len(cluster_polydatas)
        else:
            if num_of_clusters != len(cluster_polydatas):
                flag_same_num_clusters = 0
    
    if flag_same_num_clusters == 0:
        print("Error: All subjects should have the number of clusters. (Hint: there may be folders that are not the fiber clustering results.)")
        exit()
    
    # Read number of fibers per cluster per subject
    print(f"<{os.path.basename(__file__)}> calculate the number of fibers per cluster per subject.")
    num_fibers_per_subject = np.zeros([num_of_subjects, num_of_clusters])
    for sidx in range(0, num_of_subjects):
        sub = subject_list[sidx]
        print(f"   loading {sub}")
        sub_dir = os.path.join(args.inputDirectory, sub)
        input_mask = f"{sub_dir}/cluster_*.vtk"
        input_mask2 = f"{sub_dir}/cluster_*.vtp"
        cluster_polydatas = glob.glob(input_mask) + glob.glob(input_mask2)
        cluster_polydatas = sorted(cluster_polydatas)
        for cidx in range(0, num_of_clusters):
            fname = cluster_polydatas[cidx]
            pd = wma.io.read_polydata(fname)
            num_fibers_per_subject[sidx, cidx] = pd.GetNumberOfLines()
    
    subjects_per_cluster = np.sum(num_fibers_per_subject > 0, axis=0)
    #clusters_per_subject = np.sum(num_fibers_per_subject > 0, axis=1)
    
    percent_subjects_per_cluster = np.divide(subjects_per_cluster, float(num_of_subjects))
    
    clusters_qc_fname = os.path.join(output_dir, 'cluster_quality_control.txt')
    print(f"<{os.path.basename(__file__)}> Saving cluster quality control information file.")
    clusters_qc_file = open(clusters_qc_fname, 'w')
    print('cluster_idx','\t', 'number_subjects', 'percent_subjects', file=clusters_qc_file)
    for cidx in range(0, num_of_clusters):
        print(cidx + 1,'\t', subjects_per_cluster[cidx],'\t', percent_subjects_per_cluster[cidx] * 100.0, file=clusters_qc_file)
    clusters_qc_file.close()
    
    if have_mpl:
        print(f"<{os.path.basename(__file__)}> Saving subjects per cluster histogram.")
        fig, ax = plt.subplots()
        counts = np.zeros(num_of_subjects+1)
        counts[:np.max(subjects_per_cluster)+1] = np.bincount(subjects_per_cluster)
        ax.bar(list(range(num_of_subjects + 1)), counts, width=1, align='center')
        ax.set(xlim=[-1, num_of_subjects + 1])
        plt.title('Histogram of Subjects per Cluster')
        plt.xlabel('subjects per cluster')
        plt.ylabel('number of clusters')
        plt.savefig(os.path.join(output_dir, 'subjects_per_cluster_hist.pdf'))
        plt.close()

if __name__ == '__main__':
    main()

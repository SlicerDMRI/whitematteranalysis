#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import numpy as np
import scipy.stats

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Perform quality control steps (check if all subjects have matching data) for .csv or txt measurement files in the input directory.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
    parser.add_argument(
        'inputDirectory',
        help='A directory of .csv or txt measurement files from Slicer tractography measurements.')
    #parser.add_argument(
    #    'outputDirectory',
    #    help='Quality control information will be stored in the output directory, which will be created if it does not exist.')
    
    parser.add_argument(
        '-outlier_std', action="store", dest="OutlierStandardDeviation", type=float, default=2.0,
        help='Standard deviation for outlier detection of subjects.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputDirectory):
        print(f"Error: Input directory {args.inputDirectory} does not exist.")
        exit()
    
    print(f"\n\n== Testing measurement files in directory: {args.inputDirectory}")
    
    ## output_dir = args.outputDirectory
    ## if not os.path.exists(output_dir):
    ##     print(f"<{os.path.basename(__file__)}> Output directory {output_dir} does not exist, creating it.")
    ##     os.makedirs(output_dir)
    
    measurement_list = wma.tract_measurement.load_measurement_in_folder(args.inputDirectory, hierarchy = 'Column', separator = 'Tab')
    
    print(f"\n== Number of subjects data found: {len(measurement_list)}")
    
    if len(measurement_list) < 1:
        print(f"ERROR, no measurement files found in directory {args.inputDirectory}")
        exit()
    
    # print the header
    print("\n==Measurement header found. Measurements are:")
    header = measurement_list[0].measurement_header
    for variable in header:
        print(variable)
    
    # sanity check all measured variables are the same
    print("\n== Testing if all subjects have the same measurement header.")
    test = True
    for subject_measured in measurement_list:
        all_same = all(subject_measured.measurement_header == header)
        if not all_same:
            print(f"ERROR: Subject does not have same measurement header: {subject_measured.case_id}")
            test = False
    if test:
        print("Passed. All subjects have the same measurement header.")
    
    # Make subject id list for reporting of any potential outliers
    subject_id_list = []
    for subject_measured in measurement_list:
        #fname = subject_measured.cluster_path[0]
        #subject_id_list.append(os.path.basename(os.path.split(fname)[0]))
        subject_id_list.append(subject_measured.case_id)
    subject_id_list = np.array(subject_id_list)
    
    # sanity check number of clusters is the same for all subjects
    print("\n== Testing if all subjects have the same number of clusters.")
    test = True
    number_of_clusters = measurement_list[0].measurement_matrix[:,0].shape[0]
    print(f"Group number of clusters (from first subject): {number_of_clusters}")
    for subject_measured, id in zip(measurement_list, subject_id_list):
        nc = subject_measured.measurement_matrix[:,0].shape[0]
        if not nc == number_of_clusters:
            print(nc, " : ", id)
            test = False
    if test:
        print(f"Passed. All subjects have the same number of clusters ({number_of_clusters}).")
    else:
        print("ERROR: All subjects do not have the same number of clusters. There was an earlier error in clustering or measurement that must be fixed before analysis.")
    
    # sanity check numbers of fibers are ok for all subjects
    print("\n== Testing if all subjects have reasonable mean fibers per cluster.")
    print(f"Will print any subjects with more than {args.OutlierStandardDeviation} standard deviations below group mean.")
    test = True
    vidx = list(header).index('Num_Fibers')
    mean_fibers_list = []
    for subject_measured in measurement_list:
        # mean number of fibers per cluster
        mean_fibers_list.append(np.mean(subject_measured.measurement_matrix[:,vidx]))
    # print lowest numbers of fibers
    mean_fibers_list = np.array(mean_fibers_list)
    mean_fibers = np.mean(mean_fibers_list)
    std_fibers = np.std(mean_fibers_list)
    print(f"Mean and standard deviation of mean fibers per cluster in group: {mean_fibers} +/- {std_fibers}")
    threshold = mean_fibers - args.OutlierStandardDeviation * std_fibers
    ## mean_sorted_idx = np.argsort(np.array(mean_fibers_list))
    ## for idx in mean_sorted_idx:
    ##     print mean_fibers_list[idx], "  :  ", subject_id_list[idx]
    sorted_idx = np.argsort(mean_fibers_list)
    for mf, sp in zip(mean_fibers_list[sorted_idx], subject_id_list[sorted_idx]):
        if mf < threshold:
            if test:
                print(f"Subject(s) found with mean fibers more than {args.OutlierStandardDeviation} standard deviations below group mean.")
            print(f"{mf}   :   {sp}")
            test = False
    if test:
        print("Passed. All subjects have reasonable mean fibers per cluster.")
            
    # sanity check number of empty clusters is not too large
    print("\n== Testing if all subjects have reasonable numbers of empty clusters.")
    print(f"Will print any subjects with more than {args.OutlierStandardDeviation} standard deviations above group mean.")
    test = True
    empty_clusters_list = []
    for subject_measured in measurement_list:
        # check if number of fibers per cluster is 0 to indicate empty cluster
        empty_clusters_list.append(np.sum(subject_measured.measurement_matrix[:,vidx] == 0))
    # print highest numbers of empty clusters
    ## empty_sorted_idx = np.argsort(np.array(empty_clusters_list))
    ## for idx in empty_sorted_idx[::-1]:
    ##     print empty_clusters_list[idx], "  :  ", subject_id_list[idx]
    empty_clusters_list = np.array(empty_clusters_list)
    mean_clusters = np.mean(empty_clusters_list)
    std_clusters = np.std(empty_clusters_list)
    print(f"Mean and standard deviation of empty cluster number in group: {mean_clusters} +/- {std_clusters}")
    threshold = mean_clusters + args.OutlierStandardDeviation * std_clusters
    sorted_idx = np.argsort(empty_clusters_list)
    for mf, sp in zip(empty_clusters_list[sorted_idx], subject_id_list[sorted_idx]):
        if mf > threshold:
            if test:
                print(f"Subject(s) found with empty clusters more than {args.OutlierStandardDeviation} standard deviations above group mean.")
            print(f"{mf}   :   {sp}")
            test = False
    if test:
        print("Passed. All subjects have reasonable numbers of empty clusters.")
    
    print("\nPlease use reasonable judgement and visual inspection of data to decide if any subject datasets are outliers. Remember, 5% of a normally-distributed dataset will fall 2 standard deviations away from the mean, so most datasets will have some ""outliers"". Also try the parameter -outlier_std 3 to see if any datasets are 3 std from the mean.")

if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import csv
import os

import numpy as np

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Export all subjects selected measures into one excel-readable file for use in statistics programs. Please verify your input by running wm_quality_control_cluster measurements before running this program.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
    parser.add_argument(
        'inputDirectory',
        help='A directory of .csv or txt measurement files from Slicer tractography measurements.')
    parser.add_argument(
        '-m', action="store", nargs='+', dest="measures", default='tensor1.FractionalAnisotropy',
        help='Chosen measure(s) for export. Examples include tensor1.FractionalAnisotropy, tensor2.FractionalAnisotropy, FreeWater, Num_Fibers, tensor1.Trace, tensor1.MaxEigenvalue, etc. Run wm_quality_control_cluster_measurement.py to print available measures in the dataset. ')
    parser.add_argument(
        'subjectsInformationFile',
        help='A file in excel format: column 1 must be subject ID.')
    parser.add_argument(
        'outputFileName',
        help='A file name for output (ending in .txt or .csv).')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputDirectory):
        print(f"<{os.path.basename(__file__)}> Error: Input directory {args.inputDirectory} does not exist.")
        exit()

    if not os.path.exists(args.subjectsInformationFile):
        print(f"<{os.path.basename(__file__)}> Error: Input file {args.subjectsInformationFile} does not exist.")
        exit()

    measurements = args.measures

    # Read and check data
    measurement_list = wma.tract_measurement.load_measurement_in_folder(args.inputDirectory, hierarchy = 'Column', separator = 'Tab')
    print(f"Measurement directory {args.inputDirectory}")
    number_of_subjects = len(measurement_list)
    print(f"Number of subjects data found: {number_of_subjects}")
    if number_of_subjects < 1:
        print(f"ERROR, no measurement files found in directory {args.inputDirectory}")
        exit()
    header = measurement_list[0].measurement_header
    region_list = measurement_list[0].cluster_path
    #print "Measurement header is:", header
    #print "Clusters found:",  measurement_list[0].cluster_path
    number_of_clusters = measurement_list[0].measurement_matrix[:,0].shape[0]
    print(f"Number of measurement regions (clusters and hierarchy groups) is: {number_of_clusters}")

    # Read and check subject ID list
    dg = wma.tract_measurement.load_demographics(args.subjectsInformationFile)
    case_id_list = dg.case_id_list

    # Check subject IDs from excel versus the input directory.
    print("Checking subject IDs match in excel file and measurement directory.")
    if len(case_id_list) ==  number_of_subjects:
            print("Subject counts in excel subject ID list and measurement directory match.")
    else:
        print("ERROR: Subject counts in excel subject ID list and measurement directory don't match:", end=' ')
        print(len(case_id_list), number_of_subjects)
        exit()
    for subject_measured, subject_id in zip(measurement_list, case_id_list):
        if not str(subject_id) in subject_measured.case_id:
            print("ERROR: id list and input data mismatch.")
            print(f"ERROR at: {subject_measured.case_id} {subject_id}")
            exit()
    print("Dataset passed. Subject IDs in subject information excel file match subject IDs in measurement directory.")

    # get data for export
    vidx_list = []; # index of values of interest
    for measure in measurements:
        vidx_list.append(list(header).index(measure))
    print(f"Column indices of measures for export: {vidx_list}")

    # reformat this information for export.
    # export format for header is subject ID, region1.measure1, region1.measure2, ..., regionN.measureN
    # export format for row is subject ID, measures...

    output_header = []
    output_header.append('SubjectID')
    for region in region_list:
        print(region)
        for vidx in vidx_list:
            measure_name = header[vidx]
            output_header.append(region+'.'+measure_name)

    print(output_header)
    with open(args.outputFileName, 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter='\t')
        spamwriter.writerow(output_header)
        for subject_id, subject_measured in zip(case_id_list, measurement_list):
            output_data = []
            output_data.append(subject_id)
            region_id = 0
            for region in region_list:
                for vidx in vidx_list:
                    output_data.append(subject_measured.measurement_matrix[region_id,vidx])
                region_id += 1
            spamwriter.writerow(output_data)


if __name__ == "__main__":
    main()

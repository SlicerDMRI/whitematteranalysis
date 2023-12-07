#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import os

import numpy as np
import pandas


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Append diffusion measure files from multiple subjects",
        epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu.")
    parser.add_argument(
        'inputfolder',
        help='A folder contains WMA output by subjects. Each subject should be one folder.')
    parser.add_argument(
        'appenedmeasurefile',
        help='Output csv file.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    subject_folders = sorted(glob.glob(os.path.join(args.inputfolder, '*')))

    subject_folders = [d for d in subject_folders if os.path.isdir(d)]

    subject_IDs = [os.path.basename(folder) for folder in subject_folders]
    print()
    print(f'Subject IDs: (n={len(subject_IDs)}): ')
    print(subject_IDs)

    # anatomical tracts
    stats = pandas.read_table(os.path.join(subject_folders[0], 'AnatomicalTracts/diffusion_measurements_anatomical_tracts.csv'), delimiter=' , ')

    fields = [col.strip() for col in stats.columns]
    fields = fields[1:]
    print()
    print(f'Tract diffusion measures per subject: (n={len(fields)}): ')
    print(fields)

    tracts = stats.to_numpy()[:, 0]
    tracts = [os.path.basename(filepath).replace('.vtp', '').replace('.vtk', '').strip() for filepath in tracts]
    print()
    print(f'White matter tracts per subject: (n={len(tracts)}): ')
    print(tracts)

    appended_fields = ['subjectkey']
    for tract in tracts:
        for field in fields:
            appended_fields.append(f"{tract}.{field}")

    print()
    print(f'Appended tract diffusion measure fields per subject: (n={len(appended_fields)}): ')
    print(appended_fields[:10], ' ... ')

    data = []
    for s_idx, subject_folder in enumerate(subject_folders):
        print(' * loading tract feature for subject #%04d (subject ID - %s):' %(s_idx, subject_IDs[s_idx]))
        csv = os.path.join(subject_folder, 'AnatomicalTracts/diffusion_measurements_anatomical_tracts.csv')
        stats = pandas.read_table(csv, delimiter=' , ')

        stats_data = stats.to_numpy()[:, 1:]
        stats_data_vec = stats_data.flatten()

        if stats_data_vec.shape[0] != len(appended_fields) - 1:
            print("Error: Check if the diffusion measure file has the same rows and columns with other subjects!")
            exit()

        data.append(stats_data_vec)

    data = np.array(data)
    data = np.concatenate([np.array(subject_IDs)[:, np.newaxis], data], axis = 1)

    df = pandas.DataFrame(data, columns=appended_fields)

    print()
    print('Appended tract diffusion measures:')
    print(df)

    df.to_csv(os.path.abspath(args.appenedmeasurefile).replace('.csv', '_anatomical_tracts.csv'), index=False)
    print()
    print('Output file at:', os.path.abspath(args.appenedmeasurefile).replace('.csv', '_anatomical_tracts.csv'))

    # Fiber clusters
    stats = pandas.read_table(os.path.join(subject_folders[0], 'FiberClustering/SeparatedClusters/diffusion_measurements_commissural.csv'), delimiter=' , ')

    fields = [col.strip() for col in stats.columns]
    fields = fields[1:]
    print()
    print(f'Fiber cluster diffusion measures per subject: (n={len(fields)}): ')
    print(fields)

    clusters = stats.to_numpy()[:, 0]
    clusters = [os.path.basename(filepath).replace('.vtp', '').replace('.vtk', '').strip() for filepath in clusters]
    print()
    print(f'Fiber clusters per subject: (n={len(clusters)}): ')

    if len(clusters) == 800 and len(tracts) == 73:
        comm_clusters = [3, 8, 33, 40, 46, 52, 56, 57, 62, 68, 69, 73, 86, 91, 109, 110, 114, 142, 144, 145, 146, 159, 163, 250, 251, 252, 257, 262, 263, 268, 271, 305, 311, 314, 322, 330, 334, 338, 342, 350, 363, 371, 375, 403, 410, 437, 440, 448, 456, 465, 468, 475, 484, 485, 488, 519, 522, 525, 543, 545, 549, 557, 576, 577, 582, 587, 591, 597, 601, 614, 620, 623, 627, 633, 645, 653, 658, 663, 670, 677, 683, 702, 770, 781]
        comm_clusters = np.array(comm_clusters)
        hemi_clusters = np.setdiff1d(np.arange(800), comm_clusters)
    elif len(clusters) == 800 and len(tracts) == 74: # CPC separated
        comm_clusters = [3, 8, 33, 40, 46, 52, 56, 57, 62, 68, 69, 73, 86, 91, 109, 110, 114, 142, 144, 146, 163, 250, 251, 252, 257, 262, 263, 268, 271, 305, 311, 314, 322, 330, 334, 338, 342, 350, 363, 371, 375, 403, 410, 437, 440, 448, 456, 465, 468, 475, 484, 485, 488, 519, 522, 525, 543, 545, 549, 576, 577, 582, 587, 591, 597, 601, 614, 620, 623, 627, 633, 645, 653, 658, 663, 670, 683, 702, 781]
        comm_clusters = np.array(comm_clusters)
        hemi_clusters = np.setdiff1d(np.arange(800), comm_clusters)
    else:
        comm_clusters = None
        hemi_clusters = None

    locations = ['left_hemisphere', 'right_hemisphere', 'commissural']
    appended_fields = ['subjectkey']
    clusters = np.array(clusters)
    for loc in locations:
        if loc == 'left_hemisphere' or loc == 'right_hemisphere':
            clusters_loc = clusters[hemi_clusters]
        else:
            clusters_loc = clusters[comm_clusters]

        for cluster in clusters_loc:
            for field in fields:
                appended_fields.append(f"{loc}.{cluster}.{field}")

    print()
    print(f'Appended cluster diffusion measure fields per subject: (n={len(appended_fields)}): ')

    data = []
    for s_idx, subject_folder in enumerate(subject_folders):
        print(' * loading cluster feature for subject #%04d (subject ID - %s):' %(s_idx, subject_IDs[s_idx]))
        csv = os.path.join(subject_folder, 'FiberClustering/SeparatedClusters/diffusion_measurements_left_hemisphere.csv')
        stats = pandas.read_table(csv, delimiter=' , ')

        stats_data = stats.to_numpy()[:, 1:]
        stats_data = stats_data[hemi_clusters, :]
        l_stats_data_vec = stats_data.flatten()

        csv = os.path.join(subject_folder, 'FiberClustering/SeparatedClusters/diffusion_measurements_right_hemisphere.csv')
        stats = pandas.read_table(csv, delimiter=' , ')

        stats_data = stats.to_numpy()[:, 1:]
        stats_data = stats_data[hemi_clusters, :]
        r_stats_data_vec = stats_data.flatten()

        csv = os.path.join(subject_folder, 'FiberClustering/SeparatedClusters/diffusion_measurements_commissural.csv')
        stats = pandas.read_table(csv, delimiter=' , ')

        stats_data = stats.to_numpy()[:, 1:]
        stats_data = stats_data[comm_clusters, :]
        c_stats_data_vec = stats_data.flatten()

        stats_data_vec = np.concatenate([l_stats_data_vec, r_stats_data_vec, c_stats_data_vec])

        if stats_data_vec.shape[0] != len(appended_fields) - 1:
            print("Error: Check if the diffusion measure file has the same rows and columns with other subjects!")
            exit()

        data.append(stats_data_vec)

    data = np.array(data)
    data = np.concatenate([np.array(subject_IDs)[:, np.newaxis], data], axis = 1)

    df = pandas.DataFrame(data, columns=appended_fields)

    print()
    print('Appended cluster diffusion measures:')
    print(df)

    df.to_csv(os.path.abspath(args.appenedmeasurefile).replace('.csv', '_fiber_clusters.csv'), index=False)
    print()
    print('Output file at:', os.path.abspath(args.appenedmeasurefile).replace('.csv', '_fiber_clusters.csv'))


if __name__ == "__main__":
    main()

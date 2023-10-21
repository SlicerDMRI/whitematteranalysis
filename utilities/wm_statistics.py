#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import statsmodels.sandbox.stats.multicomp

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Compute requested statistics comparing groups in the input groups file. Please verify your input by running wm_quality_control_cluster measurements before running this program.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
    parser.add_argument(
        'inputDirectory',
        help='A directory of .csv or txt measurement files from Slicer tractography measurements.')
    parser.add_argument(
        'subjectsInformationFile',
        help='A file in excel format: column 1 must be subject ID, column 2 must be group label, and columns 3, 4, etc. contain other relevant information such as age or covariates.')
    parser.add_argument(
        '-atlas',  action="store", dest="atlasDirectory", default=None,
        help='A directory containing a whitematteranalysis atlas (atlas.vtp, clusters, etc.). This allows output of a MRML file if whole-brain cluster statistics were performed.')
    parser.add_argument(
        '-m', action="store", dest="measure", default='tensor1.FractionalAnisotropy',
        help='Chosen measure(s) for statistics. Examples include tensor1.FractionalAnisotropy, tensor2.FractionalAnisotropy, FreeWater, Num_Fibers, tensor1.Trace, tensor1.MaxEigenvalue, etc. Run wm_quality_control_cluster_measurement.py to print available measures in the dataset. ')
    parser.add_argument(
        '-alpha', action="store", dest="FDR_alpha", type=float, default=0.1,
        help='For FDR multiple comparison correction, the allowable percentage on average of false positives. A number from 0 to 1. Values of 0.05 or 0.10 are useful.')
    parser.add_argument('--hierarchy', dest='hierarchy', action='store_true')
    parser.add_argument('--no-hierarchy', dest='hierarchy', action='store_false')
    parser.set_defaults(hierarchy=True)
    parser.add_argument('-mode', dest='mode', action='store', default='all', help='all, clusters, toplevel, sublevel')
    parser.add_argument('--onetail', dest='one_tailed', action='store_true')
    parser.set_defaults(one_tailed=False)

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
        print(f"<{os.path.basename(__file__)}> Error: Input file", args.subjectsInformationFile, "does not exist.")
        exit()

    if args.atlasDirectory:
        if not os.path.isdir(args.atlasDirectory):
            print(f"<{os.path.basename(__file__)}> Error: Input directory {args.atlasDirectory} does not exist.")
            exit()

    measurement = args.measure
    fdr_q = args.FDR_alpha
    analyze_hierarchy = args.hierarchy
    mode = args.mode
    atlas_directory = args.atlasDirectory

    # Read and check data
    measurement_list = wma.tract_measurement.load_measurement_in_folder(args.inputDirectory, hierarchy = 'Column', separator = 'Tab')
    print(f"Measurement directory {args.inputDirectory}")
    number_of_subjects = len(measurement_list)
    print(f"Number of subjects data found: {number_of_subjects}")
    if number_of_subjects < 1:
        print(f"ERROR, no measurement files found in directory {args.inputDirectory}")
        exit()
    header = measurement_list[0].measurement_header
    #print "Measurement header is:", header
    #print "Clusters found:",  measurement_list[0].cluster_path
    number_of_clusters = measurement_list[0].measurement_matrix[:,0].shape[0]
    print(f"Number of measurement regions (clusters and hierarchy groups) is: {number_of_clusters}")

    # Read and check subject ID and other information
    dg = wma.tract_measurement.load_demographics(args.subjectsInformationFile)
    case_id_list = dg.case_id_list
    group_id_list = dg.group_id_list
    age_list = dg.get_demographics_by_header('Age')
    age_list = list(map(float, age_list))

    # Check subject IDs from excel versus the input directory.
    print("Checking subject information in input file and directory.")
    if len(case_id_list) == len(group_id_list) & len(group_id_list) == number_of_subjects:
            print("Subject counts in excel subject ID list, group list, and measurement directory match.")
    else:
        print("ERROR: Subject counts in excel subject ID list, group list, and measurement directory don't match:", end=' ')
        print(len(case_id_list), len(group_id_list), number_of_subjects)
        exit()
    for subject_measured, subj_id, group in zip(measurement_list, case_id_list, group_id_list):
        #print subject_measured.case_id, subj_id, group
        if not str(subj_id) in subject_measured.case_id:
            print("ERROR: id list and input data mismatch.")
            print(f"ERROR at: {subject_measured.case_id} {subj_id} {group}")
            exit()
    print("Dataset passed. Subject IDs in subject information file match subject IDs in measurement directory.")

    # Figure out what the groups are
    study_groups = list(set(group_id_list))
    print("Study groups found:", study_groups)
    if len(study_groups) > 2:
        print("ERROR: this code can currently only handle two study groups.")
        exit()
    print(f"Subjects per group: {np.sum(np.array(group_id_list) == study_groups[0])} {np.sum(np.array(group_id_list) == study_groups[1])}")

    print(f"Performing statistics for measure: {measurement}")

    # Note this hierarchy code will be removed/simplified when the Slicer measurement file is simplified
    # ------
    # Now analyze only the hierarchy or only the clusters. Find the appropriate indices.
    measurement_list[0].measurement_matrix[:,0].shape[0]
    hierarchy_cluster_list = []
    hierarchy_region_list = []
    hierarchy_toplevel_list = []
    hierarchy_cluster_name_list = []
    hierarchy_region_name_list = []
    hierarchy_toplevel_name_list = []
    label_list = []

    hierarchy_depth = []
    idx = 0
    for region in measurement_list[0].cluster_path:
        # if it has the word :cluster_, it's a cluster leaf of a hierarchy...
        if ":cluster_" in region:
            hierarchy_cluster_list.append(idx)
            hierarchy_cluster_name_list.append(region)
        elif ":" in region:
            if analyze_hierarchy:
                # this is a sublevel of the hierarchy, temporarily assume it's what we want, and all are on same level
                hierarchy_region_list.append(idx)
                hierarchy_region_name_list.append(region)
        else:
            if analyze_hierarchy:
                # grab only the top level for the stats
                # this will be either a hierarchy node or a regular cluster with no hierarchy
                hierarchy_toplevel_list.append(idx)
                hierarchy_toplevel_name_list.append(region)
        # figure out how many levels in this hierarchy. Stats should be performed at each level only
        hierarchy_depth.append(region.count(":"))
        idx += 1

    #print hierarchy_region_name_list
    #print hierarchy_region_list

    print("Toplevel:")
    #print hierarchy_toplevel_name_list

    # TODO: use this only if hierarchy
    #mode = "sublevel"
    #mode = "toplevel"
    #mode = "clusters"
    #mode = "all"
    # sublevels
    if mode == "sublevel":
        regions_for_stats = np.array(hierarchy_region_list)
        names_for_stats = hierarchy_region_name_list
    elif mode == "toplevel":
        # top level
        regions_for_stats = np.array(hierarchy_toplevel_list)
        names_for_stats = hierarchy_toplevel_name_list
    elif mode == "clusters":
        # clusters
        regions_for_stats = np.array(hierarchy_cluster_list)
        names_for_stats = hierarchy_cluster_name_list
    elif mode == "all":
        # we are measuring using cluster_*.vtp not a hierarchy
        regions_for_stats = np.array(list(range(number_of_clusters)))
        names_for_stats = []
        for fname in measurement_list[0].cluster_path:
            #print os.path.splitext(fname)
            #print os.path.split(os.path.split(os.path.splitext(fname)[0])[1])[1]
            names_for_stats.append(os.path.split(os.path.split(os.path.splitext(fname)[0])[1])[1])
        #names_for_stats = measurement_list[0].cluster_path
    number_of_tests = len(regions_for_stats)
    # Note this hierarchy code will be removed/simplified when the Slicer measurement file is simplified
    # ------


    # plot the measurement of interest across all subjects and clusters
    fname = f'plot_all_data_{measurement}.pdf'
    vidx = list(header).index(measurement)
    plt.figure()
    for subject_measured, group in zip(measurement_list, group_id_list):
        value_list = subject_measured.measurement_matrix[:,vidx]
        if group == study_groups[0]:
            color = 'b'
        else:
            color = 'r'
        #plt.plot(np.sort(value_list), color)
        plt.plot(value_list, color+'.')
    plt.title(f'{measurement} measurements by group')
    plt.savefig(fname)
    plt.close()
    print(f"Saved data plot: {fname}")

    # get data by group to perform stats on measurement of interest
    vidx = list(header).index(measurement)
    data_0 = []
    data_1 = []
    for subject_measured, group in zip(measurement_list, group_id_list):
        #data = subject_measured.measurement_matrix[:,vidx]
        data = subject_measured.measurement_matrix[regions_for_stats,vidx]
        if group == study_groups[0]:
            data_0.append(data)
        elif group == study_groups[1]:
            data_1.append(data)
    data_0 = np.array(data_0)
    data_1 = np.array(data_1)

    # Statistical test in each cluster or input region
    print("Doing t-tests in each cluster or region.")
    p_val_list = []
    data_list = []
    for c in range(number_of_tests):
        #t, p = scipy.stats.ttest_ind(data_0[:,c], data_1[:,c])
        # ignore nan clusters for missing data. this also prevents nan p-values
        c_data_0 = data_0[:,c]
        c_data_1 = data_1[:,c]
        shape_before = np.array([c_data_0.shape[0], c_data_1.shape[0]])
        c_data_0 = c_data_0[~np.isnan(c_data_0)]
        c_data_1 = c_data_1[~np.isnan(c_data_1)]
        data_list.append(c_data_0)
        data_list.append(c_data_1)
        shape_after = np.array([c_data_0.shape[0], c_data_1.shape[0]])
        # warn if any cluster is totally absent
        if len(c_data_0) == 0 | len(c_data_1) == 0:
            print(f"Empty cluster across all subjects indicates data issue: {c}")
        # warn if empty clusters
        nan_count = shape_before - shape_after
        if np.sum(nan_count) > 0:
            print(f"\tCluster {c} : Warning. Empty/nan found in : {nan_count} subjects.")
        t, p = scipy.stats.ttest_ind(c_data_0, c_data_1)
        #print c_data_0.shape, c_data_1.shape
        if args.one_tailed:
            p_val_list.append(p / 2.0)
        else:
            p_val_list.append(p)

    uncorrected_significant = np.sum(np.array(p_val_list) < 0.05)
    print(f"Uncorrected: {uncorrected_significant} / {number_of_tests} : {100 * uncorrected_significant / float(number_of_tests)} %")

    ## if analyze_hierarchy:
    ##     for name, p_val in zip(names_for_stats, p_val_list):
    ##         print name, ":", p_val

    # bonferroni
    threshold = 0.05 / number_of_clusters
    corrected_significant = np.sum(np.array(p_val_list) < threshold)
    print(f"Bonferroni: {corrected_significant} / {number_of_tests} : {100 * corrected_significant / float(number_of_tests)} %")

    # FDR
    reject_null, corrected_p = statsmodels.sandbox.stats.multicomp.fdrcorrection0(np.array(p_val_list), alpha=fdr_q, method='indep')
    #reject_null, corrected_p = statsmodels.sandbox.stats.multicomp.fdrcorrection0(np.array(p_val_list), alpha=fdr_q, method='negcorr')
    corrected_significant = np.sum(reject_null)
    print(f"FDR at alpha/q = {fdr_q} : {corrected_significant} / {number_of_tests} : {100 * corrected_significant / float(number_of_tests)} %")

    ## # plot the data we tested, sorted by p-value
    ## cluster_order = np.argsort(np.array(p_val_list))
    ## plt.figure()
    ## plt.plot(data_1.T[cluster_order,:],'ro',markersize=2)
    ## plt.plot(data_0.T[cluster_order,:],'bo',markersize=2)
    ## plt.title(measurement)
    ## plt.savefig('tested_'+measurement+'.pdf')
    ## plt.close()

    fname = f'plot_region_means_{measurement}.pdf'
    print(f"Plotting mean values per cluster: {fname}")
    cluster_mean_0 = np.nanmean(data_0, axis=0)
    cluster_mean_1 = np.nanmean(data_1, axis=0)
    # Plot the mean values of the groups against each other in each cluster
    # and show which ones were significant
    #print "MEAN SHAPES:", cluster_mean_0.shape, cluster_mean_1.shape
    plt.figure()
    # plot the line if they were equal
    xmin = np.min(cluster_mean_0)
    xmax = np.max(cluster_mean_0)
    plt.plot([xmin,xmax],[xmin,xmax])
    # plot the means
    #markerfacecolor='none'
    #plt.plot(cluster_mean_0, cluster_mean_1, facecolors='none', edgecolors='k', markersize=3)
    plt.plot(cluster_mean_0, cluster_mean_1, 'o', markerfacecolor='none', markersize=3)
    plt.xlabel(study_groups[0])
    plt.ylabel(study_groups[1])
    # plot the significant ones on top
    #plt.plot(cluster_mean_0[reject_null], cluster_mean_1[reject_null],'ro',markersize=3)
    plt.plot(cluster_mean_0[reject_null], cluster_mean_1[reject_null],'ro',markersize=4)
    plt.title(f'{measurement} region means')
    plt.axis('equal')
    plt.savefig(fname)
    plt.close()

    fname = f'plot_boxplot_{measurement}.pdf'
    print(f"Plotting complete boxplot of all data: {fname}")
    # Get information to plot all significant and non-significant data in a boxplot.
    label_list = []
    text_list = []
    text_pos_list = []
    for p_val, region, c in zip(corrected_p, names_for_stats, list(range(number_of_tests))):
        #print region
        label_list.append(region[-4:]+'_0')
        label_list.append(region[-4:]+'_1')
        if p_val < fdr_q:
            text_list.append('*')
            text_pos_list.append(c*2 + 0.5)
        else:
            text_list.append('')
            text_pos_list.append(c*2 + 0.5)
    # Make the boxplot
    # try to make it big enough to fit all the data
    plt.figure(figsize=(20,20))
    x_pos_list = list(range(len(data_list)))
    #plt.boxplot(data_list, labels=label_list, showmeans=True, positions=x_pos_list)
    plt.boxplot(data_list, showmeans=True, positions=x_pos_list)
    plt.title(f'{measurement} all regions')
    ylim = plt.gca().get_ylim()[1]
    for text, pos in zip(text_list, text_pos_list):
        plt.text(pos, ylim*0.95, text, bbox=dict(facecolor='red', alpha=0.5))
    plt.xticks(x_pos_list, label_list, rotation='vertical')
    plt.savefig(fname)
    plt.close()

    fname = f'plot_boxplot_sig_{measurement}.pdf'
    print(f"Plotting significant boxplot: {fname}")
    # Now make a boxplot of only the significant ones
    # The data list length is number of clusters * number of groups, which must be 2 for now
    significant_data_list = []
    significant_labels_list = []
    significant_groups_list = []
    idx = 0
    for c, sig in zip(list(range(number_of_clusters)), reject_null):
        for g in study_groups:
            if sig:
                significant_data_list.append(data_list[idx])
                significant_labels_list.append(label_list[idx])
                significant_groups_list.append(g)
            idx += 1
    # try to make it big enough to fit all the data
    plt.figure(figsize=(20,20))
    if len(significant_data_list) > 0:
            positions = list(range(len(significant_data_list)))
            # regular boxplot
            plt.boxplot(significant_data_list, positions=positions)
            plt.xticks(positions, significant_labels_list, rotation='vertical')
            # also plot means in per-group color
            print(len(significant_data_list), len(significant_groups_list), len(positions))
            for d, g, p in zip(significant_data_list, significant_groups_list, positions):
                if g == study_groups[0]:
                    color = 'bo'
                else:
                    color ='ro'
                plt.plot(p, np.mean(d), color)
            # Pad margins so that markers don't get clipped by the axes
            plt.margins(0.2)
    plt.title(f'{measurement} significant regions')
    plt.savefig(fname)
    plt.close()

    # Choose to see original p values in MRML by uncommenting this line
    #output_p_val_list = p_val_list
    output_p_val_list = corrected_p

    if atlas_directory:
        # save a MRML file with tracts colored by p-value
        fname = f'./test_{measurement}.mrml'
        print("Saving MRML visualization:", fname)
        # find clusters in subject and atlas input directories
        input_mask = f"{atlas_directory}/cluster_*.vtp"
        atlas_clusters = sorted(glob.glob(input_mask))
        number_of_files = len(atlas_clusters)
        if number_of_files == number_of_clusters:
            colors = []
            significant_clusters = []
            for p_val, cluster in zip(output_p_val_list, atlas_clusters):
                if p_val < fdr_q:
                    #if p_val < 0.05:
                    #colors.append([(1-p_val)*255, p_val*(255/2), p_val*(255/2)])
                    colors.append([(1-10*p_val)*200+55, (1-10*p_val)*200+55, 0])
                    significant_clusters.append(cluster)
                ## elif p_val < 0.10:
                ##     colors.append([(1-p_val)*(255/2), p_val*(255/2), p_val*(255/2)])
                ##     significant_clusters.append(cluster)
                ## elif p_val < 0.15:
                ##     colors.append([50,50,50])
                ##     significant_clusters.append(cluster)
                #else:
                #colors.append([50,50,50])
            colors = np.array(colors)
            #print colors
            #wma.mrml.write(atlas_clusters, colors, './test_'+str(vidx)+'.mrml', ratio=1.0)
            wma.mrml.write(significant_clusters, colors, fname, ratio=1.0)
        else:
            print(f"Error: atlas directory and measurements have different cluster numbers: {number_of_files} {number_of_clusters}")


if __name__ == "__main__":
    main()

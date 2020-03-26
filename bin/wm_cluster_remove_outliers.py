#!/usr/bin/env python
import numpy
import argparse
import os
import glob
import shutil

import vtk

try:
    import whitematteranalysis as wma
except:
    print("<wm_cluster_from_atlas.py> Error importing white matter analysis package\n")
    raise

def main():
    #-----------------
    # Parse arguments
    #-----------------
    parser = argparse.ArgumentParser(
        description="Removes outliers in a subject dataset that was clustered from a cluster atlas. This script uses the atlas to identifies and remove outliers in each cluster of the subject. The atlas must be the same one used to cluster the subject dataset",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
    parser.add_argument("-v", "--version",
        action="version", default=argparse.SUPPRESS,
        version='1.0',
        help="Show program's version number and exit")
    parser.add_argument(
        'inputDirectory',
        help='A directory containing subject clusters (.vtp). Please note that file separator should not be provided at the end of the input so that the program can automatically recognize subject ID from the input folder name.')
    parser.add_argument(
        'atlasDirectory',
        help='The directory where the atlas is stored. Must contain atlas.p and atlas.vtp, as well as atlas clusters (.vtp)')
    parser.add_argument(
        'outputDirectory',
        help='The output directory will be created if it does not exist. Outlier-removed subject clusters will be stored in a subdirectory of the output directory. The subdirectory will have the same subject id as contained in the input directory.')
    parser.add_argument(
        '-verbose', action='store_true', dest="flag_verbose",
        help='Verbose. Run with -verbose for more text output.')
    parser.add_argument(
        '-cluster_outlier_std', action='store', dest="clusterOutlierStandardDeviation", type=float, default=2.0,
        help='After clustering, reject fiber outliers whose fiber probability (within their cluster) is more than this number of standard deviations below the mean. For more strict rejection, enter a smaller number such as 1.75. This probability is measured in an accurate way using pairwise comparison of all fibers in each cluster, to all fibers in the atlas. The purpose of this is to remove outlier fibers accurately within each cluster. This parameter can be tuned by the user depending on the amount of outlier fibers present in the tractography, and depending on how strictly outliers were removed in the atlas creation (which will affect estimation of the standard deviation of fiber probabilities in the atlas).')
    parser.add_argument(
        '-advanced_only_outlier_sigma', action='store', dest="outlierSigma", type=float, default=20.0,
        help='(Advanced parameter that probably should not be changed.) Local sigma used to compute fiber probability in cluster-based outlier removal. The default is 20mm. For stricter clustering, this may be reduced to 15mm.')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.inputDirectory):
        print("<wm_cluster_from_atlas.py> Error: Input subject directory", args.inputDirectory, "does not exist.")
        exit()
    
    if not os.path.isdir(args.atlasDirectory):
        print("<wm_cluster_from_atlas.py> Error: Atlas directory", args.atlasDirectory, "does not exist.")
        exit()
    
    outdir = args.outputDirectory
    if not os.path.exists(outdir):
        print("<wm_cluster_from_atlas.py> Output directory", outdir, "does not exist, creating it.")
        os.makedirs(outdir)
        
    subject_id = os.path.basename(args.inputDirectory)
    path_split = os.path.split(args.inputDirectory)
    if path_split[1] == '':
        subject_id = os.path.split(path_split[0])[1]
    else:
        subject_id = path_split[1]
    
    outdir = os.path.join(outdir, subject_id + '_outlier_removed')
    if not os.path.exists(outdir):
        print("<wm_cluster_from_atlas.py> Output directory", outdir, "does not exist, creating it.")
        os.makedirs(outdir)
    
    print("\n==========================")
    print("input directory:", args.inputDirectory)
    print("atlas directory:", args.atlasDirectory)
    print("output directory:", outdir)
    
    if args.flag_verbose:
        print("Verbose ON.")
    else:
        print("Verbose OFF.")
    verbose = args.flag_verbose
    
    print("==========================\n")
      
    # =======================================================================
    # Above this line is argument parsing. Below this line is the pipeline.
    # =======================================================================
    
    # Copy all MRML files to the new subject directory
    input_mask = "{0}/clustered_tracts*.mrml".format(args.inputDirectory)
    mrml_files = sorted(glob.glob(input_mask))
    for fname_src in mrml_files:
        fname_base = os.path.basename(fname_src)
        fname_dest = os.path.join(outdir, fname_base)
        shutil.copyfile(fname_src, fname_dest)
    
    # Overall progress/outlier removal info file
    fname_progress = os.path.join(outdir, 'cluster_log.txt')
    log_file = open(fname_progress, 'w')
    print('cluster','\t', 'input_fibers','\t', 'output_fibers','\t', 'number_removed','\t', 'percentage_removed','\t', 'mean_distance_before','\t', 'mean_distance_after','\t', 'mean_probability_before','\t', 'mean_probability_after', file=log_file)
    log_file.close()
    
    # read atlas
    print("<wm_cluster_from_atlas.py> Loading input atlas:", args.atlasDirectory)
    atlas = wma.cluster.load_atlas(args.atlasDirectory, 'atlas')
    bilateral = atlas.bilateral
    
    # To Do: test if can use stats saved in atlas, if they are the same as what I have computed here for the atlas cluster.
    ## # Add outlier information to the atlas
    ## atlas.cluster_outlier_std_threshold = cluster_outlier_std_threshold
    ## atlas.cluster_mean_similarity = cluster_mean_similarity
    ## atlas.cluster_std_similarity = cluster_std_similarity
    ## atlas.brain_mean_similarity = brain_mean_sim
    ## atlas.brain_std_similarity = brain_std_sim
    
    # find clusters in subject and atlas input directories
    input_mask = "{0}/cluster_*.vtp".format(args.atlasDirectory)
    atlas_clusters = sorted(glob.glob(input_mask))
    input_mask = "{0}/cluster_*.vtp".format(args.inputDirectory)
    subject_clusters = sorted(glob.glob(input_mask))
    
    # check these lists are the same length
    number_ac = len(atlas_clusters)
    number_sc = len(subject_clusters)
    
    print("Number of atlas and subject clusters:", number_ac, number_sc)
    
    if number_ac != number_sc:
        print("<wm_cluster_from_atlas.py> Error: Cluster number mismatch. \nAtlas directory (", args.atlasDirectory, ") has", number_ac, "clusters but subject directory (", args.inputDirectory, ") has", number_sc, "clusters.")
        exit()
    
    
    # Loop over all clusters in atlas and subject. Remove outliers in the subject cluster using the same heuristics as during atlas creation.
    print("Starting local cluster outlier removal")
    distance_method = 'Mean'
    # 20mm works well for cluster-specific outlier removal in conjunction with the mean distance.
    cluster_local_sigma = args.outlierSigma
    print("Local sigma for cluster outlier removal:", cluster_local_sigma)
    
    cluster_outlier_std_threshold = args.clusterOutlierStandardDeviation
    print("Standard deviation for fiber outlier removal after clustering (accurate local probability): ", cluster_outlier_std_threshold)
    
    c = 1
    for ca, cs in zip(atlas_clusters,subject_clusters):
        # grab the cluster as polydata
        pd_atlas = wma.io.read_polydata(ca)
        pd_subject = wma.io.read_polydata(cs)
        pd_out_fname = os.path.join(outdir, os.path.basename(cs))
        
        number_fibers_in_atlas_cluster = pd_atlas.GetNumberOfLines()
        number_fibers_in_subject_cluster = pd_subject.GetNumberOfLines()
    
        # if either cluster is empty, we have to skip outlier removal
        if number_fibers_in_subject_cluster == 0:
            wma.io.write_polydata(pd_subject, pd_out_fname)
            print("cluster", c, "empty in subject")
            log_file = open(fname_progress, 'a')
            print(c,'\t', number_fibers_in_subject_cluster,'\t', 0,'\t', 0,'\t', 0, '\t', 0, '\t', 0, '\t', 0, '\t', 0, file=log_file)
            log_file.close()
            c += 1
            continue
    
        if number_fibers_in_atlas_cluster == 0:
            # this should not ever happen. Note that "outlier removed" atlas is temporary and should not be used for classification.
            # (Note the other option is to remove this cluster as it was removed at the end of the iteration for the atlas.)
            print("cluster", c, "empty in atlas")
            print("ERROR: An atlas should not contain empty clusters. Please use the initial_clusters atlas for this outlier removal script and for subject clustering.")
            log_file = open(fname_progress, 'a')
            print(c,'\t', number_fibers_in_subject_cluster,'\t', 0,'\t', 0,'\t', 0, '\t', 0, '\t', 0, '\t', 0, '\t', 0, file=log_file)
            log_file.close()
            c += 1
            continue
        
        # for the future, in case this can be integrated here and use atlas information for hemisphere separation 
        ## # figure out hemisphere labels for this cluster
        ## farray = wma.fibers.FiberArray()
        ## farray.hemispheres = True
        ## farray.hemisphere_percent_threshold = 0.90
        ## farray.convert_from_polydata(pd_subject, points_per_fiber=50)
        ## fiber_hemisphere[fiber_indices] = farray.fiber_hemisphere
        ## cluster_left_hem.append(farray.number_left_hem)
        ## cluster_right_hem.append(farray.number_right_hem)
        ## cluster_commissure.append(farray.number_commissure)
    
        # Compute distances and fiber probabilities FOR ATLAS to form a basis for comparison
        # Note that ideally this could be done once and stored, rather than computed for each subject
        if distance_method == 'StrictSimilarity':
            cluster_distances = wma.cluster._pairwise_distance_matrix(pd_atlas, 0.0, number_of_jobs=1, bilateral=bilateral, distance_method=distance_method, sigmasq = cluster_local_sigma * cluster_local_sigma)
            cluster_similarity_atlas = cluster_distances
        else:
            cluster_distances = wma.cluster._pairwise_distance_matrix(pd_atlas, 0.0, number_of_jobs=1, bilateral=bilateral, distance_method=distance_method)
            cluster_similarity_atlas = wma.similarity.distance_to_similarity(cluster_distances, cluster_local_sigma * cluster_local_sigma)
                
        # Compute distances from SUBJECT fibers to atlas fibers, and fiber probabilities
        if distance_method == 'StrictSimilarity':
            cluster_similarity_subject = wma.cluster._rectangular_distance_matrix(pd_subject, pd_atlas, 0.0, number_of_jobs=1, bilateral=bilateral, distance_method=distance_method, sigmasq = cluster_local_sigma * cluster_local_sigma)
        else:
            cluster_distances = wma.cluster._rectangular_distance_matrix(pd_subject, pd_atlas, 0.0, number_of_jobs=1, bilateral=bilateral, distance_method=distance_method)
            cluster_similarity_subject = wma.similarity.distance_to_similarity(cluster_distances, cluster_local_sigma * cluster_local_sigma)
    
        #print "SHAPE cluster similarity:",  cluster_similarity_atlas.shape, cluster_similarity_subject.shape
        
        #p(f1) = sum over all f2 of p(f1|f2) * p(f2)
        # by using sample we estimate expected value of the above
        #  get total similarity (probability) based on the atlas, and normalize by the number of fibers used in the comparison for easier comparison across clusters
        total_similarity_atlas = (numpy.sum(cluster_similarity_atlas, axis=1) - 1.0) / number_fibers_in_atlas_cluster
        total_similarity_subject = numpy.sum(cluster_similarity_subject, axis=0)/ number_fibers_in_atlas_cluster
    
        #print "SHAPE total similarity:", total_similarity_atlas.shape, total_similarity_subject.shape
    
        if verbose:
            print("cluster", c, "tsim_atlas:", numpy.min(total_similarity_atlas), numpy.mean(total_similarity_atlas), numpy.max(total_similarity_atlas), "num fibers atlas:", number_fibers_in_atlas_cluster)
            print("cluster", c, "tsim_subject:", numpy.min(total_similarity_subject), numpy.mean(total_similarity_subject), numpy.max(total_similarity_subject), "num fibers subject:", number_fibers_in_subject_cluster)
        
        # remove outliers with low similarity to their cluster
        mean_sim_atlas = numpy.mean(total_similarity_atlas)
        cluster_std_atlas = numpy.std(total_similarity_atlas)
        cutoff = mean_sim_atlas - cluster_outlier_std_threshold*cluster_std_atlas
        fiber_indices = list(range(number_fibers_in_subject_cluster))
        #print "LEN INDICES:", len(fiber_indices)
        reject_idx = []
        
        for (fidx, sim) in zip(fiber_indices, total_similarity_subject):
            if sim < cutoff:
                reject_idx.append(fidx)
    
        number_rejected =  len(reject_idx)
        
        print("cluster", c, "rejecting cluster outlier fibers:", number_rejected, "/", number_fibers_in_subject_cluster, "=", float(number_rejected)/number_fibers_in_subject_cluster)
    
        # Output a new polydata with the outliers removed.
        mask = numpy.ones([number_fibers_in_subject_cluster])
        mask[reject_idx] = 0
        #print mask, number_fibers_in_subject_cluster, mask.shape, pd_subject.GetNumberOfLines()
        pd_c = wma.filter.mask(pd_subject, mask, verbose=False, preserve_point_data=True, preserve_cell_data=True)
        wma.io.write_polydata(pd_c, pd_out_fname)
    
        cluster_distances = numpy.sqrt(cluster_distances)
        cluster_mean_distances = numpy.mean(cluster_distances, axis=0)
        mask = mask == 1
    
        log_file = open(fname_progress, 'a')
        print(c,'\t', number_fibers_in_subject_cluster,'\t', number_fibers_in_subject_cluster-number_rejected,'\t', number_rejected,'\t', float(number_rejected)/number_fibers_in_subject_cluster, '\t', numpy.mean(cluster_mean_distances), '\t', numpy.mean(cluster_mean_distances[mask]), '\t', numpy.mean(total_similarity_subject), '\t', numpy.mean(total_similarity_subject[mask]), file=log_file)
        log_file.close()
    
        c += 1

if __name__ == '__main__':
    main()

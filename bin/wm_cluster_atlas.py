#!/usr/bin/env python
import numpy
import argparse
import os
import multiprocessing
import vtk
import time

try:
    import whitematteranalysis as wma
except:
    print "<wm_cluster_atlas.py> Error importing white matter analysis package\n"
    raise


#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Runs clustering of tractography for multiple subjects to create an atlas. This tract atlas can then be applied to the complete set of fibers from individual subjects of interest. To make it possible to cluster the high number of fibers, this code uses random sampling and the Nystrom method as described in the reference below.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputDirectory',
    help='A directory of (already registered) whole-brain tractography as vtkPolyData (.vtk or .vtp).')
parser.add_argument(
    'outputDirectory',
    help='The output directory will be created if it does not exist.')
parser.add_argument(
    '-f', action="store", dest="numberOfFibers", type=int,
    help='Number of fibers to analyze from each subject. Default is 2000 (this assumes 10 or more subjects; for fewer subjects, more fibers are needed per subject).')
parser.add_argument(
    '-l', action="store", dest="fiberLength", type=int,
    help='Minimum length (in mm) of fibers to analyze. 25mm is default. This default is reasonable for single-tensor DTI tractography. For two-tensor UKF, 80mm is more reasonable (as tracts are longer in general). Run the quality control script on your data first and inspect the fiber length distribution in your dataset. Note that with too low a threshold, the clustering will be dominated by the more prevalent short fibers, such as u-fibers, instead of longer association fibers.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')
parser.add_argument(
    '-verbose', action='store_true', dest="flag_verbose",
    help='Verbose. Run with -verbose for more text output in the terminal window.')
parser.add_argument(
    '-k', action="store", dest="numberOfClusters", type=int,
    help='Number of clusters to find. Default is 250, which is reasonable for single-tensor DTI tractography. Useful range is from 200 to 600+. For two-tensor UKF, 400 or more is a reasonable number.')
parser.add_argument(
    '-thresh', action="store", dest="distanceThreshold", type=float,
    help='Threshold (in mm) below which fiber points are considered in the same position. Default is 2mm for cross-subject atlas clustering. This helps find correspondences across subjects by ignoring very small differences.')
parser.add_argument(
    '-nystrom_sample', action="store", dest="sizeOfNystromSample", type=int,
    help='Number of fibers to use in the Nystrom sample. This is the random sample of fibers to which every input fiber is compared, to structure the clustering problem. 2000 is default. Must be >1500. Increase for larger datasets. Reduce to limit memory use. One rule of thumb is that this number should be similar to the number of fibers used per subject.')
parser.add_argument(
    '-sigma', action="store", dest="sigma", type=float,
    help='Sigma for kernel. Controls distance over which fibers are considered similar. 60mm is default. Reduce for stricter clustering with ample data, or for physically smaller problems like a subset of the brain.')
parser.add_argument(
    '-mrml_fibers', action="store", dest="showNFibersInSlicer", type=float,
    help='Approximate upper limit on number of fibers to show when MRML scene of clusters is loaded into slicer. Default is 10000 fibers; increase for computers with more memory. Note this can be edited later in the MRML file by searching for SubsamplingRatio and editing that number throughout the file. Be sure to use a text editor program (save as plain text format). An extra MRML file will be saved for visualizing 100%% of fibers.')
#parser.add_argument(
#    '-remove_outliers', action='store_true', dest="flag_remove_outliers",
#    help='Define outlier clusters (default, ones with less than 90 percent of subjects present). These will be segmented separately.')
#parser.add_argument(
#    '-subject_percent', action="store", dest="subjectPercent", type=float,
#    help='Threshold for defining good vs outlier clusters. Default is 90 percent of subjects must be present to keep a cluster. Lower this to 70 or 80 if data are more variable or the number of clusters is very high.')
parser.add_argument(
    '-bilateral_off', action='store_true', dest="flag_bilateral_off",
    help='Turn off bilateral clustering. In general, anatomy is better and more stably represented with bilateral clusters, so that is the default. The bilateral clusters can be split at the midline later for analyses.')
parser.add_argument(
    '-outlier_std', action='store', dest="outlierStandardDeviation",type=float,
    help='Reject fiber outliers whose total fiber similarity is more than this number of standard deviations below the mean. Default is 2.0 standard deviations. For more strict rejection, enter a smaller number such as 1.75. To turn off outlier rejection, enter a large number such as 100.')
parser.add_argument(
    '-norender', action='store_true', dest="flag_norender",
    help='No Render. Prevents rendering of images that would require an X connection.')
parser.add_argument(
    '-advanced_only_testing_distance', action="store", dest="distanceMethod", type=str,
    help='(Advanced parameter for testing only.) Distance method for pairwise fiber comparison. Default is Mean, which is the average distance between points on the fibers. Other options are Hausdorff (the worst-case distance), StrictSimilarity (multiplies all pointwise similarities along fiber).')
parser.add_argument(
    '-advanced_only_debug_nystrom_off', action='store_true', dest="flag_nystrom_off",
    help='(Advanced parameter for testing only.) Turn off the Nystrom method, e.g. perform clustering by computing the complete distance matrix of pairwise distances of all input fibers. This will not create an output atlas. It is for small datasets and code testing only.')
parser.add_argument(
    '-advanced_only_testing_on', action='store_true', dest="flag_testing_on",
    help='(Advanced parameter for testing only.) Turn on saving of additional output for code testing only.')
parser.add_argument(
    '-advanced_only_random_seed', action='store', dest="randomSeed", type=int,
    help='(Advanced parameter for testing only.) Set random seed for reproducible clustering in software tests.')
parser.add_argument(
    '-advanced_only_force_pos_def_off', action='store_true', dest="flag_pos_def_off",
    help='(Advanced parameter for testing only.) Turn off replacing A matrix by a close positive definite matrix.')
args = parser.parse_args()


if not os.path.isdir(args.inputDirectory):
    print "<wm_cluster_atlas.py> Error: Input directory", args.inputDirectory, "does not exist or is not a directory."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "<wm_cluster_atlas.py> Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)
    
print "\n=========================="
print "<wm_cluster_atlas.py> Clustering parameters"
print "input directory:\n", args.inputDirectory
print "output directory:\n", args.outputDirectory

if args.numberOfFibers is not None:
    print "fibers to analyze per subject: ", args.numberOfFibers
    number_of_fibers_per_subject = args.numberOfFibers
else:
    number_of_fibers_per_subject = 2000
    print "fibers to analyze per subject: Setting to default", number_of_fibers_per_subject

if args.fiberLength is not None:
    fiber_length = args.fiberLength
else:
    fiber_length = 25.0
print "minimum length of fibers to analyze (in mm): ", args.fiberLength

if args.numberOfJobs is not None:
    number_of_jobs = args.numberOfJobs
else:
    # default to 1 job for usage on compute clusters. Also, the
    # multiprocessing is not used efficiently in our code and should
    # be avoided for large clustering problems, for now.
    number_of_jobs = 1
print 'Using N jobs:', number_of_jobs

if args.flag_verbose:
    print "Verbose ON."
else:
    print "Verbose OFF."
verbose = args.flag_verbose

if args.flag_norender:
    print "No rendering (for compute servers without X connection)."
else:
    print "Rendering. After clustering, will create colorful jpg images of the group."
render = not args.flag_norender

print "RENDER:", render
if args.numberOfClusters is not None:
    number_of_clusters = args.numberOfClusters
else:
    number_of_clusters = 250
print "Number of clusters to find: ", number_of_clusters

if args.outlierStandardDeviation is not None:
    outlier_std_threshold = args.outlierStandardDeviation 
else:
    outlier_std_threshold = 2.0
print "Standard deviation for fiber outlier removal: ", outlier_std_threshold

if args.sizeOfNystromSample is not None:
    number_of_sampled_fibers = args.sizeOfNystromSample
else:
    number_of_sampled_fibers = 2000
print "Size of Nystrom sample: ", number_of_sampled_fibers

if args.sigma is not None:
    sigma = args.sigma
else:
    sigma = 60
print "Sigma in mm: ", sigma

if args.showNFibersInSlicer is not None:
    number_of_fibers_to_display = args.showNFibersInSlicer
else:
    number_of_fibers_to_display = 10000.0
print "Maximum total number of fibers to display in MRML/Slicer: ", number_of_fibers_to_display

#if args.flag_remove_outliers:
#    if args.subjectPercent is not None:
#        fraction_to_keep_cluster = args.subjectPercent / 100.0
#    else:
#        fraction_to_keep_cluster = 0.9
#    print "Separating outlier clusters with fewer than: ", fraction_to_keep_cluster * 100.0, "percent of subjects."
#else:
#    print "Not removing outlier clusters."

if args.flag_bilateral_off:
    bilateral = False
    print "Bilateral clustering OFF."
else:
    bilateral = True
    print "Bilateral clustering ON."

if args.distanceThreshold is not None:
    threshold = args.distanceThreshold
else:
    # for across-subjects matching
    threshold = 2.0
print "Threshold (in mm) for fiber distances: ", threshold


if args.distanceMethod is not None:
    distance_method = args.distanceMethod
else:
    # for across-subjects matching
    distance_method = 'Mean'
print "Fiber distance or comparison method: ", distance_method

if args.flag_nystrom_off:
    use_nystrom = False
    print "Nystrom method is OFF (for testing of software with small datasets only)."
else:
    use_nystrom=True
    print "Nystrom method is ON (default)."

if args.flag_testing_on:
    testing = True
    print "Saving pickle files for software testing purposes."
else:
    testing = False

if args.randomSeed is not None:
    print "Setting random seed to: ", args.randomSeed
random_seed = args.randomSeed

if args.flag_pos_def_off:
    print "Positive definite approximation for A is OFF (testing only)."
else:
    print "Positive definite approximation for A is ON (default)."
    
pos_def_approx = ~args.flag_pos_def_off

# default clustering parameters that probably don't need to be changed
# from TMI 2007 paper
#use_nystrom=True
#distance_method = 'Mean'
use_normalized_cuts = True
number_of_eigenvectors = 10


# Print input parameter information to a file as well as to the terminal
#print args
#f = open(os.path.join(outdir, 'parameters.txt'), 'w+')
#print >> f, args
#f.close()

# another option. was not used in TMI 2007 paper. would need a different sigma.
#distance_method ='Hausdorff'
# This was used in the TMI paper but 10 eigenvectors
# contain almost as much information and reduce noise for single subject clustering
#number_of_eigenvectors = 20

# =======================================================================
# Above this line is argument parsing. Below this line is the pipeline.
# =======================================================================

input_polydatas = wma.io.list_vtk_files(args.inputDirectory)
number_of_subjects = len(input_polydatas)

print "<wm_cluster_atlas.py> Found ", number_of_subjects, "subjects in input directory:", args.inputDirectory
if number_of_subjects < 1:
    print "\n<wm_cluster_atlas.py> Error: No .vtk or .vtp files were found in the input directory.\n"
    exit()

total_number_of_fibers = number_of_fibers_per_subject * number_of_subjects

print "Input number of subjects (number of vtk/vtp files): ", number_of_subjects
print "==========================\n"
print "<wm_cluster_atlas.py> Starting file I/O and computation."

# output summary file to save information about what was run
readme_fname = os.path.join(outdir, 'README.txt')
readme_file = open(readme_fname, 'w')
outstr = "Group (Atlas) Clustering Summary\n"
outstr += '----------------------\n'
outstr += '\n'
outstr += "Input Directory: "
outstr += args.inputDirectory
outstr += '\n'
outstr += "Output Directory: "
outstr += args.outputDirectory
outstr += '\n'
outstr += "Number of Subjects: "
outstr += str(number_of_subjects)
outstr += '\n'
outstr += '\n'
outstr +=  "Current date: "  + time.strftime("%x")
outstr += '\n'
outstr +=  "Current time: " + time.strftime("%X")
outstr += '\n'
outstr += '\n'
outstr += "Path to Script: " + os.path.realpath(__file__)
outstr += '\n'
outstr += "Working Directory: " + os.getcwd()
outstr += '\n'
outstr += '\n'
outstr += "Description of Outputs\n"
outstr += '---------------------\n'
outstr += 'atlas.p, atlas.vtp: The information needed to cluster a new subject.\n'
outstr += 'cluster_*.vtp: Clustered fibers from all subjects for visualization.\n'
outstr += 'clustered_tracts.mrml:  Slicer scene for loading all fiber files in color.\n'
outstr += 'cluster_quality_control.txt: Measures from each cluster including mean and variability across subjects.\n'
outstr += 'input_subjects.txt:  List of subject index, ID, and full path to input file.\n'
outstr += 'README.txt:  This summary file.\n'
outstr += 'subjects_per_cluster_hist.pdf:  Histogram showing the number of subjects present in each cluster. Ideally, most clusters should contain all subjects.\n'
outstr += 'view_*.jpg: Images of the clustered brains for visual quality control. Colors should be bright and look related to the anatomy.\n'
outstr += '\n'
outstr += '\n'
outstr += "Command Line Arguments\n"
outstr += '----------------------\n'
outstr += str(args)
outstr += '\n'
outstr += '\n'
outstr += "Input Fiber Files\n"
outstr += '-----------------\n'
for pd in input_polydatas:
    outstr += pd
    outstr += '\n'
readme_file.write(outstr)
readme_file.close()


# read in data
input_pds = list()
for fname in input_polydatas:
    # read data
    print "<wm_cluster_atlas.py> Reading input file:", fname
    pd = wma.io.read_polydata(fname)
    # preprocessing step: minimum length
    #print "<wm_cluster_atlas.py> Preprocessing by length:", fiber_length, "mm."
    pd2 = wma.filter.preprocess(pd, fiber_length,verbose=verbose)
    # preprocessing step: fibers to analyze
    if number_of_fibers_per_subject is not None:
        print "<wm_cluster_atlas.py> Downsampling to", number_of_fibers_per_subject, "fibers from",  pd2.GetNumberOfLines(),"fibers over length", fiber_length, "."
        pd3 = wma.filter.downsample(pd2, number_of_fibers_per_subject, verbose=verbose, random_seed=random_seed)
        if pd3.GetNumberOfLines() != number_of_fibers_per_subject:
            print "<wm_cluster_atlas.py> Fibers found:", pd3.GetNumberOfLines(), "Fibers requested:", number_of_fibers_per_subject
            print "\n<wm_cluster_atlas.py> ERROR: too few fibers over length threshold in subject:", fname
            exit()
    else:
        pd3 = pd2
    input_pds.append(pd3)
    del pd
    del pd2
    # safe because list has a reference to pd3
    del pd3

# append into one polydata object for clustering
appender = vtk.vtkAppendPolyData()
for pd in input_pds:
    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        appender.AddInputData(pd)
    else:
        appender.AddInput(pd)
appender.Update()
input_data = appender.GetOutput()
del input_pds

# figure out which subject each fiber was from in the input to the clustering
subject_fiber_list = list()
for sidx in range(number_of_subjects):
    for fidx in range(number_of_fibers_per_subject):
        subject_fiber_list.append(sidx)
subject_fiber_list = numpy.array(subject_fiber_list)

#-----------------
# Run clustering
#-----------------

# Check there are enough fibers for requested analysis
if number_of_sampled_fibers >= input_data.GetNumberOfLines():
    print "<wm_cluster_atlas.py> Error: Nystrom sample size is larger than number of fibers available."
    print "number_of_subjects:", number_of_subjects
    print "number_of_fibers_per_subject:", number_of_fibers_per_subject
    print "total_number_of_fibers:", total_number_of_fibers
    print "requested subsample size from above total:", number_of_sampled_fibers
    exit()

# Set up random seed
if random_seed is not None:
    print "<wm_cluster_atlas.py> Setting random seed to", random_seed
    numpy.random.seed(seed=random_seed)

# Run clustering several times, removing outliers each time.
#cluster_iterations = 5
cluster_iterations = 10
# almost no effect
#cluster_std_threshold = 2.0
#cluster_std_threshold = outlier_std_threshold
# too much?
#cluster_std_threshold = 0.5
#cluster_std_threshold = 1.0
# compute pairwise distances and affinities with local sigma for outlier removal
#cluster_local_sigma = 15
#cluster_local_sigma = 30
cluster_local_sigmas = [30, 20, 15, 10, 5, 5, 5, 5, 5, 5, 5, 5, 5]
#cluster_std_thresholds = [2.0, 2.0, 2.5, 2.5, 3.0, ]
cluster_std_threshold = 2.0
        
for iter in range(cluster_iterations):
    cluster_local_sigma = cluster_local_sigmas[iter]
    # make a directory for the current iteration
    dirname = "iteration_%05d" % (iter)
    outdir_current = os.path.join(outdir, dirname)

    # Calculate indices of random sample for Nystrom method
    nystrom_mask = numpy.random.permutation(input_data.GetNumberOfLines()) < number_of_sampled_fibers

    print "BEFORE cluster: Polydata size:", input_data.GetNumberOfLines(), "Subject list for fibers:", subject_fiber_list.shape

    # Run clustering on the polydata
    print '<wm_cluster_atlas.py> Starting clustering.'
    output_polydata_s, cluster_numbers_s, color, embed, distortion, atlas, reject_idx = \
        wma.cluster.spectral(input_data, number_of_clusters=number_of_clusters, \
                                 number_of_jobs=number_of_jobs, use_nystrom=use_nystrom, \
                                 nystrom_mask = nystrom_mask, \
                                 number_of_eigenvectors=number_of_eigenvectors, \
                                 sigma=sigma, distance_method=distance_method, \
                                 normalized_cuts=use_normalized_cuts, threshold=threshold, \
                                 outlier_std_threshold=outlier_std_threshold, \
                                 pos_def_approx=pos_def_approx, \
                                 bilateral=bilateral)

    # If any fibers were rejected, delete the corresponding entry in this list
    subject_fiber_list = numpy.delete(subject_fiber_list, reject_idx)
    print "After cluster: Polydata size:", output_polydata_s.GetNumberOfLines(), "Subject list for fibers:", subject_fiber_list.shape
    
    # Save the output in our atlas format for automatic labeling of full brain datasets.
    # This is the data used to label a new subject.
    # Also write the polydata with cluster indices saved as cell data. This is a one-file output option for clusters.
    # Finally, save some quality control metrics and save the atlas clusters as individual polydatas. This is used to 
    # set up a mrml hierarchy file and to visualize the output in Slicer. This data is not used to label
    # a new subject.
    
    outdir1 = os.path.join(outdir_current, 'initial_clusters')
    if not os.path.exists(outdir1):
        print "<wm_cluster_atlas.py> Output directory", outdir1, "does not exist, creating it."
        os.makedirs(outdir1)    
    print '<wm_cluster_atlas.py> Saving output files in directory:', outdir1
    wma.cluster.output_and_quality_control_cluster_atlas(atlas, output_polydata_s, subject_fiber_list, input_polydatas, number_of_subjects, outdir1, cluster_numbers_s, color, embed, number_of_fibers_to_display, testing=testing, verbose=verbose, render_images=render)

    # Remove outliers from this iteration and save atlas again                                                 
    print "Starting local cluster outlier removal"

    # outlier cluster index is going to be -1                                                     
    reject_idx = list() 
    cluster_indices = range(atlas.centroids.shape[0])
    fiber_mean_sim = numpy.zeros(cluster_numbers_s.shape)

    for c in cluster_indices:
        mask = cluster_numbers_s == c

        # reject the whole cluster if there are few subjects (tractography error)
        subjects_threshold = 5
        subjects_per_cluster = len(set(subject_fiber_list[mask]))

        # grab the cluster as polydata
        pd_c = wma.filter.mask(output_polydata_s, mask,verbose=False)
        cluster_distances = wma.cluster._pairwise_distance_matrix(pd_c, 0.0, number_of_jobs=1, bilateral=bilateral, distance_method=distance_method)
        cluster_similarity = wma.similarity.distance_to_similarity(cluster_distances, cluster_local_sigma * cluster_local_sigma)

        tmp = numpy.sqrt(cluster_distances)
        if verbose:
            print "cluster", c, "dist:", numpy.min(tmp), numpy.mean(tmp), numpy.max(tmp)
            print "cluster", c, "sim :", numpy.min(cluster_similarity), numpy.mean(cluster_similarity), numpy.max(cluster_similarity)

        # get total similarity for each fiber
        total_similarity = numpy.sum(cluster_similarity, axis=1) - 1.0
        if verbose:
            print "cluster", c, "tsim:", numpy.min(total_similarity), numpy.mean(total_similarity), numpy.max(total_similarity), "num fibers:", numpy.sum(mask), "num subjects:", subjects_per_cluster

        # remove outliers with low similarity to their cluster
        number_fibers_in_cluster = numpy.sum(mask)
        mean_sim = numpy.mean(total_similarity)
        cluster_std = numpy.std(total_similarity)
        count = 0
        for fidx in range(len(cluster_numbers_s)):
            if cluster_numbers_s[fidx] == c:
                fiber_mean_sim[fidx] = total_similarity[count] / number_fibers_in_cluster
                if total_similarity[count] < mean_sim - cluster_std_threshold*cluster_std:
                    #print fidx, cluster_numbers_s[fidx], dists[count]
                    cluster_numbers_s[fidx] = -1
                    reject_idx.append(fidx)
                elif subjects_per_cluster < subjects_threshold:
                    cluster_numbers_s[fidx] = -1
                    reject_idx.append(fidx)
                count += 1
    print "Rejecting cluster outlier fibers:", len(reject_idx)

    # in a second pass, also remove outliers whose average fiber
    # similarity to their cluster is too low compared to the whole brain.
    # This can prune fibers from variable clusters that might be missed above
    brain_mean_sim = numpy.mean(fiber_mean_sim)
    brain_std_sim = numpy.std(fiber_mean_sim)

    for fidx in range(len(cluster_numbers_s)):
        if fiber_mean_sim[fidx] < brain_mean_sim - cluster_std_threshold*brain_std_sim:
            #print fidx, cluster_numbers_s[fidx], fiber_mean_sim[count]
            if cluster_numbers_s[fidx] >= 0:
                cluster_numbers_s[fidx] = -1
                reject_idx.append(fidx)

    reject_idx = numpy.array(reject_idx)

    print "Rejecting whole-brain cluster outlier fibers:", len(reject_idx)

    outdir2 = os.path.join(outdir_current, 'remove_outliers')
    if not os.path.exists(outdir2):
        print "<wm_cluster_atlas.py> Output directory", outdir2, "does not exist, creating it."
        os.makedirs(outdir2)
    print '<wm_cluster_atlas.py> Saving output files in directory:', outdir2

    mask = cluster_numbers_s == -1
    fname_c = os.path.join(outdir2, 'outlier.vtp')
    pd_c = wma.filter.mask(output_polydata_s, mask,verbose=False)
    wma.io.write_polydata(pd_c, fname_c)

    print "Before save Polydata size:", output_polydata_s.GetNumberOfLines(), "Subject list for fibers:", subject_fiber_list.shape

    wma.cluster.output_and_quality_control_cluster_atlas(atlas, output_polydata_s, subject_fiber_list, input_polydatas, number_of_subjects, outdir2, cluster_numbers_s, color, embed, number_of_fibers_to_display, testing=testing, verbose=verbose, render_images=render)

    test = subject_fiber_list[numpy.nonzero(mask)]

    # Remove outliers for the next iteration
    # If any fibers were rejected, delete the corresponding entry in this list
    subject_fiber_list = numpy.delete(subject_fiber_list, reject_idx)
    # Also delete the fiber in the polydata
    mask = cluster_numbers_s >= 0

    input_data = wma.filter.mask(output_polydata_s, mask, verbose=verbose)

    print "End iter Polydata size:", input_data.GetNumberOfLines(), "Subject list for fibers:", subject_fiber_list.shape, "TEST:", test.shape

print "==========================\n"
print '<wm_cluster_atlas.py> Done clustering atlas. See output in directory:\n ', outdir, '\n'



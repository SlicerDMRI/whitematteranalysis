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

HAVE_PLT = 1
try:
    import matplotlib.pyplot as plt
except:
    print "<wm_cluster_atlas.py> Error importing matplotlib.pyplot package, can't plot quality control data.\n"
    HAVE_PLT = 0    

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
    help='Verbose. Run with -verbose for more text output.')
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
    help='Approximate upper limit on number of fibers to show when MRML scene of clusters is loaded into slicer. Default is 10000 fibers; increase for computers with more memory. Note this can be edited later in the MRML file by searching for SubsamplingRatio and editing that number throughout the file. Be sure to use a text editor program (save as plain text format). An extra MRML file will be saved for visualizing 100% of fibers.')
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
    '-advanced_only_testing_distance', action="store", dest="distanceMethod", type=str,
    help='(Advanced parameter for testing only.) Distance method for pairwise fiber comparison. Default is Mean, which is the average distance between points on the fibers. Other options are Hausdorff (the worst-case distance), StrictSimilarity (multiplies all pointwise similarities along fiber).')
parser.add_argument(
    '-advanced_only_debug_nystrom_off', action='store_true', dest="flag_nystrom_off",
    help='(Advanced parameter for testing only.) Turn off the Nystrom method, e.g. perform clustering by computing the complete distance matrix of pairwise distances of all input fibers. This will not create an output atlas. It is for small datasets and code testing only.')
    
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
    print 'CPUs detected:', multiprocessing.cpu_count()
    number_of_jobs = multiprocessing.cpu_count()
print 'Using N jobs:', number_of_jobs

if args.flag_verbose:
    print "Verbose ON."
else:
    print "Verbose OFF."
verbose = args.flag_verbose
if args.numberOfClusters is not None:
    number_of_clusters = args.numberOfClusters
else:
    number_of_clusters = 250
print "Number of clusters to find: ", number_of_clusters

if args.distanceThreshold is not None:
    threshold = args.distanceThreshold
else:
    # for across-subjects matching
    threshold = 2.0
print "Threshold (in mm) for fiber distances: ", threshold

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
    show_fibers = args.showNFibersInSlicer
else:
    show_fibers = 10000.0
print "Maximum total number of fibers to display in MRML/Slicer: ", show_fibers

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
outstr += 'view_*.png: Images of the clustered brains for visual quality control. Colors should be bright and look related to the anatomy.\n'
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
    pd2 = wma.filter.preprocess(pd, fiber_length,verbose=False)
    # preprocessing step: fibers to analyze
    if number_of_fibers_per_subject is not None:
        print "<wm_cluster_atlas.py> Downsampling to ", number_of_fibers_per_subject, "fibers."
        pd3 = wma.filter.downsample(pd2, number_of_fibers_per_subject,verbose=False)
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

#-----------------
# Run clustering
#-----------------

# Check there are enough fibers for requested analysis
if number_of_sampled_fibers >= input_data.GetNumberOfLines():
    print "<wm_cluster_atlas.py>Error: Nystrom sample size is larger than number of fibers available."
    print "number_of_subjects:", number_of_subjects
    print "number_of_fibers_per_subject:", number_of_fibers_per_subject
    print "total_number_of_fibers:", total_number_of_fibers
    print "requested subsample size from above total:", number_of_sampled_fibers
    exit()

# Calculate indices of random sample for Nystrom method
nystrom_mask = numpy.random.permutation(input_data.GetNumberOfLines()) < number_of_sampled_fibers

# Run clustering on the polydata
output_polydata_s, cluster_numbers_s, color, embed, distortion, atlas = \
    wma.cluster.spectral(input_data, number_of_clusters=number_of_clusters, \
                             number_of_jobs=number_of_jobs, use_nystrom=use_nystrom, \
                             nystrom_mask = nystrom_mask, \
                             number_of_eigenvectors=number_of_eigenvectors, \
                             sigma=sigma, distance_method=distance_method, \
                             normalized_cuts=use_normalized_cuts, threshold=threshold, \
                             bilateral=bilateral)

print '<wm_cluster_atlas.py>Saving output files in directory:', outdir

# Save the output in our atlas format for automatic labeling of full brain datasets.
# This is the data used to label a new subject
atlas.save(outdir,'atlas')

# Write the polydata with cluster indices saved as cell data
fname_output = os.path.join(outdir, 'clustered_whole_brain.vtp')
wma.io.write_polydata(output_polydata_s, fname_output)

# figure out which subject each fiber was from in the input to the clustering
subject_fiber_list = list()
for sidx in range(number_of_subjects):
    for fidx in range(number_of_fibers_per_subject):
        subject_fiber_list.append(sidx)
subject_fiber_list = numpy.array(subject_fiber_list)

def output_and_quality_control_cluster_atlas(atlas, output_polydata_s, subject_fiber_list, input_polydatas, outdir):

    # output summary file to save information about all subjects
    subjects_qc_fname = os.path.join(outdir, 'input_subjects.txt')
    subjects_qc_file = open(subjects_qc_fname, 'w')
    outstr = "Subject_idx\tSubject_ID\tfilename\n"
    subjects_qc_file.write(outstr)
    idx = 1
    for fname in input_polydatas:
        subject_id = os.path.splitext(os.path.basename(fname))[0]
        outstr =  str(idx) + '\t' + str(subject_id) + '\t' + str(fname) + '\n'
        subjects_qc_file.write(outstr)
        idx += 1
    subjects_qc_file.close()

    # output summary file to save information about all clusters
    clusters_qc_fname = os.path.join(outdir, 'cluster_quality_control.txt')
    clusters_qc_file = open(clusters_qc_fname, 'w')

    # Figure out how many subjects in each cluster (ideally, most subjects in most clusters)
    subjects_per_cluster = list()
    percent_subjects_per_cluster = list()
    fibers_per_cluster = list()
    mean_fiber_len_per_cluster = list()
    std_fiber_len_per_cluster = list()
    mean_fibers_per_subject_per_cluster = list()
    std_fibers_per_subject_per_cluster = list()

    # find out length of each fiber
    fiber_length = list()
    cell_idx = 0
    ptids = vtk.vtkIdList()
    inpoints = output_polydata_s.GetPoints()
    # loop over lines
    output_polydata_s.GetLines().InitTraversal()
    num_lines = output_polydata_s.GetNumberOfLines()
    for lidx in range(0, num_lines):
        output_polydata_s.GetLines().GetNextCell(ptids)
        # compute step size (assume it's fixed along line length)
        if ptids.GetNumberOfIds() >= 2:
            point0 = inpoints.GetPoint(ptids.GetId(0))
            point1 = inpoints.GetPoint(ptids.GetId(1))
            step_size = numpy.sqrt(numpy.sum(numpy.power(
                        numpy.subtract(point0, point1), 2)))
        else:
            step_size = 0.0
        fiber_length.append(ptids.GetNumberOfIds() * step_size)
    fiber_length = numpy.array(fiber_length)

    # loop over each cluster and compute quality control metrics
    cluster_indices = range(atlas.centroids.shape[0])
    for cidx in cluster_indices:
        cluster_mask = (cluster_numbers_s==cidx) 
        subjects_per_cluster.append(len(set(subject_fiber_list[cluster_mask])))
        fibers_per_subject = list()
        for sidx in range(number_of_subjects):
            fibers_per_subject.append(list(subject_fiber_list[cluster_mask]).count(sidx))
        mean_fibers_per_subject_per_cluster.append(numpy.mean(numpy.array(fibers_per_subject)))
        std_fibers_per_subject_per_cluster.append(numpy.std(numpy.array(fibers_per_subject)))
        mean_fiber_len_per_cluster.append(numpy.mean(fiber_length[cluster_mask]))
        std_fiber_len_per_cluster.append(numpy.std(fiber_length[cluster_mask]))

    percent_subjects_per_cluster = numpy.divide(numpy.array(subjects_per_cluster),float(number_of_subjects))

    # Save output quality control information
    clusters_qc_file = open(clusters_qc_fname, 'w')
    print >> clusters_qc_file, 'cluster_idx','\t', 'number_subjects','\t', 'percent_subjects','\t', 'mean_length','\t', 'std_length','\t', 'mean_fibers_per_subject','\t', 'std_fibers_per_subject'
    for cidx in cluster_indices:
        print >> clusters_qc_file, cidx + 1,'\t', subjects_per_cluster[cidx],'\t', percent_subjects_per_cluster[cidx] * 100.0,'\t', \
            mean_fiber_len_per_cluster[cidx],'\t', std_fiber_len_per_cluster[cidx],'\t', \
            mean_fibers_per_subject_per_cluster[cidx],'\t', std_fibers_per_subject_per_cluster[cidx]

    clusters_qc_file.close()

    if HAVE_PLT:
        plt.figure()
        plt.hist(subjects_per_cluster, number_of_subjects)
        plt.title('Histogram of Subjects per Cluster')
        plt.xlabel('subjects per cluster')
        plt.ylabel('number of clusters')
        plt.savefig( os.path.join(outdir, 'subjects_per_cluster_hist.pdf'))
        plt.close()
        
    # Save the entire combined atlas as individual clusters for visualization
    # and labeling/naming of structures. This will include all of the data
    # that was clustered to make the atlas.

    # Figure out file name and mean color for each cluster, and write the individual polydatas
    fnames = list()
    cluster_colors = list()
    cluster_sizes = list()
    cluster_fnames = list()
    for c in cluster_indices:
        mask = cluster_numbers_s == c
        cluster_size = numpy.sum(mask)
        cluster_sizes.append(cluster_size)
        # color by subject so in theory we can see which one it came from
        # but this is cell data and may not be correctly shown in Slicer.
        #colors = subject_fiber_list
        pd_c = wma.filter.mask(output_polydata_s, mask,verbose=False)
        # The clusters are stored starting with 1, not 0, for user friendliness.
        fname_c = 'cluster_{0:05d}.vtp'.format(c+1)
        # save the filename for writing into the MRML file
        fnames.append(fname_c)
        # prepend the output directory
        fname_c = os.path.join(outdir, fname_c)
        #print fname_c
        wma.io.write_polydata(pd_c, fname_c)
        cluster_fnames.append(fname_c)
        color_c = color[mask,:]
        cluster_colors.append(numpy.mean(color_c,0))

    # Notify user if some clusters empty
    print "<wm_cluster_atlas.py> Checking for empty clusters (should not happen in atlas clustering)."
    for sz, fname in zip(cluster_sizes,cluster_fnames):
        if sz == 0:
            print sz, ":", fname

    cluster_sizes = numpy.array(cluster_sizes)
    print "<wm_cluster_atlas.py> Mean number of fibers per cluster:", numpy.mean(cluster_sizes), "Range:", numpy.min(cluster_sizes), "..", numpy.max(cluster_sizes)

    # Estimate subsampling ratio to display approximately show_fibers total fibers in 3D Slicer
    number_fibers = len(cluster_numbers_s)
    if number_fibers < show_fibers:
        ratio = 1.0
    else:
        ratio = show_fibers / number_fibers
    #print "<wm_cluster_atlas.py> Total fibers:", number_fibers, "Fibers to show by default:", show_fibers
    print "<wm_cluster_atlas.py> Subsampling ratio for display of", show_fibers, "total fibers estimated as:", ratio

    # Write the MRML file into the directory where the polydatas were already stored
    fname = os.path.join(outdir, 'clustered_tracts.mrml')
    wma.mrml.write(fnames, numpy.around(numpy.array(cluster_colors), decimals=3), fname, ratio=ratio)

    # Also write one with 100% of fibers displayed
    fname = os.path.join(outdir, 'clustered_tracts_display_100_percent.mrml')
    wma.mrml.write(fnames, numpy.around(numpy.array(cluster_colors), decimals=3), fname, ratio=1.0)
    
    # View the whole thing in png format for quality control
    print '<wm_cluster_atlas.py> Rendering and saving images of cluster atlas.'
    ren = wma.render.render(output_polydata_s, 1000, data_mode='Cell', data_name='EmbeddingColor', verbose=False)
    ren.save_views(outdir, verbose=False)
    del ren



# Save some quality control metrics and save the atlas as individual polydata. This is used to 
# set up a mrml hierarchy file and to visualize the output. This data is not used to label
# a new subject.
output_and_quality_control_cluster_atlas(atlas, output_polydata_s, subject_fiber_list, input_polydatas, outdir)

print "==========================\n"
print '<wm_cluster_atlas.py> Done clustering atlas. See output in directory:\n ', outdir, '\n'



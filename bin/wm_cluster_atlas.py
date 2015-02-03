#!/usr/bin/env python
import numpy
import argparse
import os
import multiprocessing
import vtk

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
    description="Runs clustering of tractography for multiple subjects to create an atlas.",
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
    help='Number of fibers to analyze from each subject.')
parser.add_argument(
    '-l', action="store", dest="fiberLength", type=int,
    help='Minimum length (in mm) of fibers to analyze. 25mm is default.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')
parser.add_argument(
    '-verbose', action='store_true', dest="flag_verbose",
    help='Verbose. Run with -verbose for more text output.')
parser.add_argument(
    '-k', action="store", dest="numberOfClusters", type=int,
    help='Number of clusters to find. Default is 250. Useful range is from 200 to 600+.')
parser.add_argument(
    '-thresh', action="store", dest="distanceThreshold", type=float,
    help='Threshold (in mm) below which fiber points are considered in the same position. Default is 2mm for cross-subject atlas clustering.')
parser.add_argument(
    '-nystrom_sample', action="store", dest="sizeOfNystromSample", type=int,
    help='Number of fibers to use in the Nystrom sample. 2000 is default. Must be >1500. Increase for larger datasets. Reduce to limit memory use.')
parser.add_argument(
    '-sigma', action="store", dest="sigma", type=float,
    help='Sigma for kernel. Controls distance over which fibers are considered similar. 60mm is default. Reduce for stricter clustering with ample data, or for physically smaller problems like a subset of the brain.')
parser.add_argument(
    '-mrml_fibers', action="store", dest="showNFibersInSlicer", type=float,
    help='Approximate upper limit on number of fibers to show when MRML scene of clusters is loaded into slicer')
parser.add_argument(
    '-remove_outliers', action='store_true', dest="flag_remove_outliers",
    help='Define outlier clusters (default, ones with less than 90 percent of subjects present). These will be segmented separately.')
parser.add_argument(
    '-subject_percent', action="store", dest="subjectPercent", type=float,
    help='Threshold for defining good vs outlier clusters. Default is 90 percent of subjects must be present to keep a cluster. Lower this to 70 or 80 if data are more variable or the number of clusters is very high.')
parser.add_argument(
    '-bilateral_off', action='store_true', dest="flag_bilateral_off",
    help='Turn off bilateral clustering. In general, anatomy is better and more stably represented with bilateral clusters. They can be split at the midline later if needed')

args = parser.parse_args()


if not os.path.isdir(args.inputDirectory):
    print "<wm_cluster_atlas.py> Error: Input directory", args.inputDirectory, "does not exist or is not a directory."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "<wm_cluster_atlas.py> Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

print "<wm_cluster_atlas.py> Starting computation."
print ""
print "=====input directory ======\n", args.inputDirectory
print "=====output directory =====\n", args.outputDirectory
print "=========================="
print ""


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
    show_fibers = 5000.0
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

    
# default clustering parameters that probably don't need to be changed
# from TMI 2007 paper
use_nystrom=True
distance_method = 'Mean'
use_normalized_cuts = True
number_of_eigenvectors = 10


# Print input parameter information to a file as well as to the terminal
print args
f = open(os.path.join(outdir, 'parameters.txt'), 'w+')
print >> f, args
f.close()

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

print "<wm_cluster_atlas.py> Input number of vtk/vtp files: ", number_of_subjects

# read in data
input_pds = list()
for fname in input_polydatas:
    print fname
    # read data
    print "<wm_cluster_atlas.py> Reading input file:", fname
    pd = wma.io.read_polydata(fname)
    # preprocessing step: minimum length
    print "<wm_cluster_atlas.py> Preprocessing by length:", fiber_length, "mm."
    pd2 = wma.filter.preprocess(pd, fiber_length)
    # preprocessing step: fibers to analyze
    if number_of_fibers_per_subject is not None:
        print "<wm_cluster_atlas.py> Downsampling to ", number_of_fibers_per_subject, "fibers."
        pd3 = wma.filter.downsample(pd2, number_of_fibers_per_subject)
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
    print "<wm_cluster_atlas.py>Error Nystrom sample size is larger than number of fibers available."
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

# Save the output in our atlas format for automatic labeling of full brain datasets.
# This is the data used to label a new subject
atlas.save(outdir,'atlas')

# Write the polydata with cluster indices saved as cell data
#fname_output = os.path.join(outdir, 'clustered_whole_brain.vtp')
#wma.io.write_polydata(output_polydata_s, fname_output)

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
        print >> clusters_qc_file, cidx,'\t', subjects_per_cluster[cidx],'\t', percent_subjects_per_cluster[cidx],'\t', \
            mean_fiber_len_per_cluster[cidx],'\t', std_fiber_len_per_cluster[cidx],'\t', \
            mean_fibers_per_subject_per_cluster[cidx],'\t', std_fibers_per_subject_per_cluster[cidx]

    clusters_qc_file.close()

    if HAVE_PLT:
        plt.figure()
        plt.hist(subjects_per_cluster, number_of_subjects)
        plt.savefig( os.path.join(outdir, 'subjects_per_cluster_hist.pdf'))
        plt.close()
        
    # Save the entire combined atlas as individual clusters for visualization
    # and labeling/naming of structures. This will include all of the data
    # that was clustered to make the atlas.

    # Figure out file name and mean color for each cluster, and write the individual polydatas
    fnames = list()
    cluster_colors = list()
    for c in cluster_indices:
        mask = cluster_numbers_s == c
        # color by subject so in theory we can see which one it came from
        # but this is cell data and may not be correctly shown in Slicer.
        #colors = subject_fiber_list
        pd_c = wma.filter.mask(output_polydata_s, mask)
        # The clusters are stored starting with 1, not 0, for user friendliness.
        fname_c = 'cluster_{0:05d}.vtp'.format(c+1)
        # save the filename for writing into the MRML file
        fnames.append(fname_c)
        # prepend the output directory
        fname_c = os.path.join(outdir, fname_c)
        print fname_c
        wma.io.write_polydata(pd_c, fname_c)
        color_c = color[mask,:]
        cluster_colors.append(numpy.mean(color_c,0))
    
    # Estimate subsampling ratio to display approximately show_fibers total fibers in 3D Slicer
    number_fibers = len(cluster_numbers_s)
    if number_fibers < show_fibers:
        ratio = 1.0
    else:
        ratio = show_fibers / number_fibers
    #print "<wm_cluster_atlas.py> Total fibers:", number_fibers, "Fibers to show by default:", show_fibers
    #print "<wm_cluster_atlas.py> Subsampling ratio estimated as:", ratio

    # Write the MRML file into the directory where the polydatas were already stored
    fname = os.path.join(outdir, 'clustered_tracts.mrml')
    wma.mrml.write(fnames, numpy.around(numpy.array(cluster_colors), decimals=3), fname, ratio=ratio)

    # View the whole thing in png format for quality control
    print 'Rendering and saving images of cluster atlas.'
    ren = wma.render.render(output_polydata_s, 1000, data_mode='Cell', data_name='EmbeddingColor')
    ren.save_views(outdir)
    del ren



# Save some quality control metrics and save the atlas as individual polydata. This is used to 
# set up a mrml hierarchy file and to visualize the output. This data is not used to label
# a new subject.
output_and_quality_control_cluster_atlas(atlas, output_polydata_s, subject_fiber_list, input_polydatas, outdir)

print 'Done clustering atlas.'


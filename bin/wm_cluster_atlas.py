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

#TESTING
bilateral = True
#TESTING

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
    help='Number of clusters to find. Default is 200. Useful range is from 200 to 600+.')
parser.add_argument(
    '-thresh', action="store", dest="distanceThreshold", type=float,
    help='Threshold (in mm) below which fiber points are considered in the same position. 0mm is default. Set to 2mm for cross-subject clustering.')
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
    help='Define outlier clusters (default, ones with less than 90% of subjects present). These will be segmented separately.')
parser.add_argument(
    '-subject_percent', action="store", dest="subjectPercent", type=float,
    help='Threshold for defining good vs outlier clusters. Default is 90% of subjects must be present to keep a cluster. Lower this to 70 or 80 if data are more variable or the number of clusters is very high.')

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
    number_of_clusters = 200
print "Number of clusters to find: ", number_of_clusters

if args.distanceThreshold is not None:
    threshold = args.distanceThreshold
else:
    threshold = 0.0
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

if args.flag_remove_outliers:
    if args.subjectPercent is not None:
        fraction_to_keep_cluster = args.subjectPercent / 100.0
    else:
        fraction_to_keep_cluster = 0.9
    print "Separating outlier clusters with fewer than: ", fraction_to_keep_cluster * 100.0, "percent of subjects."
else:
    print "Not removing outlier clusters."


    
# default clustering parameters that probably don't need to be changed
# from TMI 2007 paper
use_nystrom=True
distance_method = 'Mean'
use_normalized_cuts = True
number_of_eigenvectors = 10


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

# Save the output in our atlas format for automatic labeling of full brain datasets
atlas.save(outdir,'atlas')

# Write the polydata with cluster indices saved as cell data
#fname_output = os.path.join(outdir, 'clustered_whole_brain.vtp')
#wma.io.write_polydata(output_polydata_s, fname_output)


# Save some quality control metrics.
# figure out which subject each fiber was from 
subject_fiber_list = list()
for sidx in range(number_of_subjects):
    for fidx in range(number_of_fibers_per_subject):
        subject_fiber_list.append(sidx)
subject_fiber_list = numpy.array(subject_fiber_list)
# Figure out how many subjects in each cluster (ideally, most subjects in most clusters)
subjects_per_cluster = list()
for cidx in range(atlas.centroids.shape[0]):
    cluster_mask = (cluster_numbers_s==cidx) 
    subjects_per_cluster.append(len(set(subject_fiber_list[cluster_mask])))
if HAVE_PLT:
    plt.figure()
    plt.hist(subjects_per_cluster, number_of_subjects)
    plt.savefig( os.path.join(outdir, 'subjects_per_cluster_hist.pdf'))
    plt.close()


# Figure out what clusters we have, and keep the consistent ones if requested.
number_of_clusters = numpy.max(cluster_numbers_s)
        
if args.flag_remove_outliers:
    subjects_to_keep_cluster = numpy.round(fraction_to_keep_cluster * number_of_subjects)
    print "Separating outlier clusters with fewer than: ", subjects_to_keep_cluster, "/", number_of_subjects, "subjects."   
    # experiment
    cluster_indices = numpy.where(numpy.array(subjects_per_cluster) >= subjects_to_keep_cluster)[0]
    outlier_indices = numpy.where(numpy.array(subjects_per_cluster) < subjects_to_keep_cluster)[0]
    print "Keeping", len(cluster_indices), "/", number_of_clusters, "clusters."
else:
    first_cluster = numpy.min(cluster_numbers_s)
    print "Cluster indices range from:", first_cluster, "to", number_of_clusters
    cluster_indices = range(first_cluster,number_of_clusters)

# Also save the entire combined atlas as individual clusters for visualization
# and labeling/naming of structures. This will be the subset of the data
# that was analyzed to make the atlas.

# Figure out file name and mean color for each cluster, and write the individual polydatas
fnames = list()
cluster_colors = list()
for c in cluster_indices:
    mask = cluster_numbers_s == c
    pd_c = wma.filter.mask(output_polydata_s, mask)
    fname_c = 'cluster_{0:05d}.vtp'.format(c)
    # save the filename for writing into the MRML file
    fnames.append(fname_c)
    # prepend the output directory
    fname_c = os.path.join(outdir, fname_c)
    print fname_c
    wma.io.write_polydata(pd_c, fname_c)
    color_c = color[mask,:]
    cluster_colors.append(numpy.mean(color_c,0))
    
# Estimate subsampling ratio to display approx. show_fibers total fibers
number_fibers = len(cluster_numbers_s)
if number_fibers < show_fibers:
    ratio = 1.0
else:
    ratio = show_fibers / number_fibers
print "<wm_cluster_atlas.py> Total fibers:", number_fibers, "Fibers to show by default:", show_fibers
print "<wm_cluster_atlas.py> Subsampling ratio estimated as:", ratio

# Write the MRML file into the directory where the polydatas were already stored
fname = os.path.join(outdir, 'clustered_tracts.mrml')
wma.mrml.write(fnames, numpy.around(numpy.array(cluster_colors), decimals=3), fname, ratio=ratio)

# Also save the outlier clusters for quality control
if args.flag_remove_outliers:
    outdir_outlier = os.path.join(outdir, 'outliers')
    if not os.path.exists(outdir_outlier):
        os.makedirs(outdir_outlier)
    # Figure out file name and mean color for each cluster, and write the individual polydatas
    fnames = list()
    cluster_colors = list()
    for c in outlier_indices:
        mask = cluster_numbers_s == c
        pd_c = wma.filter.mask(output_polydata_s, mask)
        fname_c = 'cluster_{0:05d}.vtp'.format(c)
        # save the filename for writing into the MRML file
        fnames.append(fname_c)
        # prepend the output directory
        fname_c = os.path.join(outdir_outlier, fname_c)
        print fname_c
        wma.io.write_polydata(pd_c, fname_c)
        color_c = color[mask,:]
        cluster_colors.append(numpy.mean(color_c,0))
        
    # Write the MRML file into the directory where the polydatas were already stored
    fname = os.path.join(outdir_outlier, 'outlier_tracts.mrml')
    wma.mrml.write(fnames, numpy.around(numpy.array(cluster_colors), decimals=3), fname, ratio=ratio)

    # Save the outlier information into the atlas...
    print "LAUREN SAVE OUTLIER and qc number of subjects INFORMATION IN ATLAS AND MOVE THAT CODE TO cluster.py"
    
# View the whole thing in png format for quality control
print '<wm_cluster_atlas.py> Rendering and saving image'
ren = wma.render.render(output_polydata_s, 1000)
ren.save_views(outdir)
del ren

print 'Done clustering atlas.'

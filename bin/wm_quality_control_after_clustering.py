#!/usr/bin/env python
import argparse
import os
import numpy

try:
    import whitematteranalysis as wma
except:
    print "<wm_quality_control_after_clustering> Error importing white matter analysis package\n"
    raise

HAVE_PLT = 1
try:
    import matplotlib
    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except:
    print "<wm_cluster.py> Error importing matplotlib.pyplot package, can't plot quality control data.\n"
    HAVE_PLT = 0

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Quality control of fiber clustering results across multiple subjects.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputDirectory',
    help='Directory of fiber clustering results obtained by <wm_cluster_from_altas.py> of multiple subjects. Make sure only the fiber clustering results are stored in this folder, making one subdirectory corresponding to one subject.')
parser.add_argument(
    'outputDirectory',
    help='Quality control information will be stored in the output directory, which will be created if it does not exist.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')

args = parser.parse_args()

if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

output_dir = args.outputDirectory
if not os.path.exists(output_dir):
    print "<register> Output directory", output_dir, "does not exist, creating it."
    os.makedirs(output_dir)

subject_list = os.listdir(args.inputDirectory)
print "<wm_quality_control_after_clustering> found", len(subject_list), "subjects."

# Check if all subjects have the same number of clusters.
# This can also help find the sub folders that are not the fiber clustering results.
num_of_subjects = len(subject_list)
flag_same_num_clusters = 1
for sidx in range(0, num_of_subjects):
    sub = subject_list[sidx]
    cluster_polydatas = wma.io.list_vtk_files(os.path.join(args.inputDirectory, sub))
    print "  ", sub, "has", len(cluster_polydatas), "clusters."
    if sidx == 0:
        num_of_clusters = len(cluster_polydatas)
    else:
        if num_of_clusters != len(cluster_polydatas):
            flag_same_num_clusters = 0

if flag_same_num_clusters == 0:
    print "Error: All subjects should have the number of clusters. (Hint: there may be folders that are not the fiber clustering results.)"
    exit()

# Read number of fibers per cluster per subject
print "<wm_quality_control_after_clustering.py> calculate the number of fibers per cluster per subject."
num_fibers_per_subject = numpy.zeros([num_of_subjects, num_of_clusters])
for sidx in range(0, num_of_subjects):
    sub = subject_list[sidx]
    print "   loading", sub
    cluster_polydatas = wma.io.list_vtk_files(os.path.join(args.inputDirectory, sub))
    for cidx in range(0, num_of_clusters):
        fname = cluster_polydatas[cidx]
        pd = wma.io.read_polydata(fname)
        num_fibers_per_subject[sidx, cidx] = pd.GetNumberOfLines()

subjects_per_cluster = numpy.sum(num_fibers_per_subject > 0, axis=0)
#clusters_per_subject = numpy.sum(num_fibers_per_subject > 0, axis=1)

percent_subjects_per_cluster = numpy.divide(subjects_per_cluster, float(num_of_subjects))

clusters_qc_fname = os.path.join(output_dir, 'cluster_quality_control.txt')
print "<wm_quality_control_after_clustering.py> Saving cluster quality control information file."
clusters_qc_file = open(clusters_qc_fname, 'w')
print >> clusters_qc_file, 'cluster_idx','\t', 'number_subjects', 'percent_subjects'
for cidx in range(0, num_of_clusters):
    print >> clusters_qc_file, cidx + 1,'\t', subjects_per_cluster[cidx],'\t', percent_subjects_per_cluster[cidx] * 100.0
clusters_qc_file.close()

if HAVE_PLT:
    print "<wm_quality_control_after_clustering.py> Saving subjects per cluster histogram."
    plt.figure()
    plt.hist(subjects_per_cluster, num_of_subjects)
    plt.title('Histogram of Subjects per Cluster')
    plt.xlabel('subjects per cluster')
    plt.ylabel('number of clusters')
    plt.savefig(os.path.join(output_dir, 'subjects_per_cluster_hist.pdf'))
    plt.close()
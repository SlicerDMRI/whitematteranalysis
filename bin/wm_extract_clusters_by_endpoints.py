#!/usr/bin/env python
import argparse
import os
import numpy

try:
    import whitematteranalysis as wma
except:
    print "<wm_endpoint_analysis> Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Extract fiber clusters that connect to one particular region, such as one cortical parcel or one fMRI functional area. This is based on the results from <wm_measure_endpoint_overlap.py>.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputDirectory',
    help='Contains endpoint overlap percentage measurement as txt files, which are computed from <wm_measure_endpoint_overlap.py>.')
parser.add_argument(
    '-r', type=str,  dest="region",
    help='Label of region-of-interest as in the header of the measurement files.')
parser.add_argument(
    '-p', type=int,  dest="percentage_threshold", default='10',
    help='Percentage threshold [0, 100] that is used to decide if fiber cluster connect to the input region in individual subject. For each subject, if over p/100 of a fiber cluster\'s endpoints overlap with the region-of-interest, we consider this cluster and the input region are connected in this subject.')
parser.add_argument(
    '-s', type=int,  dest="subject_number_threshold", default='1',
    help='Subject number threshold that is used to decide if fiber cluster connect to the input region in the population. Across all subjects, if a cluster is connected to the input region in at least s subjects, we consider this cluster as one output result.')
parser.add_argument(
    '-fc_folder', type=str, dest="fiber_cluster_folder",
    help='Contains the fiber clustering result, in which each sub folder represents one subject. If specified, a mrml file that displays all output clusters will be generated. ')

args = parser.parse_args()

if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

percentage_threshold = args.percentage_threshold
if percentage_threshold < 0 or percentage_threshold > 100:
    print "Error: Percentage threshold should be in 0 to 100:"
    exit()

measurement_list = wma.tract_measurement.load_measurement_in_folder(args.inputDirectory, hierarchy = 'Column', separator = 'Tab')
if len(measurement_list) < 1:
    print "Error: No measurement files found in directory:", args.inputDirectory
    exit()

number_of_subjects = len(measurement_list)
print "\n== Number of subjects data found:", number_of_subjects

# print the header
print "\n== Region labels found in the measurement file:"
header = measurement_list[0].measurement_header
print header

region = args.region
print "\n== Input region label:", region

region_exist = False
region_index = 0
for variable in header:
    if region == variable:
        region_exist = True
        break
    region_index = region_index + 1

if not region_exist:
    print "\nError: region", region, "does not exist."
    exit()
else:
    print "Region found at position", region_index+1

subject_number_threshold = args.subject_number_threshold
if subject_number_threshold < 0 or subject_number_threshold > number_of_subjects:
    print "Error: Subject number threshold should be in 0 to number of subject ("+ str(number_of_subjects)+ ")"
    exit()

# sanity check all measured variables are the same
print "\n== Testing if all subjects have the same measurement header."
test = True
for subject_measured in measurement_list:
    all_same = all(subject_measured.measurement_header == header)
    if not all_same:
        print "ERROR: Subject does not have same measurement header:", subject_measured.case_id
        test = False
if test:
    print "Passed. All subjects have the same measurement header."

# Make subject id list for reporting of any potential outliers
subject_id_list = []
for subject_measured in measurement_list:
    subject_id_list.append(subject_measured.case_id)
subject_id_list = numpy.array(subject_id_list)

# sanity check number of clusters is the same for all subjects
print "\n== Testing if all subjects have the same number of clusters."
test = True
number_of_clusters = measurement_list[0].measurement_matrix[:,0].shape[0]
print "Group number of clusters (from first subject):", number_of_clusters
for subject_measured, id in zip(measurement_list, subject_id_list):
    nc = subject_measured.measurement_matrix[:,0].shape[0]
    if not nc == number_of_clusters:
        print nc, " : ", id
        test = False
if test:
    print "Passed. All subjects have the same number of clusters ("+ str(number_of_clusters) + ")."
else:
    print "ERROR: All subjects do not have the same number of clusters. There was an earlier error in clustering or measurement that must be fixed before analysis."

print "\n== Histogram of the overlap percentage of each cluster to the input region."
percentage_range = range(0, 100, 5)
print_title = False
for sub_id, subject_measured in zip(subject_id_list, measurement_list):
    subject_percentage_distribution = subject_measured.measurement_matrix[:, region_index]

    percent_str = sub_id
    title_str = ("{:<"+str(len(sub_id))+"}").format('Percentage:')
    for p in percentage_range:
        percent_str = percent_str + ', ' + "{:<4}".format(str(sum(subject_percentage_distribution > p)))
        title_str = title_str + ', ' + "{:<4}".format(str(p))

    if not print_title:
        print title_str
        print_title = True

    print percent_str

print "\n== Number of clusters whose percentages are over the percentage threshold."
connected_matrix = []
for sub_id, subject_measured in zip(subject_id_list, measurement_list):
    subject_percentage_distribution = subject_measured.measurement_matrix[:, region_index]
    connected_matrix.append(subject_percentage_distribution > percentage_threshold)
    print sub_id, "has", str(sum(subject_percentage_distribution > percentage_threshold)), "/", number_of_clusters, " clusters connected to the input region."

print "\n== Number of clusters that connect to the input region across all subjects."
num_subjects_per_cluster = sum(connected_matrix)
for num_idx in range(number_of_subjects+1):
    print "Clusters that are connected in", num_idx, '(subject number threshold) subjects: ', sum(num_subjects_per_cluster == num_idx)

print "\n== Result clusters that are detected by at least",subject_number_threshold,'subjects'
output_cluster_idx = numpy.where(num_subjects_per_cluster >= subject_number_threshold)[0]
output_cluster_idx = output_cluster_idx + 1
print "Total", len(output_cluster_idx), "clusters."
print output_cluster_idx

if not (args.fiber_cluster_folder is None):
    print "\n== Output mrml to the fiber cluster folder"
    sub_folder = os.listdir(args.fiber_cluster_folder)[0]
    pd_cluster_3 = wma.io.list_vtk_files(os.path.join(args.fiber_cluster_folder, sub_folder))[2]
    suffix = os.path.split(pd_cluster_3)[1][13:-4]

    cluster_polydatas = []
    for c_idx in output_cluster_idx:
        cluster_polydatas.append("cluster_"+str(c_idx).zfill(5)+suffix+".vtp")

    number_of_files = len(cluster_polydatas)

    step = int(100 * 255.0 / (number_of_files - 1))
    R = numpy.array(range(0, 100 * 255 + 1, step)) / 100.0
    G = numpy.abs(range(100 * -127, 100 * 128 + 1, step)) * 2.0 / 100.0
    B = numpy.array(range(100 * 255 + 1, 0, -step)) / 100.0

    colors = list()
    idx = 0
    for pd in cluster_polydatas:
        colors.append([R[idx], G[idx], B[idx]])
        idx += 1
    colors = numpy.array(colors)

    mrml_filename = "cluster_connecting_region_" + region + ".mrml"

    for sub_folder in os.listdir(args.fiber_cluster_folder):
        wma.mrml.write(cluster_polydatas, colors, os.path.join(args.fiber_cluster_folder+'/'+sub_folder, mrml_filename), ratio=1.0)
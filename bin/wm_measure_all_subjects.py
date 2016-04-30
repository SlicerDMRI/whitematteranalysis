#!/usr/bin/env python
import argparse
import os
from joblib import Parallel, delayed

try:
    import whitematteranalysis as wma
except:
    print "<wm_measure_all_subjects> Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Extract scalar measurement of fiber cluster.",
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
    help='Directory of output CSV files of fiber scalar measurement (computed using Slicer FiberTractMeasurements module).')
parser.add_argument(
    'Slicer',
    help='Path of 3D Slicer.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')

args = parser.parse_args()

if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

if not os.path.exists(args.Slicer):
    print "Error: 3D Slicer", args.Slicer, "does not exist."
    exit()

module_FTSM = args.Slicer + ' --launch FiberTractMeasurements '

if args.numberOfJobs is not None:
    number_of_jobs = args.numberOfJobs
else:
    number_of_jobs = 1

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

print "<wm_measure_all_subjects>. Starting scalar measurement extraction."
print ""
print "=====input directory======\n", args.inputDirectory
print "=====output directory=====\n", args.outputDirectory
print "=====3D Slicer====\n", args.Slicer
print '=====Using N jobs:', number_of_jobs, "====\n"
print "=========================="

sub_dirs = os.listdir(args.inputDirectory)

print "<wm_measure_all_subjects> found", len(sub_dirs), "sub directories."

subject_list = []
for dir in sub_dirs:
    # consider one subdirectory as one subject
    subject_list.append(dir)

print "<wm_measure_all_subjects> found", len(subject_list), "subjects."

def extract_measures(sub, subject_list, module_FTSM, args):
    sub_name = sub

    count = subject_list.index(sub)
    print " -", count+1, "/", len(subject_list), " subject id:", sub_name

    os.system(module_FTSM + \
              ' --inputtype Fibers_File_Folder --format Column_Hierarchy --separator Tab ' + \
              ' --inputdirectory ' + os.path.join(args.inputDirectory, sub) + \
              ' --outputfile ' + os.path.join(args.outputDirectory, sub_name+'.txt ') + \
              ' > ' + os.path.join(args.outputDirectory, 'tmp'+sub_name))
    os.remove(os.path.join(args.outputDirectory, 'tmp'+sub_name))


Parallel(n_jobs=number_of_jobs, verbose=1)(
                        delayed(extract_measures)(sub, subject_list, module_FTSM, args)
                        for sub in subject_list)

print "<wm_measure_all_subjects> Measurements from", len(subject_list), "subjects were extracted."

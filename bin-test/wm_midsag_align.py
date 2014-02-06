#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/ipython

import os
import glob
import matplotlib.pyplot as plt
import numpy
import time
import argparse
import multiprocessing

import vtk

try:
    import whitematteranalysis as wma
except:
    print "<wm_laterality.py> Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Optimizes symmetry of brain across midsagittal plane. Preprocess for laterality analysis.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputFileName',
    help='A file of whole-brain tractography as vtkPolyData (.vtk or .vtp).')
parser.add_argument(
    'outputDirectory',
    help='The output directory will be created if it does not exist.')
parser.add_argument(
    '-f', action="store", dest="numberOfFibers", type=int,
    help='Number of fibers to analyze from each dataset. 2000 is reasonable.')
parser.add_argument(
    '-l', action="store", dest="fiberLength", type=int,
    help='Minimum length (in mm) of fibers to analyze. 75mm is default.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')
parser.add_argument(
    '-verbose', action='store_true', dest="flag_verbose",
    help='Verbose. Run with -verbose to store images of intermediate and final polydatas.')

args = parser.parse_args()

if not os.path.isfile(args.inputFileName):
    print "Error: Input file name", args.inputFileName, "does not exist."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)


print "wm_midsag_align. Starting symmetric alignment computation."
print ""
print "=====input file name======\n", args.inputFileName
print "=====output directory =====\n", args.outputDirectory
print "=========================="

if args.numberOfFibers is not None:
    print "fibers to analyze per subject: ", args.numberOfFibers
else:
    print "fibers to analyze per subject: ALL"

number_of_fibers = args.numberOfFibers

if args.fiberLength is not None:
    print "minimum length of fibers to analyze (in mm): ", args.fiberLength
    fiber_length = args.fiberLength
else:
    fiber_length = 75

if args.flag_verbose:
    print "Verbose display and intermediate image saving ON."
else:
    print "Verbose display and intermediate image saving OFF."
verbose = args.flag_verbose

print 'CPUs detected:', multiprocessing.cpu_count()
if args.numberOfJobs is not None:
    parallel_jobs = args.numberOfJobs
else:
    parallel_jobs = multiprocessing.cpu_count()
print 'Using N jobs:', parallel_jobs

input_poly_data = args.inputFileName

print input_poly_data


## run the registration ONCE and output result to disk
points_per_fiber = 5
number_of_fibers_per_step = [100, 200, 200, 250]
# small sigmas only, this is a minor adjustment 
sigma_per_step = [10, 10, 5, 5]
# cobyla: number of objective function calculations per iteration
maxfun = 30


# registration
wma.registration_functions.run_midsag_align(input_poly_data, outdir,
                                            number_of_fibers=number_of_fibers,
                                            points_per_fiber=points_per_fiber,
                                            parallel_jobs=parallel_jobs,
                                            number_of_fibers_per_step=number_of_fibers_per_step,
                                            sigma_per_step=sigma_per_step,
                                            maxfun=maxfun,
                                            fiber_length=fiber_length,
                                            verbose=verbose)






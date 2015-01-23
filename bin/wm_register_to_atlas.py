#!/usr/bin/env python
import os
import glob
import matplotlib.pyplot as plt
import numpy
import time
import multiprocessing
import argparse

import vtk

try:
    import whitematteranalysis as wma
except:
    print "<wm_register.py> Error importing white matter analysis package\n"
    raise

# defaults for two subject registration
# (may be added as parameters later)
fiber_sample_fractions = [.4, .5, .6, .8]
#sigma_per_scale = [30, 10, 10, 5]
sigma_per_scale = [10, 10, 10, 5]
#steps_per_scale=[10, 3, 2, 2]
steps_per_scale=[5, 2, 2, 2]
fibers_rendered = 100

# figure out how many cobyla iterations are needed
# this is reasonable for two subjects
#maxfun_per_scale = [20, 40, 60, 80]


#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Registers all whole-brain vtk tractography files in one directory to another vtk tractography file (an atlas).",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"Unbiased Groupwise Registration of White Matter Tractography. LJ O'Donnell,  WM Wells III, Golby AJ, CF Westin. Med Image Comput Comput Assist Interv. 2012;15(Pt 3):123-30.\"",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='A directory of whole-brain tractography as vtkPolyData (.vtk or .vtp), for registration to the atlas.')
parser.add_argument(
    'inputAtlas',
    help='An atlas, one file containing whole-brain tractography as vtkPolyData (.vtk or .vtp).')
parser.add_argument(
    'outputDirectory',
    help='The output directory will be created if it does not exist.')
parser.add_argument(
    '-f', action="store", dest="numberOfFibers", type=int,
    help='Number of fibers to analyze from each dataset. 300-600 or more is reasonable. Depends on total number of datasets and desired run time/memory use.')
parser.add_argument(
    '-l', action="store", dest="fiberLength", type=int,
    help='Minimum length (in mm) of fibers to analyze. 60mm is default.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')
parser.add_argument(
    '-verbose', action='store_true', dest="flag_verbose",
    help='Verbose. Run with -verbose to store images of intermediate and final polydatas.')
parser.add_argument(
    '-pf', action="store", dest="pointsPerFiber", type=int,
    help='Number of points for fiber representation during registration. 5 is reasonable, or more.')

 
args = parser.parse_args()

print "\n\n<register> =========GROUP REGISTRATION============"
print "<register> Performing unbiased group registration."
print "<register> Input  directory: ", args.inputDirectory
print "<register> Output directory: ", args.outputDirectory
print "\n<register> ============PARAMETERS================="

if not os.path.isdir(args.inputDirectory):
    print "<register> Error: Input directory", args.inputDirectory, "does not exist."
    exit()

if not os.path.isfile(args.inputAtlas):
    print "<register> Error: Input atlas", args.inputAtlas, "does not exist."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "<register> Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

if args.numberOfFibers is not None:
    number_of_fibers = args.numberOfFibers
else:
    number_of_fibers = 300
print "<register> Number of fibers to analyze per subject: ", number_of_fibers

if args.fiberLength is not None:
    fiber_length = args.fiberLength
else:
    fiber_length = 75
print "<register> Minimum length of fibers to analyze (in mm): ", fiber_length
    
if args.numberOfJobs is not None:
    parallel_jobs = args.numberOfJobs
else:
    parallel_jobs = multiprocessing.cpu_count()
print "<register> CPUs detected:", multiprocessing.cpu_count(), ". Number of jobs to use:", parallel_jobs

if args.flag_verbose:
    print "<register> Verbose display and intermediate image saving ON."
else:
    print "<register> Verbose display and intermediate image saving OFF."
verbose = args.flag_verbose

if args.pointsPerFiber is not None:
    points_per_fiber = args.pointsPerFiber
else:
    points_per_fiber = 5
print "<register> Number of points for fiber representation: ", points_per_fiber


print "\n<register> Starting registration...\n"



## run the registration ONCE and output result to disk
elapsed = wma.registration_functions.run_atlas_registration(args.inputDirectory,
                                                            args.inputAtlas,
                                                            outdir,
                                                            number_of_fibers=number_of_fibers,
                                                            points_per_fiber=points_per_fiber,
                                                            parallel_jobs=parallel_jobs,
                                                            fiber_sample_fractions=fiber_sample_fractions,
                                                            sigma_per_scale=sigma_per_scale,
                                                            verbose=verbose,
                                                            fiber_length=fiber_length,
                                                            fibers_rendered=fibers_rendered,
                                                            steps_per_scale=steps_per_scale)

print "TIME:", elapsed

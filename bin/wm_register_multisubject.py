#!/usr/bin/env python
#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/ipython

# Run registration on the test dataset.

import argparse
import os
import multiprocessing

try:
    import whitematteranalysis as wma
except:
    print "<wm_register.py> Error importing white matter analysis package\n"
    raise

# defaults for multiple subject registration
# (may be added as parameters later)
fiber_sample_fractions = [.10, .20, .30, .40]
#sigma_per_scale = [30, 10, 10, 5]
# Note Dec 15 2014. Scale 3 is not decreasing the objective function much
# on two-tensor tractography datasets with 10-12 subjects. Perhaps it 
# has already converged at sigma of 10 and more computation is a waste of time. 
# So try decreasing sigma to 7.5 for the third scale space step.
# This seems to do a better, finer registration.
sigma_per_scale = [30, 10, 7.5, 5]
# on a 27-subject dataset, 12 or 14 processors, this takes nearly 24 hours.
# the implementation is not fast, but this is a long time to wait.
#steps_per_scale=[10, 3, 2, 2]
# the big gain in objective from the 1st scale is early,
# and actually after the fine (3rd scale) registration
# the result looks good enough. So, save output after third step,
# and also reduce the time spent in first and last scales.
# note: this is still more computation than in the publication.
# Trying to decide how much computation is needed for desired result quality.
steps_per_scale=[5, 3, 2, 1]

fibers_rendered = 100

# figure out how many cobyla iterations are needed
# this is reasonable for two subjects
#maxfun_per_scale = [20, 40, 60, 80]


#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Runs multisubject unbiased group registration of tractography.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"Unbiased Groupwise Registration of White Matter Tractography. LJ O'Donnell,  WM Wells III, Golby AJ, CF Westin. Med Image Comput Comput Assist Interv. 2012;15(Pt 3):123-30.\"",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='A directory of whole-brain tractography as vtkPolyData (.vtk or .vtp).')
parser.add_argument(
    'outputDirectory',
    help='The output directory will be created if it does not exist.')
parser.add_argument(
    '-f', action="store", dest="numberOfFibers", type=int,
    help='Number of fibers to analyze from each dataset. 300-2000 or more is reasonable. Depends on total number of datasets and desired run time/memory use. Default setting is for a fast test run: 300 fibers per subject. A good setting for a paper is 1000-2000 per subject, if greater than 10 subjects.')
parser.add_argument(
    '-l', action="store", dest="fiberLength", type=int,
    help='Minimum length (in mm) of fibers to analyze. 60mm is default (good for DTI single-tensor tractography which is shorter in general). Use a higher value such as 80 or 100 for two-tensor or other advanced tractography. This parameter removes short, noisy fibers and focuses on larger structures that can be registered well.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')
parser.add_argument(
    '-verbose', action='store_true', dest="flag_verbose",
    help='Verbose. Run with -verbose to store more files and images of intermediate and final polydatas.')
parser.add_argument(
    '-pf', action="store", dest="pointsPerFiber", type=int,
    help='Number of points for fiber representation during registration. The default of 5 is reasonable.')
parser.add_argument(
    '-norender', action='store_true', dest="flag_norender",
    help='No Render. Prevents rendering of images that would require an X connection.')
parser.add_argument(
    '-distance_method', action="store", dest="distance_method", type=str,
    help='Distance method to use. Options are Hausdorff, Mean, and MeanSquared. Changing this from the default of Hausdorff is not recommended.')
parser.add_argument(
    '-midsag_symmetric', action="store_true", dest="flag_midsag_symmetric",
    help='Register all subjects including reflected copies of input subjects, for a symmetric registration.')

 
 
args = parser.parse_args()

print "\n\n<register> =========GROUP REGISTRATION============"
print "<register> Performing unbiased group registration."
print "<register> Input  directory: ", args.inputDirectory
print "<register> Output directory: ", args.outputDirectory
print "\n<register> ============PARAMETERS================="

if not os.path.isdir(args.inputDirectory):
    print "<register> Error: Input directory", args.inputDirectory, "does not exist."
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

if args.flag_norender:
    print "<register> No rendering (for compute servers without X connection)."
else:
    print "<register> Rendering. For intermediate image saving to check progress."
no_render = args.flag_norender

if args.distance_method is not None:
    distance_method = args.distance_method
else:
    distance_method = 'Hausdorff'
    
print "<register> Distance method: ", distance_method

if args.flag_midsag_symmetric:
    print "<register> Midsag_Symmetric registration ON."
else:
    print "<register> Midsag_Symmetric registration OFF."
midsag_symmetric = args.flag_midsag_symmetric

print "\n<register> Starting registration...\n"



## run the registration ONCE and output result to disk
register, elapsed = wma.registration_functions.run_multisubject_registration(args.inputDirectory,
                                                                             outdir,
                                                                             number_of_fibers=number_of_fibers,
                                                                             points_per_fiber=points_per_fiber,
                                                                             parallel_jobs=parallel_jobs,
                                                                             fiber_sample_fractions=fiber_sample_fractions,
                                                                             sigma_per_scale=sigma_per_scale,
                                                                             verbose=verbose,
                                                                             fiber_length=fiber_length,
                                                                             fibers_rendered=fibers_rendered,
                                                                             steps_per_scale=steps_per_scale,
                                                                             no_render=no_render,
                                                                             distance_method=distance_method,
                                                                             midsag_symmetric=midsag_symmetric)

print "TIME:", elapsed

#!/usr/bin/env python
#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/ipython

import argparse
import os
import multiprocessing

import numpy
import vtk

try:
    import whitematteranalysis as wma
except:
    print "<wm_register.py> Error importing white matter analysis package\n"
    raise

try:
    import scipy.optimize
except ImportError:
    print ""
    print "<wm_register_multisubject.py> ERROR: Failed to import scipy.optimize, cannot run registration."
    print "Please install scipy."
    print ""
    exit()


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
    '-mode', action="store", dest="mode", type=str, default="affine",
    help='The mode can be affine or nonlinear. Affine is the default. It should be run first before nonlinear.')
parser.add_argument(
    '-f', action="store", dest="numberOfFibers", type=int, default=20000,
    help='Total number of fibers to analyze from each dataset. During registration, at each iteration fibers are randomly sampled from within this data. 20000 is the default number of total fibers.')
parser.add_argument(
    '-l', action="store", dest="fiberLength", type=int, default=80,
    help='Minimum length (in mm) of fibers to analyze. 60mm is reasonable for DTI single-tensor tractography which is shorter in general. Use a higher value such as 80 or 100 for two-tensor or other advanced tractography. This parameter removes short, noisy fibers and focuses on larger structures that can be registered well. For neonate data, a value of 40mm is suggested. The default is 80mm.')
parser.add_argument(
    '-lmax', action="store", dest="fiberLengthMax", type=int, default=150,
    help='Maximum length (in mm) of fibers to analyze. This parameter can be used to remove extremely long fibers that may have traversed several structures. For example, a value of 150 will avoid sampling the tail end of the fiber length distribution. The default is 150 mm.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')
parser.add_argument(
    '-verbose', action='store_true', dest="flag_verbose",
    help='Verbose. Run with -verbose to store more files and images of intermediate and final polydatas.')
#parser.add_argument(
#    '-pf', action="store", dest="pointsPerFiber", type=int, default=15,
#    help='Number of points for fiber representation during registration. The default of 15 is reasonable.')
parser.add_argument(
    '-norender', action='store_true', dest="flag_norender",
    help='No Render. Prevents rendering of images that would require an X connection.')
parser.add_argument(
    '-midsag_symmetric', action="store_true", dest="flag_midsag_symmetric",
    help='Register all subjects including reflected copies of input subjects, for a symmetric registration.')
parser.add_argument(
    '-advanced_only_random_seed', action='store', dest="randomSeed", type=int,
    help='(Advanced parameter for testing only.) Set random seed for reproducible clustering in software tests.')
 
 
args = parser.parse_args()

print "\n\n<register> =========GROUP REGISTRATION============"
print "<register> Performing unbiased group registration."
print "<register> Input  directory: ", args.inputDirectory
print "<register> Output directory: ", args.outputDirectory
print "\n<register> ============PARAMETERS================="

mode = args.mode
print "<register> Registration mode:", mode

if not os.path.isdir(args.inputDirectory):
    print "<register> Error: Input directory", args.inputDirectory, "does not exist."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "<register> Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

number_of_fibers = args.numberOfFibers
print "<register> Number of fibers to analyze per subject: ", number_of_fibers

fiber_length = args.fiberLength
print "<register> Minimum length of fibers to analyze (in mm): ", fiber_length

fiber_length_max = args.fiberLengthMax
print "<register> Maximum  length of fibers to analyze (in mm): ", fiber_length_max

if args.numberOfJobs is not None:
    parallel_jobs = args.numberOfJobs
    print "<register> Number of jobs to use:", parallel_jobs
else:
    parallel_jobs = multiprocessing.cpu_count()
    print "<register> CPUs detected:", multiprocessing.cpu_count(), ". Number of jobs to use:", parallel_jobs

if args.flag_verbose:
    print "<register> Verbose display and intermediate image saving ON."
else:
    print "<register> Verbose display and intermediate image saving OFF."
verbose = args.flag_verbose


#points_per_fiber = args.pointsPerFiber
#print "<register> Number of points for fiber representation: ", points_per_fiber

if args.flag_norender:
    print "<register> No rendering (for compute servers without X connection)."
else:
    print "<register> Rendering. For intermediate image saving to check progress."
no_render = args.flag_norender

if args.flag_midsag_symmetric:
    print "<register> Midsag_Symmetric registration ON."
else:
    print "<register> Midsag_Symmetric registration OFF."
midsag_symmetric = args.flag_midsag_symmetric

if args.randomSeed is not None:
    print "<register> Setting random seed to: ", args.randomSeed
random_seed = args.randomSeed

# Test the input files exist
input_polydatas = wma.io.list_vtk_files(args.inputDirectory)
number_of_subjects = len(input_polydatas)
print "<register> Found ", number_of_subjects, "subjects in input directory:", args.inputDirectory
if number_of_subjects < 1:
    print "\n<register> Error: No .vtk or .vtp files were found in the input directory.\n"
    exit()

# Get input data
input_pds, subject_ids = wma.io.read_and_preprocess_polydata_directory(args.inputDirectory, fiber_length, number_of_fibers, random_seed, fiber_length_max)
        
# If we are registering for symmetry, include reflected copy of each brain
if midsag_symmetric:
    input_pds2 = []
    subject_ids2 = []
    for pd, id in zip(input_pds, subject_ids):
        # make the reflected copy
        trans = vtk.vtkTransform()
        trans.Scale(-1,1,1)
        transformer = vtk.vtkTransformPolyDataFilter()
        transformer.SetTransform(trans)
        if (vtk.vtkVersion().GetVTKMajorVersion() > 5.0):
            transformer.SetInputData(pd)
        else:
            transformer.SetInput(pd)
        transformer.Update()
        pd_reflect = transformer.GetOutput()
        # append input and its reflection to list of data to register
        input_pds2.append(pd)
        input_pds2.append(pd_reflect)
        # create subject id for reflected
        subject_ids2.append(id)
        subject_ids2.append(str(id) + '_reflect')
    input_pds = input_pds2
    subject_ids = subject_ids2

print "\n<register> Starting registration...\n"

register = wma.congeal_multisubject.MultiSubjectRegistration()
register.input_directory = args.inputDirectory
register.output_directory = args.outputDirectory

register.parallel_jobs = parallel_jobs
register.render = not no_render

# -------------
# SETTINGS
# -------------

points_per_fiber = 10

if mode == "affine":
    # default affine for adults with fewer fibers than neonate brains (lower variability)
    sigma_per_scale = [30, 10, 7.5, 5]
    iterations_per_scale=[2, 3, 3, 3]
    maxfun_per_scale = [45, 60, 75, 90]
    mean_brain_size_per_scale = [1000, 3000, 4000, 5000]
    subject_brain_size_per_scale = [250, 1500, 1750, 2000]
    initial_step_per_scale = [10, 5, 5, 5]
    final_step_per_scale = [5, 2, 2, 2]
    register.nonlinear = False
    
elif mode == "affine_neonate":
    # Try to improve previous affine registration using more subjects
    # Try smaller sigma for neonates
    sigma_per_scale = [20, 10, 7.5, 5]
    iterations_per_scale=[2, 3, 3, 3]
    maxfun_per_scale = [50, 80, 80, 200]
    mean_brain_size_per_scale = [1000, 3000, 5000, 7500]
    subject_brain_size_per_scale = [250, 1500, 2000, 2500]
    initial_step_per_scale = [10, 5, 5, 5]
    final_step_per_scale = [5, 2, 2, 2]
    register.nonlinear = False
    
elif mode == "nonlinear":
    # this is in mm space.
    initial_step_per_scale = [5, 3, 2]
    final_step_per_scale = [2, 1, 1]
    # use only very local information (small sigma)
    sigma_per_scale = [3, 2, 1]
    # how many times to repeat the process at each scale
    iterations_per_scale = [2, 2, 2]
    # this takes twice as long--still testing what is optimal for this parameter
    #iterations_per_scale = [4, 4, 4]
    # these are small samples to go (relatively) quickly: the goal is just to improve the mean brain each time
    mean_brain_size_per_scale = [2500, 2750, 3000]
    subject_brain_size_per_scale = [500, 750, 900]
    # stop computation early. no need to ever converge, just improve objective as quickly as possible
    # These settings are for a 5x5x5 grid, 125*3 = 375 parameter space.
    maxfun_per_scale = [500, 500, 750]
    # fiber representation for computation.
    points_per_fiber = 15
    register.nonlinear = True
    register.nonlinear_grid_resolution = 5
    
elif mode == "nonlinear_fine":
    # this is in mm space.
    initial_step_per_scale = [2, 2, 2]
    final_step_per_scale = [1, 1, 1]
    # use only very local information (small sigma)
    sigma_per_scale = [3, 2, 1]
    # how many times to repeat the process at each scale
    iterations_per_scale = [2, 2, 2]
    # this takes twice as long--still testing what is optimal for this parameter
    #iterations_per_scale = [4, 4, 4]
    # these are small samples to go (relatively) quickly: the goal is just to improve the mean brain each time
    mean_brain_size_per_scale = [2500, 2750, 3000]
    subject_brain_size_per_scale = [500, 750, 900]
    # stop computation early. no need to ever converge, just improve objective as quickly as possible
    # These settings are for a 6x6x6 grid, 216*3 = 648 parameter space.
    maxfun_per_scale = [700, 1300, 2000]
    # fiber representation for computation.
    points_per_fiber = 15
    register.nonlinear = True
    register.nonlinear_grid_resolution = 6
    register.initialize_nonlinear_grid()
    
elif mode == "affineTEST":
    # very quick test if software is working
    sigma_per_scale = [30, 10, 7.5]
    iterations_per_scale=[1, 1, 1]
    maxfun_per_scale = [60, 80, 100]
    mean_brain_size_per_scale = [1500, 2000, 3000]
    subject_brain_size_per_scale = [100, 500, 1000]
    initial_step_per_scale = [5, 5, 5, 5]
    final_step_per_scale = [2, 2, 2, 2]
    register.nonlinear = False
    
elif mode == "nonlinearTEST":
    # very quick test if software is working
    initial_step_per_scale = [5, 3, 1]
    final_step_per_scale = [2, 1, 0.05]
    sigma_per_scale = [3, 2, 1]
    iterations_per_scale = [1, 1, 1]
    mean_brain_size_per_scale = [1500, 2000, 3000]
    subject_brain_size_per_scale = [500, 750, 1000]
    # stop computation: this is just a quick test the software is working
    maxfun_per_scale = [10, 10, 10]
    points_per_fiber = 15
    register.nonlinear = True

else:
    print "\n<register> Error: Unknown registration mode:", mode
    exit()

# We have to add polydatas after setting nonlinear in the register object
for (pd, id) in zip(input_pds, subject_ids):
    register.add_polydata(pd, id)

register.points_per_fiber = points_per_fiber

# -------------
# Done SETTINGS. Below is computation
# -------------
total_iterations = numpy.sum(numpy.array(iterations_per_scale))
iteration = 1
# estimate percentage complete based on number of fibers compared,
# because the times cobyla calls the objective function are approx
# constant per scale (except first scale where they are cut short)
total_comparisons = numpy.multiply(iterations_per_scale,numpy.multiply(numpy.array(mean_brain_size_per_scale), numpy.array(subject_brain_size_per_scale)))
total_comparisons = numpy.sum(total_comparisons)
comparisons_so_far = 0

do_scales = range(len(sigma_per_scale))

for scale in do_scales:
    register.sigma = sigma_per_scale[scale]
    register.initial_step = initial_step_per_scale[scale]
    register.final_step = final_step_per_scale[scale]
    register.maxfun = maxfun_per_scale[scale]
    register.mean_brain_size = mean_brain_size_per_scale[scale]
    register.subject_brain_size = subject_brain_size_per_scale[scale]
    
    for idx in range(0,iterations_per_scale[scale]):
        register.iterate()
        comparisons_this_scale = mean_brain_size_per_scale[scale]*subject_brain_size_per_scale[scale]
        comparisons_so_far += comparisons_this_scale
        percent = 100*(float(comparisons_so_far)/total_comparisons)
        print "Done iteration", iteration, "/", total_iterations, ". Percent finished approx:", "%.2f" % percent
        iteration += 1
        # Intermediate save. For testing only.
        if verbose:
            register.save_transformed_polydatas(intermediate_save=True, midsag_symmetric=midsag_symmetric)

# Final save when we are done
register.save_transformed_polydatas(midsag_symmetric=midsag_symmetric)

print "Done registering."


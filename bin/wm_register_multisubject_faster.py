#!/usr/bin/env python
#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/ipython

import argparse
import os
import multiprocessing
import time

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
    help='The mode can be affine or nonrigid. Affine is the default. It should be run first before nonrigid.')
parser.add_argument(
    '-f', action="store", dest="numberOfFibers", type=int, default=20000,
    help='Total number of fibers to analyze from each dataset. During registration, at each iteration fibers are randomly sampled from within this data. 20000 is the default number of total fibers.')
parser.add_argument(
    '-l', action="store", dest="fiberLength", type=int, default=80,
    help='Minimum length (in mm) of fibers to analyze. 60mm is reasonable for DTI single-tensor tractography which is shorter in general. Use a higher value such as 80 or 100 for two-tensor or other advanced tractography. This parameter removes short, noisy fibers and focuses on larger structures that can be registered well. For neonate data, a value of 40mm is suggested. The default is 80mm.')
parser.add_argument(
    '-lmax', action="store", dest="fiberLengthMax", type=int, default=260,
    help='Maximum length (in mm) of fibers to analyze. This parameter can be used to remove extremely long fibers that may have traversed several structures. For example, a value of 200 will avoid sampling the tail end of the fiber length distribution. The default is 260 mm, which is a safe value that should have little effect, as there are few to no fibers expected to be longer than 260mm.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')
parser.add_argument(
    '-verbose', action='store_true', dest="flag_verbose",
    help='Verbose. Run with -verbose to store more files and images of intermediate and final polydatas.')
parser.add_argument(
    '-pf', action="store", dest="pointsPerFiber", type=int, default=15,
    help='Number of points for fiber representation during registration. The default of 15 is reasonable.')
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


points_per_fiber = args.pointsPerFiber
print "<register> Number of points for fiber representation: ", points_per_fiber

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

register.verbose = verbose

register.parallel_jobs = parallel_jobs
register.render = not no_render

# -------------
# SETTINGS
# -------------

# This is the fast default mode using Powell's method for optimization
# and preferring fast iterations with a smaller sample of data
# instead of slower iterations.
if mode == "affine":
    sigma_per_scale = [10, 10, 5]
    iterations_per_scale=[3, 3, 3]
    # Powell seems not to pay much attention to requested max.
    maxfun_per_scale = [20, 40, 80]
    mean_brain_size_per_scale = [3000, 3000, 3500]
    subject_brain_size_per_scale = [1500, 1500, 1500]
    initial_step_per_scale = [10, 5, 3]
    final_step_per_scale = [8, 4, 2]
    # 10ppf is very fast here. <3 minutes per iteration.
    # We override user input here (which was for testing)
    points_per_fiber = 10
    register.nonrigid = False

# This is not particularly better than the above though it is much slower
elif mode == "affine_powellSLOWERTEST":
    sigma_per_scale = [10, 10, 5]
    iterations_per_scale=[3, 1, 3]
    #maxfun_per_scale = [25, 50, 80]
    maxfun_per_scale = [20, 40, 80]
    mean_brain_size_per_scale = [3000, 5000, 5000]
    subject_brain_size_per_scale = [1500, 2500, 2500]
    initial_step_per_scale = [10, 5, 2]
    #final_step_per_scale = [5, 2, 1]
    # try quicker convergence
    final_step_per_scale = [8, 4, 1.5]
    #points_per_fiber = 10
    register.nonrigid = False

# This is now the default fast nonrigid method.
# Points per fiber is set to 5.
# This does multiscale registration using the 5x5x5 fast grid
elif mode == "nonrigid_bspline":
    # testing a fast nonrigid mode to follow after improved affine registration.
    grid_resolution_per_scale = [6, 6]
    # this is in mm space.
    initial_step_per_scale = [2, 1.5]
    final_step_per_scale = [1.5, 1]
    # use only very local information (small sigma)
    sigma_per_scale = [2, 2]
    # how many times to repeat the process at each scale
    iterations_per_scale = [1, 2]
    # Data samples from group model and each brain
    mean_brain_size_per_scale = [3000, 3000]
    subject_brain_size_per_scale = [1000, 1000]
    # 3x3x3 grid, 27*3 = 81 parameter space.
    # 4x4x4 grid, 64*3 = 192 parameter space.
    # 5x5x5 grid, 125*3 = 375 parameter space.
    # 6x6x6 grid, 216*3 = 648 parameter space.
    # Inspection of output pdfs shows that objective decreases steadily for all subjects,
    # so stop the optimizer early and create a better current model.
    #maxfun_per_scale = [400]
    maxfun_per_scale = [700, 1400]
    # fiber representation for computation.
    points_per_fiber = 5
    register.nonrigid = True


# This is now the default fast nonrigid method.
# Points per fiber is set to 5.
# This does multiscale registration using the 5x5x5 fast grid
elif mode == "nonrigid":
    # testing a fast nonrigid mode to follow after improved affine registration.
    grid_resolution_per_scale = [5, 5, 5]
    # this is in mm space.
    initial_step_per_scale = [2, 2, 1]
    final_step_per_scale = [1.5, 1, 0.5]
    # use only very local information (small sigma)
    sigma_per_scale = [2, 2, 1.75]
    # how many times to repeat the process at each scale
    iterations_per_scale = [2, 2, 2]
    # Data samples from group model and each brain
    mean_brain_size_per_scale = [3000, 3000, 4000]
    subject_brain_size_per_scale = [1000, 1500, 2000]
    # 3x3x3 grid, 27*3 = 81 parameter space.
    # 4x4x4 grid, 64*3 = 192 parameter space.
    # 5x5x5 grid, 125*3 = 375 parameter space.
    # 6x6x6 grid, 216*3 = 648 parameter space.
    # Inspection of output pdfs shows that objective decreases steadily for all subjects,
    # so stop the optimizer early and create a better current model.
    maxfun_per_scale = [390, 390, 390]
    # fiber representation for computation.
    points_per_fiber = 5
    register.nonrigid = True

elif mode == "nonrigid_bfgs":
    # this is the fastest possible test if the optimizer is working
    grid_resolution_per_scale = [3,3,3]
    # this is in mm space.
    initial_step_per_scale = [2, 2, 1]
    final_step_per_scale = [1.5, 1, 0.5]
    # use only very local information (small sigma)
    sigma_per_scale = [10, 5, 3]
    #sigma_per_scale = [2, 2, 1.75]
    # how many times to repeat the process at each scale
    iterations_per_scale = [2, 2, 2]
    # Data samples from group model and each brain
    mean_brain_size_per_scale = [2000, 2000, 2000]
    subject_brain_size_per_scale = [750, 750, 750]
    # 3x3x3 grid, 27*3 = 81 parameter space.
    # 4x4x4 grid, 64*3 = 192 parameter space.
    # 5x5x5 grid, 125*3 = 375 parameter space.
    # 6x6x6 grid, 216*3 = 648 parameter space.
    # Inspection of output pdfs shows that objective decreases steadily for all subjects,
    # so stop the optimizer early and create a better current model.
    #maxfun_per_scale = [95, 95, 95]
    #maxfun_per_scale = [200, 200, 200]
    # this is max iterations, where each makes a function call for every dimension in parameter space
    maxfun_per_scale = [2, 2, 2]
    # fiber representation for computation.
    points_per_fiber = 5
    register.nonrigid = True

# This mode also uses the slower-to-compute 6x6x6 grid
elif mode == "nonrigid_finer":
    # testing a fast nonrigid mode to follow after improved affine registration.
    grid_resolution_per_scale = [4, 5, 6, 6]
    # this is in mm space.
    initial_step_per_scale = [2, 2, 1, 1]
    final_step_per_scale = [1.5, 1, 0.8, 0.8]
    # use only very local information (small sigma)
    # going to 1mm so quickly seems not to help. It's not clear 1mm works well. 2mm seems more stable.
    #sigma_per_scale = [2, 2, 1, 1]
    sigma_per_scale = [3, 2, 1.75, 1.5]
    # how many times to repeat the process at each scale
    iterations_per_scale = [3, 3, 2, 2]
    # These are a bit lower than the totals in the affine because
    # this computation is expensive
    mean_brain_size_per_scale = [3000, 3000, 4000, 5000]
    subject_brain_size_per_scale = [1000, 1000, 1500, 2000]
    # 3x3x3 grid, 27*3 = 81 parameter space.
    # 4x4x4 grid, 64*3 = 192 parameter space.
    # 5x5x5 grid, 125*3 = 375 parameter space.
    # 6x6x6 grid, 216*3 = 648 parameter space.
    # Inspection of output pdfs shows that objective decreases steadily for all subjects,
    # so stop the optimizer early and create a better current model.
    maxfun_per_scale = [205, 390, 700, 700]
    # fiber representation for computation.
    points_per_fiber = 5
    register.nonrigid = True

elif mode == "affine_powellOLD":
    sigma_per_scale = [10, 10, 10]
    iterations_per_scale=[3, 3, 3]
    maxfun_per_scale = [25, 50, 80]
    #mean_brain_size_per_scale = [3000, 5000, 6000]
    #subject_brain_size_per_scale = [1500, 2500, 3000]
    mean_brain_size_per_scale = [3000, 5000, 5000]
    subject_brain_size_per_scale = [1500, 2500, 2500]
    initial_step_per_scale = [10, 5, 2]
    #final_step_per_scale = [5, 2, 1]
    # try quicker convergence
    final_step_per_scale = [8, 4, 1.5]
    #points_per_fiber = 10
    register.nonrigid = False

elif mode == "affine_overlapping_cobylaOLD":
    # Mode so we don't waste time computing the first three iterations for connectome data that is rigidly aligned already.
    # Default parameters for affine registration, optimized for speed by stopping computation when most subjects should have improved
    sigma_per_scale = [10, 5, 3, 2]
    iterations_per_scale=[3, 3, 3, 6]
    maxfun_per_scale = [50, 50, 60, 100]
    mean_brain_size_per_scale = [3000, 5000, 6000, 6000]
    subject_brain_size_per_scale = [1500, 2000, 2500, 3000]
    initial_step_per_scale = [5, 2, 1, 0.5]
    final_step_per_scale = [2, 1, 0.5, 0.4]
    register.nonrigid = False

elif mode == "affine_cobylaOLD":
    # Default parameters for affine registration, optimized for speed by stopping computation when most subjects should have improved
    # Note: sigma 1mm is too small for affine: it makes the anterior commissure register worse
    sigma_per_scale = [20, 10, 5, 3, 2]
    iterations_per_scale=[3, 3, 3, 3, 6]
    maxfun_per_scale = [25, 50, 50, 60, 100]
    mean_brain_size_per_scale = [2000, 3000, 5000, 6000, 6000]
    subject_brain_size_per_scale = [500, 1500, 2000, 2500, 3000]
    initial_step_per_scale = [10, 5, 2, 1, 0.5]
    final_step_per_scale = [5, 2, 1, 0.5, 0.4]
    register.nonrigid = False
    #points_per_fiber = 10

elif mode == "nonrigid_powell":
    # testing a fast nonrigid mode to follow after improved affine registration.
    grid_resolution_per_scale = [3, 4, 5, 6]
    # this is in mm space.
    initial_step_per_scale = [2, 2, 1, 0.75]
    final_step_per_scale = [1.5, 1, 0.75, 0.5]
    # use only very local information (small sigma)
    # going to 1mm so quickly seems not to help. It's not clear 1mm works well. 2mm seems more stable.
    #sigma_per_scale = [2, 2, 1, 1]
    #sigma_per_scale = [2, 2, 1.75, 1.5]
    # test larger because powell is more sensitive
    # test smaller because powell is more sensitive
    sigma_per_scale = [2, 2, 2, 2]
    #sigma_per_scale = [5, 3, 2, 1]
    # how many times to repeat the process at each scale
    #iterations_per_scale = [3, 6, 3, 3]
    iterations_per_scale = [1, 1, 1, 1]
    # These are a bit lower than the totals in the affine because
    # this computation is expensive
    mean_brain_size_per_scale = [3000, 3000, 3000, 3000]
    subject_brain_size_per_scale = [1000, 1000, 1000, 1000]
    #mean_brain_size_per_scale = [4500, 5000, 5000, 5000]
    #subject_brain_size_per_scale = [1000, 2000, 2000, 2000]
    # 3x3x3 grid, 27*3 = 81 parameter space.
    # 4x4x4 grid, 64*3 = 192 parameter space.
    # 5x5x5 grid, 125*3 = 375 parameter space.
    # 6x6x6 grid, 216*3 = 648 parameter space.
    # Inspection of output pdfs shows that objective decreases steadily for all subjects,
    # so stop the optimizer early and create a better current model.
    # this is maxiter, which seems to actually stop it (still not)
    maxfun_per_scale = [2, 2, 2, 2]
    #maxfun_per_scale = [100, 210, 390, 700]
    # fiber representation for computation.
    #points_per_fiber = 15
    register.nonrigid = True

elif mode == "nonrigid_powellOLD":
    # testing a fast nonrigid mode to follow after improved affine registration.
    grid_resolution_per_scale = [3, 4, 5, 6]
    # this is in mm space.
    initial_step_per_scale = [2, 2, 1, 0.75]
    final_step_per_scale = [1.5, 1, 0.75, 0.5]
    # use only very local information (small sigma)
    # going to 1mm so quickly seems not to help. It's not clear 1mm works well. 2mm seems more stable.
    #sigma_per_scale = [2, 2, 1, 1]
    #sigma_per_scale = [2, 2, 1.75, 1.5]
    # test larger because powell is more sensitive
    sigma_per_scale = [5, 3, 2, 1]
    # how many times to repeat the process at each scale
    iterations_per_scale = [3, 6, 3, 3]
    # These are a bit lower than the totals in the affine because
    # this computation is expensive
    mean_brain_size_per_scale = [4500, 5000, 5000, 5000]
    subject_brain_size_per_scale = [1000, 2000, 2000, 2000]
    # 3x3x3 grid, 27*3 = 81 parameter space.
    # 4x4x4 grid, 64*3 = 192 parameter space.
    # 5x5x5 grid, 125*3 = 375 parameter space.
    # 6x6x6 grid, 216*3 = 648 parameter space.
    # Inspection of output pdfs shows that objective decreases steadily for all subjects,
    # so stop the optimizer early and create a better current model.
    maxfun_per_scale = [90, 210, 390, 700]
    #maxfun_per_scale = [100, 210, 390, 700]
    # fiber representation for computation.
    #points_per_fiber = 15
    register.nonrigid = True

elif mode == "nonrigidOLDCOBYLA":
    # testing a fast nonrigid mode to follow after improved affine registration.
    grid_resolution_per_scale = [5, 6, 6]
    # this is in mm space.
    initial_step_per_scale = [2, 2, 1]
    final_step_per_scale = [1.5, 1, 1.5]
    # use only very local information (small sigma)
    # going to 1mm so quickly seems not to help. It's not clear 1mm works well. 2mm seems more stable.
    #sigma_per_scale = [2, 2, 1, 1]
    sigma_per_scale = [2, 2, 1]
    # how many times to repeat the process at each scale
    iterations_per_scale = [2, 2, 1]
    # These are a bit lower than the totals in the affine because
    # this computation is expensive
    mean_brain_size_per_scale = [5000, 5000, 6000]
    subject_brain_size_per_scale = [2000, 2000, 3000]
    # 3x3x3 grid, 27*3 = 81 parameter space.
    # 4x4x4 grid, 64*3 = 192 parameter space.
    # 5x5x5 grid, 125*3 = 375 parameter space.
    # 6x6x6 grid, 216*3 = 648 parameter space.
    # Inspection of output pdfs shows that objective decreases steadily for all subjects,
    # so stop the optimizer early and create a better current model.
    maxfun_per_scale = [390, 700, 1200]
    # fiber representation for computation.
    #points_per_fiber = 15
    register.nonrigid = True

elif mode == "nonrigidORIGINAL":
    # testing a fast nonrigid mode to follow after improved affine registration.
    grid_resolution_per_scale = [3, 4, 5, 6]
    # this is in mm space.
    initial_step_per_scale = [2, 2, 1, 0.75]
    final_step_per_scale = [1.5, 1, 0.75, 0.5]
    # use only very local information (small sigma)
    # going to 1mm so quickly seems not to help. It's not clear 1mm works well. 2mm seems more stable.
    #sigma_per_scale = [2, 2, 1, 1]
    sigma_per_scale = [2, 2, 1.75, 1.5]
    # how many times to repeat the process at each scale
    iterations_per_scale = [3, 6, 3, 3]
    # These are a bit lower than the totals in the affine because
    # this computation is expensive
    mean_brain_size_per_scale = [4500, 5000, 5000, 5000]
    subject_brain_size_per_scale = [1000, 2000, 2000, 2000]
    # 3x3x3 grid, 27*3 = 81 parameter space.
    # 4x4x4 grid, 64*3 = 192 parameter space.
    # 5x5x5 grid, 125*3 = 375 parameter space.
    # 6x6x6 grid, 216*3 = 648 parameter space.
    # Inspection of output pdfs shows that objective decreases steadily for all subjects,
    # so stop the optimizer early and create a better current model.
    maxfun_per_scale = [100, 210, 390, 700]
    # fiber representation for computation.
    #points_per_fiber = 15
    register.nonrigid = True
    
elif mode == "affineTEST":
    # very quick test if software is working
    sigma_per_scale = [30, 10, 7.5]
    iterations_per_scale=[1, 1, 1]
    maxfun_per_scale = [60, 80, 100]
    mean_brain_size_per_scale = [1500, 2000, 3000]
    subject_brain_size_per_scale = [100, 500, 1000]
    initial_step_per_scale = [5, 5, 5, 5]
    final_step_per_scale = [2, 2, 2, 2]
    register.nonrigid = False
    points_per_fiber = 5
    
elif mode == "nonrigidTEST":
    # very quick test if software is working
    grid_resolution_per_scale = [3, 4, 5, 6, 8, 10]
    initial_step_per_scale = [5, 3, 1, 1, 1, 1]
    final_step_per_scale = [2, 1, 1, 1, 1, 1]
    sigma_per_scale = [3, 2, 1, 1, 1, 1]
    iterations_per_scale = [1, 1, 1, 1, 1, 1]
    mean_brain_size_per_scale = [1000, 1000, 1000, 1000, 1000, 1000]
    subject_brain_size_per_scale = [100, 100, 100, 100, 100, 100]
    # stop computation: this is just a quick test the software is working
    maxfun_per_scale = [10, 10, 10, 10, 10, 10]
    points_per_fiber = 15
    register.nonrigid = True

else:
    print "\n<register> Error: Unknown registration mode:", mode
    exit()

# We have to add polydatas after setting nonrigid in the register object
for (pd, id) in zip(input_pds, subject_ids):
    register.add_polydata(pd, id)

register.points_per_fiber = points_per_fiber

# output summary file to save information about what was run
readme_fname = os.path.join(args.outputDirectory, 'README.txt')
readme_file = open(readme_fname, 'w')
outstr = "Groupwise Registration Summary\n"
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
outstr += 'registration_atlas.vtk: This is the template created from groupwise registration.\n'
outstr += 'registration_performance.txt: Parameters and objective values from each iteration.\n'
outstr += 'progress.txt: Update of how far registration has progressed.\n'
outstr += 'input_subjects.txt:  List of subject index, ID, and full path to input file.\n'
outstr += 'README.txt:  This summary file.\n'
outstr += 'output_tractography/\n'
outstr += '\tThe output directory with registered tractography and corresponding transforms.\n'
outstr += 'The files inside each iteration directory are for testing purposes:\n'
outstr += '\titk_txform files output from that iteration.\n'
outstr += '\tpdf files plot objective function changes for all subjects.\n'
outstr += '\n'
outstr += '\n'
outstr += "Command Line Arguments\n"
outstr += '----------------------\n'
outstr += str(args)
outstr += '\n'
outstr += '\n'
outstr += "Parameters\n"
outstr += '----------------------\n'
outstr += "Registration mode: " + mode
outstr += '\n'
outstr += "Number of fibers to analyze per subject: " + str(number_of_fibers)
outstr += '\n'
outstr += "Minimum length of fibers to analyze (in mm): " + str(fiber_length)
outstr += '\n'
outstr += "Maximum  length of fibers to analyze (in mm): " + str(fiber_length_max)
outstr += '\n'
outstr += "Number of jobs to use: " + str(parallel_jobs)
outstr += '\n'
outstr += "verbose: " + str(verbose)
outstr += '\n'
outstr += "render: " + str(not no_render)
outstr += '\n'
outstr += "midsag_symmetric: " + str(midsag_symmetric)
outstr += '\n'
outstr += "random seed: " + str(random_seed)
outstr += '\n'
outstr += '\n'
outstr += "Input Fiber Files\n"
outstr += '-----------------\n'
for pd in input_polydatas:
    outstr += pd
    outstr += '\n'
readme_file.write(outstr)
readme_file.close()

# output summary file to save information about all subjects
subjects_qc_fname = os.path.join(args.outputDirectory, 'input_subjects.txt')
subjects_qc_file = open(subjects_qc_fname, 'w')
outstr = "Subject_idx\tSubject_ID\tInput Filename\n"
subjects_qc_file.write(outstr)
idx = 1
for fname in input_polydatas:
    subject_id = os.path.splitext(os.path.basename(fname))[0]
    outstr =  str(idx) + '\t' + str(subject_id) + '\t' + str(fname) + '\n'
    subjects_qc_file.write(outstr)
    idx += 1
subjects_qc_file.close()

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
progress_filename = os.path.join(args.outputDirectory, 'progress.txt')
progress_file = open(progress_filename, 'w')
print >> progress_file, "Beginning registration. Total iterations will be:", total_iterations
print >> progress_file,"Start date: "  + time.strftime("%x")
print >> progress_file, "Start time: " + time.strftime("%X") + '\n'
progress_file.close()
prev_time = time.time()
do_scales = range(len(sigma_per_scale))

for scale in do_scales:
    register.sigma = sigma_per_scale[scale]
    register.initial_step = initial_step_per_scale[scale]
    register.final_step = final_step_per_scale[scale]
    register.maxfun = maxfun_per_scale[scale]
    register.mean_brain_size = mean_brain_size_per_scale[scale]
    register.subject_brain_size = subject_brain_size_per_scale[scale]
    if register.nonrigid:
        register.nonrigid_grid_resolution = grid_resolution_per_scale[scale]
        register.update_nonrigid_grid()
    
    for idx in range(0,iterations_per_scale[scale]):
        register.iterate()
        comparisons_this_scale = mean_brain_size_per_scale[scale]*subject_brain_size_per_scale[scale]
        comparisons_so_far += comparisons_this_scale
        percent = 100*(float(comparisons_so_far)/total_comparisons)
        print "Done iteration", iteration, "/", total_iterations, ". Percent finished approx:", "%.2f" % percent
        progress_file = open(progress_filename, 'a')
        curr_time = time.time()
        print >> progress_file, "Done iteration", iteration, "/", total_iterations, ". Percent finished approx:", "%.2f" % percent, ". Time:", time.strftime("%X"), ". Minutes Elapsed:", (curr_time - prev_time)/60
        progress_file.close()
        prev_time = curr_time

        iteration += 1
        # Intermediate save. For testing only.
        if verbose:
            register.save_transformed_polydatas(intermediate_save=True, midsag_symmetric=midsag_symmetric)

# Final save when we are done
register.save_transformed_polydatas(midsag_symmetric=midsag_symmetric)

print "\nDone registering. For more information on the output, please read:", readme_fname, "\n"

progress_file = open(progress_filename, 'a')
print >> progress_file, "\nFinished registration."
print >> progress_file,"End date: "  + time.strftime("%x")
print >> progress_file, "End time: " + time.strftime("%X")
progress_file.close()

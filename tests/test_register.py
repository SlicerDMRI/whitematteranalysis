#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/ipython

# Run registration on the test dataset.

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

#indir1 = 'test_data'
#outdir = 'test_register_results'

# defaults that may be added as parameters later
#fiber_sample_sizes = [25, 50, 75, 100]
fiber_sample_sizes = [50, 150, 200, 200]
sigma_per_step = [30, 10, 10, 5]
# not used
#maxfun = 300
## for multiple subjects this is good. for two, it's not enough.
##minfun = number_of_datasets * 3
##maxfun_per_step = [minfun*1.5, minfun*2, minfun*5, minfun*10]

#maxfun_per_step = [50, 75, 200]
# this is reasonable for two subjects, except for shear.
# why does shear not always converge? what is rho set as?
maxfun_per_step = [20, 40, 60, 80]


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

if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

if args.numberOfFibers is not None:
    number_of_fibers = args.numberOfFibers
else:
    number_of_fibers = 300
print "fibers to analyze per subject: ", number_of_fibers


if args.fiberLength is not None:
    fiber_length = args.fiberLength
else:
    fiber_length = 75
print "minimum length of fibers to analyze (in mm): ", fiber_length
    
print 'CPUs detected:', multiprocessing.cpu_count()
if args.numberOfJobs is not None:
    parallel_jobs = args.numberOfJobs
else:
    parallel_jobs = multiprocessing.cpu_count()
print 'Using N jobs:', parallel_jobs

if args.flag_verbose:
    print "Verbose display and intermediate image saving ON."
else:
    print "Verbose display and intermediate image saving OFF."
verbose = args.flag_verbose

if args.pointsPerFiber is not None:
    print "Number of points for fiber representation: ", args.pointsPerFiber
    points_per_fiber = args.pointsPerFiber
else:
    points_per_fiber = 5


print "<wm_register.py> Starting registration computation."
print ""
print "=====input directory ======\n", args.inputDirectory
print "=====output directory =====\n", args.outputDirectory
print "=========================="

# Find input files
inputMask1 = "{0}/*.vtk".format(args.inputDirectory)
inputMask2 = "{0}/*.vtp".format(args.inputDirectory)
input_poly_datas = glob.glob(inputMask1) + glob.glob(inputMask2)
print "<wm_register.py> Input number of files: ", len(input_poly_datas)
print input_poly_datas

def register_scale_step(register, scale_mode, n_steps):

    if scale_mode == "Coarse":
        # n = 5
        # only translation and rotation. initialization.
        for idx in range(0, n_steps):
            register.translate_only()
            register.compute()
            register.rotate_only()
            register.compute()
    elif scale_mode == "Medium":
        # n = 1
        for idx in range(0, n_steps):
            register.translate_only()
            register.compute()
            register.rotate_only()
            register.compute()
            register.scale_only()
            register.compute()
            register.shear_only()
            register.compute()
    elif scale_mode == "Fine":
        # n = 1
        for idx in range(0, n_steps):
            register.translate_only()
            register.compute()
            register.rotate_only()
            register.compute()
            register.scale_only()
            register.compute()
            register.shear_only()
            register.compute()
    elif scale_mode == "Finest":
        # n = 1
        for idx in range(0, n_steps):
            register.translate_only()
            register.compute()
            register.rotate_only()
            register.compute()
            register.scale_only()
            register.compute()
            register.shear_only()
            register.compute()
            
    
def run_registration(input_poly_datas, outdir, number_of_fibers=150,
    fiber_sample_sizes=[75, 75, 75, 100],
    parallel_jobs=2,
    points_per_fiber=5,
    sigma_per_step=[30, 10, 10, 5],
    maxfun_per_step=[10, 40, 60, 80],
    distance_method='Hausdorff', verbose=True, fiber_length=75):

    elapsed = list()

    number_of_datasets = len(input_poly_datas)

    print 'Read and preprocess'
    input_pds = list()
    for fname in input_poly_datas:
        print fname
        pd = wma.io.read_polydata(fname)
        pd2 = wma.filter.preprocess(pd, fiber_length)
        pd3 = wma.filter.downsample(pd2, number_of_fibers)
        input_pds.append(pd3)

    # create registration object and apply settings
    register = wma.congeal.CongealTractography()
    register.parallel_jobs = parallel_jobs
    register.threshold = 0
    register.points_per_fiber = points_per_fiber
    register.distance_method = distance_method
    
    # add inputs to the registration
    for pd in input_pds:
        register.add_subject(pd)

    # view output data from the initialization
    outdir_current =  os.path.join(outdir, 'iteration_0')
    if not os.path.exists(outdir_current):
        os.makedirs(outdir_current)
    output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
    ren = wma.registration_functions.view_polydatas(output_pds)
    ren.save_views(outdir_current)
    del ren
    wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir_current)
    
    scales = ["Coarse", "Medium", "Fine", "Finest"]
    #steps_per_scale = [5, 1, 1, 1]
    steps_per_scale = [10, 3, 2, 2]
    scale_idx = 0
    for scale in scales:
        start = time.time()
        # run the basic iteration of translate, rotate, scale
        register.fiber_sample_size = fiber_sample_sizes[scale_idx]
        register.sigma = sigma_per_step[scale_idx]
        register.maxfun = maxfun_per_step[scale_idx]
        scale_mode = scales[scale_idx]
        register_scale_step(register, scale_mode, steps_per_scale[scale_idx])
        elapsed.append(time.time() - start)
        scale_idx += 1
        
        # view output data from this big iteration
        if verbose | (scale_idx == 4):
            outdir_current =  os.path.join(outdir, 'iteration_'+str(scale_idx))
            if not os.path.exists(outdir_current):
                os.makedirs(outdir_current)
            output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
            ren = wma.registration_functions.view_polydatas(output_pds)
            ren.save_views(outdir_current)
            del ren
            if scale_idx == 4:
                wma.registration_functions.transform_polydatas_from_disk(input_poly_datas, register, outdir_current)
            wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir_current)
    
            plt.figure() # to avoid all results on same plot
            plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
            plt.savefig(os.path.join(outdir_current, 'objective_function.pdf'))
        
    return register, elapsed

## run the registration ONCE and output result to disk
register, elapsed = run_registration(input_poly_datas, outdir,
                        number_of_fibers=number_of_fibers,
                        points_per_fiber=points_per_fiber,
                        parallel_jobs=parallel_jobs,
                        fiber_sample_sizes=fiber_sample_sizes,
                        sigma_per_step=sigma_per_step,
                        maxfun_per_step=maxfun_per_step,
                        verbose=verbose,
                        fiber_length=fiber_length)

print "TIME:", elapsed

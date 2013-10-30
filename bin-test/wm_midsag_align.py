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
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu",
    version='1.0')

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
    help='Minimum length (in mm) of fibers to analyze. 60mm is default.')
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
    fiber_length = 60

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

def run_registration(input_poly_data, outdir, number_of_fibers=150,
    number_of_fibers_per_step=[75, 75, 75, 100],
    parallel_jobs=2,
    points_per_fiber=5,
    sigma_per_step=[30, 10, 10, 5],
    maxfun=150,
    distance_method='Hausdorff',
    fiber_length=60,
    verbose=True):
    
    elapsed = list()
    number_of_fibers_step_one = number_of_fibers_per_step[0]
    number_of_fibers_step_two = number_of_fibers_per_step[1]
    number_of_fibers_step_three = number_of_fibers_per_step[2]
    number_of_fibers_step_four = number_of_fibers_per_step[3]

    minfun = maxfun/5.0
    #maxfun_per_step = [50, 75, 200]
    maxfun_per_step = [minfun*1.5, minfun*2, minfun*3, minfun*5]

    sigma_step_one = sigma_per_step[0]
    sigma_step_two = sigma_per_step[1]
    sigma_step_three = sigma_per_step[2]
    sigma_step_four = sigma_per_step[3]

    print 'Read and preprocess'
    input_pds = list()

    pd = wma.io.read_polydata(input_poly_data)
    pd2 = wma.filter.preprocess(pd, fiber_length)
    pd3 = wma.filter.downsample(pd2, number_of_fibers)
    pd4 = wma.filter.downsample(pd2, number_of_fibers)
    trans = vtk.vtkTransform()
    trans.Scale(-1,1,1)
    transformer = vtk.vtkTransformPolyDataFilter()
    transformer.SetTransform(trans)
    transformer.SetInputData(pd4)
    transformer.Update()
    pd5 = transformer.GetOutput()
    # append input and its reflection to list of data to register
    input_pds.append(pd3)
    input_pds.append(pd5)

    # create registration object and apply settings
    register = wma.congeal.CongealTractography()
    register.parallel_jobs = parallel_jobs
    register.threshold = 0
    register.points_per_fiber = points_per_fiber
    register.distance_method = distance_method
    #register.maxfun = maxfun
    # make sure we take very small steps, the brains are already overlapping
    inc_rot = (0.5 / 180.0) * numpy.pi
    inc_trans = 0.5
    inc_scale = 0.001
    inc_shear = (.5 / 180.0) * numpy.pi
    register.set_rhobeg(inc_rot, inc_trans, inc_scale, inc_shear)
    
    # add inputs to the registration
    for pd in input_pds:
        register.add_subject(pd)

    # view input data from the initialization
    if verbose:
        outdir_current =  os.path.join(outdir, 'iteration_0')
        if not os.path.exists(outdir_current):
            os.makedirs(outdir_current)
        output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
        ren = wma.registration_functions.view_polydatas(output_pds, 1000)
        ren.save_views(outdir_current)
        #wma.registration_functions.transform_polydatas_from_disk(input_poly_data, register, outdir_current)
        wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir)
    
    # STEP ONE
    start = time.time()
    # run the basic iteration of translate, rotate, scale
    register.fiber_sample_size = number_of_fibers_step_one
    register.sigma = sigma_step_one
    register.maxfun = maxfun_per_step[0]
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    # Don't scale for midsagittal alignment
    #register.scale_only()
    #register.compute()
    elapsed.append(time.time() - start)

    # STEP TWO
    start = time.time()
    # run the basic iteration of translate, rotate, scale AGAIN
    register.fiber_sample_size = number_of_fibers_step_two
    register.sigma = sigma_step_two
    register.maxfun = maxfun_per_step[1]    
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    #register.scale_only()
    register.compute()
    register.shear_only()
    register.compute()
    elapsed.append(time.time() - start)

    #if 0:
    # STEP THREE
    start = time.time()
    # run the basic iteration of translate, rotate, scale AGAIN
    register.fiber_sample_size = number_of_fibers_step_three
    register.sigma = sigma_step_three
    register.maxfun = maxfun_per_step[2]
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    #register.scale_only()
    register.compute()
    register.shear_only()
    register.compute()
    elapsed.append(time.time() - start)

    # STEP FOUR
    start = time.time()
    # run the basic iteration of translate, rotate, scale AGAIN
    register.fiber_sample_size = number_of_fibers_step_four
    register.sigma = sigma_step_four
    register.maxfun = maxfun_per_step[3]
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    #register.scale_only()
    register.compute()
    register.shear_only()
    register.compute()
    elapsed.append(time.time() - start)
    
    # view output data from this big iteration
    if verbose:
        outdir_current =  os.path.join(outdir, 'iteration_4')
        if not os.path.exists(outdir_current):
            os.makedirs(outdir_current)
        output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
        ren = wma.registration_functions.view_polydatas(output_pds, 500)
        ren.save_views(outdir_current)
        #wma.registration_functions.transform_polydatas_from_disk(input_poly_data, register, outdir_current)
        wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir)
    
        plt.figure() # to avoid all results on same plot
        plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
        plt.savefig(os.path.join(outdir_current, 'objective_function.pdf'))

    return register, elapsed, input_pds

## run the registration ONCE and output result to disk
points_per_fiber = 5
number_of_fibers_per_step = [100, 200, 200, 250]

# for the cluster
#number_of_fibers_per_step = [300, 300, 500, 1000]
# small sigmas only, this is a minor adjustment 
sigma_per_step = [10, 10, 5, 5]
maxfun = 30
# output location

# registration
register, elapsed, input_pds = run_registration(input_poly_data, outdir,
                                                number_of_fibers=number_of_fibers,
                                                points_per_fiber=points_per_fiber,
                                                parallel_jobs=parallel_jobs,
                                                number_of_fibers_per_step=number_of_fibers_per_step,
                                                sigma_per_step=sigma_per_step,
                                                maxfun=maxfun,
                                                fiber_length=fiber_length,
                                                verbose=verbose)

print "TIME:", elapsed

print "Re-reading and transforming original data to aligned. Writing outputs."

# now apply the appropriate transform to the input data.
# half of transform 1 times transform 2 inverse

tx = register.convert_transforms_to_vtk()
tx[1].Inverse()
tx[0].Concatenate(tx[1])

txs = vtk.vtkTransform()
m = vtk.vtkMatrix4x4()
for i in range(0,4):
    for j in range(0,4):
        # scale is 1, and value 3,3 is 1
        if i == j:
            m.SetElement(i,j, 1.0)
        else:
            el = tx[0].GetMatrix().GetElement(i,j)
            m.SetElement(i,j, el / 2.0)

txs.SetMatrix(m)

pd = wma.io.read_polydata(input_poly_data)

trans = vtk.vtkTransformPolyDataFilter()
trans.SetTransform(txs)
trans.SetInputData(pd)
trans.Update()

outdir_current =  os.path.join(outdir, 'output')
if not os.path.exists(outdir_current):
    os.makedirs(outdir_current)

fname1 = os.path.split(input_poly_data)[1]
fname1 = os.path.splitext(fname1)[0]
fname1 = os.path.join(outdir_current, fname1+'_sym.vtp')
#wma.io.write_polydata(trans.GetOutput(), os.path.join(outdir_current,'symmetric.vtk'))
print "Writing output polydata..."
wma.io.write_polydata(trans.GetOutput(), fname1)

# make the picture. View input, reflection, and output
all_pds = list()
all_pds.append(input_pds[0])
all_pds.append(input_pds[1])
all_pds.append(trans.GetOutput())
ren = wma.registration_functions.view_polydatas(all_pds, 500)
ren.save_views(outdir_current)
    




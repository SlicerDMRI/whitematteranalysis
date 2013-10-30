#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/ipython

# Run registration on the test dataset.

import os
import glob
import matplotlib.pyplot as plt
import numpy
import time

import vtk

import whitematteranalysis as wma

import multiprocessing
parallel_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', parallel_jobs
parallel_jobs = 4
print 'Using N jobs:', parallel_jobs

indir1 = 'test_data'
outdir = 'test_register_results'

input_mask1 = "{0}/*.vtk".format(indir1)
input_poly_datas = glob.glob(input_mask1)

print input_poly_datas
    
def run_registration(input_poly_datas, outdir, number_of_fibers=150,
    number_of_fibers_per_step=[75, 75, 75, 100],
    parallel_jobs=2,
    points_per_fiber=5,
    sigma_per_step=[30, 10, 10, 5],
    maxfun=150,
    distance_method='Hausdorff'):

    elapsed = list()
    minimum_length = 60
    number_of_fibers_step_one = number_of_fibers_per_step[0]
    number_of_fibers_step_two = number_of_fibers_per_step[1]
    number_of_fibers_step_three = number_of_fibers_per_step[2]
    number_of_fibers_step_four = number_of_fibers_per_step[3]

    number_of_datasets = len(input_poly_datas)
    minfun = number_of_datasets * 3
    #maxfun_per_step = [50, 75, 200]
    maxfun_per_step = [minfun*1.5, minfun*2, minfun*5, minfun*10]

    sigma_step_one = sigma_per_step[0]
    sigma_step_two = sigma_per_step[1]
    sigma_step_three = sigma_per_step[2]
    sigma_step_four = sigma_per_step[3]

    print 'Read and preprocess'
    input_pds = list()
    for fname in input_poly_datas:
        print fname
        pd = wma.io.read_polydata(fname)
        pd2 = wma.filter.preprocess(pd, 60)
        pd3 = wma.filter.downsample(pd2, number_of_fibers)
        input_pds.append(pd3)

    # view input data
    #ren = wma.registration_functions.view_polydatas(input_pds)

    # create registration object and apply settings
    register = wma.congeal.CongealTractography()
    register.parallel_jobs = parallel_jobs
    register.threshold = 0
    register.points_per_fiber = points_per_fiber
    register.distance_method = distance_method
    register.maxfun = maxfun
    
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
    #wma.registration_functions.transform_polydatas_from_disk(input_poly_datas, register, outdir_current)
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
    # Don't scale at the first step
    #register.scale_only()
    #register.compute()
    elapsed.append(time.time() - start)

    # view output data from this big iteration
    outdir_current =  os.path.join(outdir, 'iteration_1')
    if not os.path.exists(outdir_current):
        os.makedirs(outdir_current)
    output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
    ren = wma.registration_functions.view_polydatas(output_pds)
    ren.save_views(outdir_current)
    #wma.registration_functions.transform_polydatas_from_disk(input_poly_datas, register, outdir_current)
    wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir)
    
    plt.figure() # to avoid all results on same plot
    #plt.plot(range(len(register.objective_function_values)), numpy.log(register.objective_function_values))
    plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
    plt.savefig(os.path.join(outdir_current, 'objective_function.pdf'))

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
    register.scale_only()
    register.compute()
    register.shear_only()
    register.compute()
    elapsed.append(time.time() - start)

    # view output data from this big iteration
    outdir_current =  os.path.join(outdir, 'iteration_2')
    if not os.path.exists(outdir_current):
        os.makedirs(outdir_current)
    output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
    ren = wma.registration_functions.view_polydatas(output_pds)
    ren.save_views(outdir_current)
    #wma.registration_functions.transform_polydatas_from_disk(input_poly_datas, register, outdir_current)
    wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir)
    
    plt.figure() # to avoid all results on same plot
    #plt.plot(range(len(register.objective_function_values)), numpy.log(register.objective_function_values))
    plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
    plt.savefig(os.path.join(outdir_current, 'objective_function.pdf'))

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
    register.scale_only()
    register.compute()
    register.shear_only()
    register.compute()
    elapsed.append(time.time() - start)

    # view output data from this big iteration
    outdir_current =  os.path.join(outdir, 'iteration_3')
    if not os.path.exists(outdir_current):
        os.makedirs(outdir_current)
    output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
    ren = wma.registration_functions.view_polydatas(output_pds)
    ren.save_views(outdir_current)
    #wma.registration_functions.transform_polydatas_from_disk(input_poly_datas, register, outdir_current)
    wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir)
    
    plt.figure() # to avoid all results on same plot
    #plt.plot(range(len(register.objective_function_values)), numpy.log(register.objective_function_values))
    plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
    plt.savefig(os.path.join(outdir_current, 'objective_function.pdf'))

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
    register.scale_only()
    register.compute()
    register.shear_only()
    register.compute()
    elapsed.append(time.time() - start)
    
    # view output data from this big iteration
    outdir_current =  os.path.join(outdir, 'iteration_4')
    if not os.path.exists(outdir_current):
        os.makedirs(outdir_current)
    output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
    ren = wma.registration_functions.view_polydatas(output_pds)
    ren.save_views(outdir_current)
    wma.registration_functions.transform_polydatas_from_disk(input_poly_datas, register, outdir_current)
    wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir)
    
    plt.figure() # to avoid all results on same plot
    #plt.plot(range(len(register.objective_function_values)), numpy.log(register.objective_function_values))
    plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
    plt.savefig(os.path.join(outdir_current, 'objective_function.pdf'))



    return register, elapsed

## run the registration ONCE and output result to disk
number_of_fibers = 400
points_per_fiber = 5
number_of_fibers_per_step = [25, 50, 75, 100]
sigma_per_step = [30, 10, 10, 5]
maxfun = 300
# output location
if not os.path.exists(outdir):
    os.makedirs(outdir)

# registration
register, elapsed = run_registration(input_poly_datas, outdir,
                        number_of_fibers=number_of_fibers,
                        points_per_fiber=points_per_fiber,
                        parallel_jobs=parallel_jobs,
                        number_of_fibers_per_step=number_of_fibers_per_step,
                        sigma_per_step=sigma_per_step,
                        maxfun=maxfun\
                        )

print "TIME:", elapsed

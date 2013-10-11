# create synthetic data
# load one subject, with different random samples of fibers
# apply some random transformations (only translate, rotate, scale)
# run registration
# compute errors in trans, rot, scale

import matplotlib.pyplot as plt
import numpy
import os

import vtk

import whitematteranalysis as wma

# SETTINGS
range_trans = [-15, 15]
range_scale = [0.85, 1.15]
range_rot = [-20, 20]

outdir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/MICCAI2012/synthetic_data/synth_data_feb15c'
number_of_datasets = 10
parallel_jobs = 25
print 'Using N jobs:', parallel_jobs
number_of_fibers = 300
points_per_fiber = 5
number_of_fibers_per_step = [25, 50, 75]
sigma_per_step = [30, 10, 10]
minfun = number_of_datasets * 3
#maxfun_per_step = [50, 75, 200]
maxfun_per_step = [minfun*1.5, minfun*2, minfun*5]
# output location
if not os.path.exists(outdir):
    os.makedirs(outdir)


fname = '/Users/odonnell/Dropbox/Data/TBI_FE_PNL/controls/01231-dwi-filt-Ed-DTI-tract.vtp'


# data creation
pd_original = wma.io.read_polydata(fname)
# within subject 30mm gives lower error than 60mm
#fiber_length = 30
#fiber_length = 60
fiber_length = 40
number_of_fibers = 300

synthetic_data_list = list()
for idx in range(number_of_datasets):
    synthetic_data_list.append(
    wma.filter.downsample(
        wma.filter.preprocess(pd_original, 
                              fiber_length), 
                        number_of_fibers
        ))
        

input_pds = list()
tx_trans_list = list()
tx_rot_list = list()
tx_scale_list = list()
for pd in synthetic_data_list:
    vtktrans = vtk.vtkTransform()
    vtktrans.Identity()
    tx = numpy.random.random(3)*(range_trans[1] - range_trans[0]) + range_trans[0]
    tx_trans_list.append(tx)
    vtktrans.Translate(tx[0], tx[1], tx[2])
    tx = numpy.random.random(3)*(range_rot[1] - range_rot[0]) + range_rot[0]
    tx_rot_list.append(tx)
    vtktrans.RotateX(tx[0])
    vtktrans.RotateY(tx[1])
    vtktrans.RotateZ(tx[2])
    tx = numpy.random.random(3)*(range_scale[1] - range_scale[0]) + range_scale[0]
    tx_scale_list.append(tx)
    vtktrans.Scale(tx[0], tx[1], tx[2])
    
    transformer = vtk.vtkTransformPolyDataFilter()
    transformer.SetInput(pd)
    transformer.SetTransform(vtktrans)
    transformer.Update()
    input_pds.append(transformer.GetOutput())
    
# create registration object and apply settings
register = wma.congeal.CongealTractography()
register.parallel_jobs = parallel_jobs
register.threshold = 0
register.points_per_fiber = points_per_fiber
register.distance_method = "Hausdorff"

    
# add inputs to the registration
for pd in input_pds:
    register.add_subject(pd)

# view synthetic data
outdir_current =  os.path.join(outdir, 'iteration_0')
if not os.path.exists(outdir_current):
    os.makedirs(outdir_current)
output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
ren = wma.registration_functions.view_polydatas(output_pds)
ren.save_views(outdir_current)
    
# run the registration
if 1:
    number_of_fibers_step_one = number_of_fibers_per_step[0]
    number_of_fibers_step_two = number_of_fibers_per_step[1]
    number_of_fibers_step_three = number_of_fibers_per_step[2]

    sigma_step_one = sigma_per_step[0]
    sigma_step_two = sigma_per_step[1]
    sigma_step_three = sigma_per_step[2]

    maxfun_step_one = maxfun_per_step[0]
    maxfun_step_two = maxfun_per_step[1]
    maxfun_step_three = maxfun_per_step[2]

    # STEP ONE
    # run the basic iteration of translate, rotate, scale, twice
    register.fiber_sample_size = number_of_fibers_step_one
    register.sigma = sigma_step_one
    register.maxfun = maxfun_step_one
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
    register.maxfun = maxfun_step_two    
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

    # view output data from this big iteration
    outdir_current =  os.path.join(outdir, 'iteration_1')
    if not os.path.exists(outdir_current):
        os.makedirs(outdir_current)
    output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
    ren = wma.registration_functions.view_polydatas(output_pds)
    ren.save_views(outdir_current)
    #wma.registration_functions.transform_polydatas_from_disk(input_poly_datas, register, outdir_current)
    
    plt.figure() # to avoid all results on same plot
    #plt.plot(range(len(register.objective_function_values)), numpy.log(register.objective_function_values))
    plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
    plt.savefig(os.path.join(outdir_current, 'objective_function.pdf'))

    # STEP TWO
    # run the basic iteration of translate, rotate, scale AGAIN
    register.fiber_sample_size = number_of_fibers_step_two
    register.sigma = sigma_step_two
    register.maxfun = maxfun_step_three
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.scale_only()
    register.compute()
    # run again
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.scale_only()
    register.compute()

    # view output data from this big iteration
    outdir_current =  os.path.join(outdir, 'iteration_2')
    if not os.path.exists(outdir_current):
        os.makedirs(outdir_current)
    output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
    ren = wma.registration_functions.view_polydatas(output_pds)
    ren.save_views(outdir_current)
    #wma.registration_functions.transform_polydatas_from_disk(input_poly_datas, register, outdir_current)
    
    plt.figure() # to avoid all results on same plot
    #plt.plot(range(len(register.objective_function_values)), numpy.log(register.objective_function_values))
    plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
    plt.savefig(os.path.join(outdir_current, 'objective_function.pdf'))

    #if 0:
    # STEP THREE
    # run the basic iteration of translate, rotate, scale AGAIN
    register.fiber_sample_size = number_of_fibers_step_three
    register.sigma = sigma_step_three
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.scale_only()
    register.compute()

    # view output data from this big iteration
    outdir_current =  os.path.join(outdir, 'iteration_3')
    if not os.path.exists(outdir_current):
        os.makedirs(outdir_current)
    output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
    ren = wma.registration_functions.view_polydatas(output_pds)
    ren.save_views(outdir_current)
    #wma.registration_functions.transform_polydatas_from_disk(input_poly_datas, register, outdir_current)
    
    plt.figure() # to avoid all results on same plot
    #plt.plot(range(len(register.objective_function_values)), numpy.log(register.objective_function_values))
    plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
    plt.savefig(os.path.join(outdir_current, 'objective_function.pdf'))

calculated_transforms = list()
for subj in register._subjects:
    calculated_transforms.append(subj.transform)

gold_standard = list()
for idx in range(number_of_datasets):
    rand_trans = numpy.array([tx_rot_list[idx]*numpy.pi/180.0, tx_trans_list[idx], tx_scale_list[idx]]).ravel()
    gold_standard.append(rand_trans)

calculated_transforms = numpy.array(calculated_transforms)
gold_standard = numpy.array(gold_standard)

# remove the mean transform from the gold standard
meantrans = numpy.mean(gold_standard,0)
gold_standard[:,0:6] -= meantrans[0:6]
gold_standard[:,6:9] = numpy.divide(gold_standard[:,6:9], meantrans[6:9])

#for idx in range(number_of_datasets):
#    print gold_standard[idx,0:3] + calculated_transforms[idx,0:3]

error = gold_standard + calculated_transforms
#error[:,6:9] = numpy.multiply(gold_standard[:,6:9], calculated_transforms[:,6:9])
error[:,6:9] = numpy.divide(1.0, gold_standard[:,6:9]) - calculated_transforms[:,6:9]

print numpy.round((error*100))/100.0

# above errors are pretty low. see if we can do better,
# is the information there?
register.fiber_sample_size = 200
register.sigma = sigma_step_three
register.translate_only()
register.compute()
register.rotate_only()
register.compute()
register.scale_only()
register.compute()


calculated_transforms2 = list()
for subj in register._subjects:
    calculated_transforms2.append(subj.transform)
calculated_transforms2 = numpy.array(calculated_transforms2)
error2 = gold_standard + calculated_transforms2
#error[:,6:9] = numpy.multiply(gold_standard[:,6:9], calculated_transforms[:,6:9])
error2[:,6:9] = numpy.divide(1.0, gold_standard[:,6:9]) - calculated_transforms2[:,6:9]

print numpy.round((error2*100))/100.0

# view output data from this big iteration
outdir_current =  os.path.join(outdir, 'iteration_4')
if not os.path.exists(outdir_current):
    os.makedirs(outdir_current)
output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
ren = wma.registration_functions.view_polydatas(output_pds)
ren.save_views(outdir_current)

#meantrans = numpy.mean(calculated_transforms2,0)
#calculated_transforms2[:,0:6] -= meantrans[0:6]
#calculated_transforms2[:,6:9] = numpy.divide(calculated_transforms2[:,6:9], meantrans[6:9])


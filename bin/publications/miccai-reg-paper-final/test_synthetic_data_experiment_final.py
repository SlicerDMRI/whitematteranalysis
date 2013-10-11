# create synthetic data
# load one subject, with different random samples of fibers
# apply some random transformations (only translate, rotate, scale)
# run registration
# compute errors in trans, rot, scale
import time
import matplotlib.pyplot as plt
import numpy
import os

import vtk

import whitematteranalysis as wma

# SETTINGS
#range_trans = [-15, 15]
range_trans = [-20, 20]
range_scale = [0.85, 1.15]
range_rot = [-20, 20]
#range_shear = [-10, 10]

#outdir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/MICCAI2012/synthetic_data/synth_data_feb16prob'
#outdir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/MICCAI2012/synthetic_data/synth_data_feb22_10subj_shear'
#outdir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/MICCAI2012/synthetic_data/synth_data_feb24_10subj_noshear'
outdir = '/Users/odonnell/Dropbox/Work/Publications/Results/2012_MICCAI_reg_results/synthetic_data/synth_data_jun07_10subj_noshear_5pts_25_50_75'

number_of_datasets = 10
#number_of_datasets = 4
parallel_jobs = 25
print 'Using N jobs:', parallel_jobs
points_per_fiber = 5
#number_of_fibers_per_step = [25, 50, 75, 100]
number_of_fibers_per_step = [25, 50, 50, 75]
#sigma_per_step = [30, 10, 10]
sigma_per_step = [30, 10, 10, 5]
#sigma_per_step = [20, 10, 10, 5]
minfun = number_of_datasets * 3
#maxfun_per_step = [50, 75, 200]
maxfun_per_step = [minfun*1.5, minfun*2, minfun*5, minfun*10]
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
tx_shear_list = list()
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
    tx_shear_list.append(numpy.zeros(6))
    #tx = numpy.random.random(6)*(range_shear[1] - range_shear[0]) + range_shear[0]
    #tx = tx * numpy.pi / 180.0
    #tx_shear_list.append(tx)
    #sxy = tx[0]
    #sxz = tx[1]
    #syx = tx[2]
    #syz = tx[3]
    #szx = tx[4]
    #szy = tx[5]
    #skewx = vtk.vtkMatrix4x4()
    #skewy = vtk.vtkMatrix4x4()
    #skewz = vtk.vtkMatrix4x4()
    #skewx.SetElement(2, 1, numpy.tan(szy))
    #skewx.SetElement(1, 2, numpy.tan(syz))
    #skewy.SetElement(2, 0, numpy.tan(szx))
    #skewy.SetElement(0, 2, numpy.tan(sxz))
    #skewz.SetElement(1, 0, numpy.tan(sxy))
    #skewz.SetElement(0, 1, numpy.tan(syx))
    #vtktrans.Concatenate(skewx)
    #vtktrans.Concatenate(skewy)
    #vtktrans.Concatenate(skewz)   
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
# TO USE NEW OBJECTIVE FUNCTION
#register.use_probabilistic_objective()
    
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
    number_of_fibers_step_four = number_of_fibers_per_step[3]

    sigma_step_one = sigma_per_step[0]
    sigma_step_two = sigma_per_step[1]
    sigma_step_three = sigma_per_step[2]
    sigma_step_four = sigma_per_step[3]

    maxfun_step_one = maxfun_per_step[0]
    maxfun_step_two = maxfun_per_step[1]
    maxfun_step_three = maxfun_per_step[2]
    maxfun_step_four = maxfun_per_step[3]

    # STEP ONE
    # run the basic iteration of translate, rotate, scale, twice
    start = time.time()
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
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    # Don't scale at the first step
    #register.scale_only()
    #register.compute()

    elapsed1 = (time.time() - start)

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

    # STEP TWO: add scale (and shear if present)
    # reset the time, we don't need to include rendering and saving images here
    start = time.time()
    # run the basic iteration of translate, rotate, scale AGAIN
    register.fiber_sample_size = number_of_fibers_step_two
    register.sigma = sigma_step_two
    register.maxfun = maxfun_step_two
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.scale_only()
    register.compute()
    #register.shear_only()
    #register.compute()
    # run again
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.scale_only()
    register.compute()
    #register.shear_only()
    #register.compute()
    
    elapsed2 = (time.time() - start)

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
    # reset the time, we don't need to include rendering and saving images here
    start = time.time()
    # run the basic iteration of translate, rotate, scale AGAIN
    register.fiber_sample_size = number_of_fibers_step_three
    register.sigma = sigma_step_three
    register.maxfun = maxfun_step_three
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.scale_only()
    register.compute()
    #register.shear_only()
    #register.compute()
    
    elapsed3 = (time.time() - start)

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

# Save the transforms calculated from the 3 big iterations
calculated_transforms = list()
for subj in register._subjects:
    calculated_transforms.append(subj.transform)

gold_standard = list()
for idx in range(number_of_datasets):
    rand_trans = list(tx_rot_list[idx]*numpy.pi/180.0) + \
                 list(tx_trans_list[idx]) + \
                 list(tx_scale_list[idx]) + \
                 list(tx_shear_list[idx])
    gold_standard.append(rand_trans)

calculated_transforms = numpy.array(calculated_transforms)
gold_standard = numpy.array(gold_standard)

# remove the mean transform from the gold standard
meantrans = numpy.mean(gold_standard,0)
gold_standard[:,0:6] -= meantrans[0:6]
gold_standard[:,6:9] = numpy.divide(gold_standard[:,6:9], meantrans[6:9])

#for idx in range(number_of_datasets):
#    print gold_standard[idx,0:3] + calculated_transforms[idx,0:3]
# want calculated transform to be INVERSE of gold standard (applied transform)
# that is why we do gs - - calculated
# and
# 1/gs - calculated
error = gold_standard + calculated_transforms
#error[:,6:9] = numpy.multiply(gold_standard[:,6:9], calculated_transforms[:,6:9])
error[:,6:9] = numpy.divide(1.0, gold_standard[:,6:9]) - calculated_transforms[:,6:9]

#print numpy.round((error*100))/100.0

# STEP FOUR
# this is a very brief iteration at 5mm...
# reset the time, we don't need to include rendering and saving images here
start = time.time()
# run the basic iteration of translate, rotate, scale AGAIN
# above errors are pretty low. see if we can do better,
# is the information there?
register.fiber_sample_size = number_of_fibers_step_four
register.sigma = sigma_step_four
register.maxfun = maxfun_step_four
register.translate_only()
register.compute()
register.rotate_only()
register.compute()
register.scale_only()
register.compute()
#register.shear_only()
#register.compute()
# run three times like previous step (2+3) at 10mm
register.translate_only()
register.compute()
register.rotate_only()
register.compute()
register.scale_only()
register.compute()
register.translate_only()
register.compute()
register.rotate_only()
register.compute()
register.scale_only()
register.compute()
elapsed4 = (time.time() - start)

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

print 'ERROR 10 mm---------'
print numpy.sqrt(numpy.mean(numpy.multiply(error,error),0))
print 'ERROR 5 mm---------'
print numpy.sqrt(numpy.mean(numpy.multiply(error2,error2),0))

#print "TIME ELAPSED 1 2 3 4:", elapsed1, elapsed2, elapsed3, elapsed4

print "TIME ELAPSED in minutes (iteration 1 2 3 4):"
print elapsed1 / 60.0, elapsed2 / 60.0, elapsed3 / 60.0, elapsed4 / 60.0

print "SIGMAS per step:"
print sigma_per_step

conversions = [180/numpy.pi, 180/numpy.pi, 180/numpy.pi, 1, 1, 1, 1, 1,1]

msq = numpy.sqrt(numpy.mean(numpy.multiply(error,error),0))
print msq[0:9]*conversions

msq3 = numpy.sqrt(numpy.mean(numpy.multiply(error2,error2),0))
print msq3[0:9]*conversions

# absolute errors:
abs_error = numpy.abs(error)
abs_error2 = numpy.abs(error2)

# mean and std for each component
mean_abs_error = numpy.mean(abs_error, 0)[0:9]*conversions
std_abs_error = numpy.std(abs_error,0)[0:9]*conversions

mean_abs_error2 = numpy.mean(abs_error2, 0)[0:9]*conversions
std_abs_error2 = numpy.std(abs_error2,0)[0:9]*conversions

print "ERRORS at 10mm:"
print mean_abs_error
print "+/-"
print std_abs_error

print "ERRORS at 5mm:"
print mean_abs_error2
print "+/-"
print std_abs_error2


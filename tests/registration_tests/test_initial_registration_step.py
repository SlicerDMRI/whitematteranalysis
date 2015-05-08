import matplotlib.pyplot as plt
import numpy
import os

import vtk

import whitematteranalysis as wma

#fname = '/Users/odonnell/Dropbox/Data/TBI_FE_PNL/controls/01231-dwi-filt-Ed-DTI-tract.vtp'
#pd_fixed = wma.filter.downsample(wma.filter.preprocess(wma.io.read_polydata(fname), 30), 2000)
#pd_moving = wma.filter.downsample(wma.filter.preprocess(wma.io.read_polydata(fname), 30), 2000)
#wma.io.write_polydata(pd_fixed, 'test_objective_pd_fixed.vtp')
#wma.io.write_polydata(pd_moving, 'test_objective_pd_moving.vtp')

number_of_fibers = 500
test_maxfun = 5
test_objective_index = 4

#fname1 = 'test_data/wm_0001.vtp'
#fname2 = 'test_data/wm_0001.vtp'
#pd_fixed = wma.filter.downsample(wma.filter.preprocess(wma.io.read_polydata(fname1), 30, verbose=False), number_of_fibers, verbose=False)
#pd_moving = wma.filter.downsample(wma.filter.preprocess(wma.io.read_polydata(fname2), 30, verbose=False), number_of_fibers, verbose=False)

# same as in MICCAI paper.
pd_fixed = wma.io.read_polydata('test_data_objective/test_objective_pd_fixed.vtp')
pd_moving = wma.io.read_polydata('test_data_objective/test_objective_pd_moving.vtp')

transformer = vtk.vtkTransformPolyDataFilter()
transformer.SetInput(pd_moving)
vtktrans = vtk.vtkTransform()
transformer.SetTransform(vtktrans)

# translation
objective_list = list()
all_translations = list()
all_sigmas = list()
# 33 translations from -40 to 40
translation_list = numpy.array(range(-80,85,5))*0.5
sigma_list = [5, 10, 20, 30]
for x in translation_list:
    vtktrans.Identity()
    vtktrans.Translate(x, 0, 0)
    transformer.Update()
    for sigma in sigma_list:
        register = wma.register_two_subjects.RegisterTractography()
        
        fibers = wma.fibers.FiberArray()
        fibers.convert_from_polydata(pd_fixed, 5)
        fibers = numpy.array([fibers.fiber_array_r,fibers.fiber_array_a,fibers.fiber_array_s])
        register.fixed = fibers

        fibers = wma.fibers.FiberArray()
        fibers.convert_from_polydata(transformer.GetOutput(), 5)
        fibers = numpy.array([fibers.fiber_array_r,fibers.fiber_array_a,fibers.fiber_array_s])
        register.moving = fibers

        register.maxfun = test_maxfun
        register.sigma = sigma
        register.compute()
        objective_list.append(register.objective_function_values[test_objective_index])
        all_translations.append(x)
        all_sigmas.append(sigma)

translations = numpy.array(all_translations)
sigmas = numpy.array(all_sigmas)

objectives = numpy.array(objective_list)
sigma_5 = numpy.nonzero(sigmas == 5)[0]
sigma_10 = numpy.nonzero(sigmas == 10)[0]
sigma_20 = numpy.nonzero(sigmas == 20)[0]
sigma_30 = numpy.nonzero(sigmas == 30)[0]

plt.rcParams['font.size'] = 19

plt.figure() 

plt.plot(translations[sigma_5], objectives[sigma_5], linewidth=2)
plt.plot(translations[sigma_10], objectives[sigma_10], linewidth=2)
plt.plot(translations[sigma_20], objectives[sigma_20], linewidth=2)
plt.plot(translations[sigma_30], objectives[sigma_30], linewidth=2)

plt.savefig('objective_function_trans_x_prob_one_step.pdf')

plt.figure() 
plt.plot(translations[sigma_5], objectives[sigma_5],'k', linewidth=2)
plt.plot(translations[sigma_10], objectives[sigma_10],'k', linewidth=2)
plt.plot(translations[sigma_20], objectives[sigma_20],'k', linewidth=2)
plt.plot(translations[sigma_30], objectives[sigma_30],'k', linewidth=2)

plt.savefig('objective_function_trans_x_prob_bw_one_step.pdf')

plt.figure() 
plt.plot(translations[sigma_10], objectives[sigma_10], linewidth=2)
plt.plot(translations[sigma_20], objectives[sigma_20], linewidth=2)
plt.plot(translations[sigma_30], objectives[sigma_30], linewidth=2)

plt.savefig('objective_function_trans_x_prob_sig_10_30_one_step.pdf')

plt.figure() 
plt.plot(translations[sigma_10], objectives[sigma_10],'k', linewidth=2)
plt.plot(translations[sigma_20], objectives[sigma_20],'k', linewidth=2)
plt.plot(translations[sigma_30], objectives[sigma_30],'k', linewidth=2)

plt.savefig('objective_function_trans_x_prob_bw_sig_10_30_one_step.pdf')


# rotation
objective_list = list()
all_rotations = list()
all_sigmas = list()
# 33 rotations from -40 to 40
rotation_list = numpy.array(range(-80,85,5))*0.5
sigma_list = [5, 10, 20, 30]
for r in rotation_list:
    vtktrans.Identity()
    vtktrans.RotateX(r)
    transformer.Update()
    for sigma in sigma_list:
        register = wma.register_two_subjects.RegisterTractography()
        
        fibers = wma.fibers.FiberArray()
        fibers.convert_from_polydata(pd_fixed, 5)
        fibers = numpy.array([fibers.fiber_array_r,fibers.fiber_array_a,fibers.fiber_array_s])
        register.fixed = fibers

        fibers = wma.fibers.FiberArray()
        fibers.convert_from_polydata(transformer.GetOutput(), 5)
        fibers = numpy.array([fibers.fiber_array_r,fibers.fiber_array_a,fibers.fiber_array_s])
        register.moving = fibers

        register.maxfun = test_maxfun
        register.sigma = sigma
        register.compute()
        objective_list.append(register.objective_function_values[test_objective_index])
        all_rotations.append(r)
        all_sigmas.append(sigma)

rotations = numpy.array(all_rotations)
sigmas = numpy.array(all_sigmas)

objectives = numpy.array(objective_list)
sigma_5 = numpy.nonzero(sigmas == 5)[0]
sigma_10 = numpy.nonzero(sigmas == 10)[0]
sigma_20 = numpy.nonzero(sigmas == 20)[0]
sigma_30 = numpy.nonzero(sigmas == 30)[0]

plt.rcParams['font.size'] = 19

plt.figure() 
plt.plot(rotations[sigma_5], objectives[sigma_5], linewidth=2)
plt.plot(rotations[sigma_10], objectives[sigma_10], linewidth=2)
plt.plot(rotations[sigma_20], objectives[sigma_20], linewidth=2)
plt.plot(rotations[sigma_30], objectives[sigma_30], linewidth=2)

plt.savefig('objective_function_rot_x_prob_one_step.pdf')

plt.figure() 
plt.plot(rotations[sigma_5], objectives[sigma_5],'k', linewidth=2)
plt.plot(rotations[sigma_10], objectives[sigma_10],'k', linewidth=2)
plt.plot(rotations[sigma_20], objectives[sigma_20],'k', linewidth=2)
plt.plot(rotations[sigma_30], objectives[sigma_30],'k', linewidth=2)

plt.savefig('objective_function_rot_x_prob_bw_one_step.pdf')

plt.figure() 
plt.plot(rotations[sigma_10], objectives[sigma_10], linewidth=2)
plt.plot(rotations[sigma_20], objectives[sigma_20], linewidth=2)
plt.plot(rotations[sigma_30], objectives[sigma_30], linewidth=2)

plt.savefig('objective_function_rot_x_prob_sig_10_30_one_step.pdf')

plt.figure() 
plt.plot(rotations[sigma_10], objectives[sigma_10],'k', linewidth=2)
plt.plot(rotations[sigma_20], objectives[sigma_20],'k', linewidth=2)
plt.plot(rotations[sigma_30], objectives[sigma_30],'k', linewidth=2)

plt.savefig('objective_function_rot_x_prob_bw_sig_10_30_one_step.pdf')

# scale
objective_list = list()
all_scales = list()
all_sigmas = list()
# 34 scale factors from .5 to 1.5
scale_list = numpy.array(range(0,100,3))*.01 + 0.5
sigma_list = [5, 10, 20, 30]
for x in scale_list:
    vtktrans.Identity()
    vtktrans.Scale(x, 1, 1)
    transformer.Update()
    for sigma in sigma_list:
        register = wma.register_two_subjects.RegisterTractography()
        
        fibers = wma.fibers.FiberArray()
        fibers.convert_from_polydata(pd_fixed, 5)
        fibers = numpy.array([fibers.fiber_array_r,fibers.fiber_array_a,fibers.fiber_array_s])
        register.fixed = fibers

        fibers = wma.fibers.FiberArray()
        fibers.convert_from_polydata(transformer.GetOutput(), 5)
        fibers = numpy.array([fibers.fiber_array_r,fibers.fiber_array_a,fibers.fiber_array_s])
        register.moving = fibers

        register.maxfun = test_maxfun
        register.sigma = sigma
        register.compute()
        objective_list.append(register.objective_function_values[test_objective_index])
        all_scales.append(x)
        all_sigmas.append(sigma)

scales = numpy.array(all_scales)
sigmas = numpy.array(all_sigmas)

objectives = numpy.array(objective_list)
sigma_5 = numpy.nonzero(sigmas == 5)[0]
sigma_10 = numpy.nonzero(sigmas == 10)[0]
sigma_20 = numpy.nonzero(sigmas == 20)[0]
sigma_30 = numpy.nonzero(sigmas == 30)[0]

plt.rcParams['font.size'] = 19

plt.figure() 
plt.plot(scales[sigma_5], objectives[sigma_5], linewidth=2)
plt.plot(scales[sigma_10], objectives[sigma_10], linewidth=2)
plt.plot(scales[sigma_20], objectives[sigma_20], linewidth=2)
plt.plot(scales[sigma_30], objectives[sigma_30], linewidth=2)
plt.axis('tight')
plt.savefig('objective_function_scale_x_prob_one_step.pdf')

plt.figure() 
plt.plot(scales[sigma_5], objectives[sigma_5],'k', linewidth=2)
plt.plot(scales[sigma_10], objectives[sigma_10],'k', linewidth=2)
plt.plot(scales[sigma_20], objectives[sigma_20],'k', linewidth=2)
plt.plot(scales[sigma_30], objectives[sigma_30],'k', linewidth=2)
plt.axis('tight')
plt.savefig('objective_function_scale_x_prob_bw_one_step.pdf')

plt.figure() 
plt.plot(scales[sigma_10], objectives[sigma_10], linewidth=2)
plt.plot(scales[sigma_20], objectives[sigma_20], linewidth=2)
plt.plot(scales[sigma_30], objectives[sigma_30], linewidth=2)
plt.axis('tight')
plt.savefig('objective_function_scale_x_prob_sig_10_30_one_step.pdf')

plt.figure() 
plt.plot(scales[sigma_10], objectives[sigma_10],'k', linewidth=2)
plt.plot(scales[sigma_20], objectives[sigma_20],'k', linewidth=2)
plt.plot(scales[sigma_30], objectives[sigma_30],'k', linewidth=2)
plt.axis('tight')
plt.savefig('objective_function_scale_x_prob_bw_sig_10_30_one_step.pdf')


# shear
objective_list = list()
all_rotations = list()
all_sigmas = list()
# 33 rotations from -40 to 40
rotation_list = numpy.array(range(-80,85,5))*0.5
sigma_list = [5, 10, 20, 30]
for r in rotation_list:
    vtktrans.Identity()
    skewx = vtk.vtkMatrix4x4()
    szy = syz = r * numpy.pi / 180.0
    skewx.SetElement(2, 1, numpy.tan(szy))
    skewx.SetElement(1, 2, numpy.tan(syz))
    vtktrans.Concatenate(skewx)
    del skewx
    transformer.Update()
    for sigma in sigma_list:
        register = wma.register_two_subjects.RegisterTractography()
        
        fibers = wma.fibers.FiberArray()
        fibers.convert_from_polydata(pd_fixed, 5)
        fibers = numpy.array([fibers.fiber_array_r,fibers.fiber_array_a,fibers.fiber_array_s])
        register.fixed = fibers

        fibers = wma.fibers.FiberArray()
        fibers.convert_from_polydata(transformer.GetOutput(), 5)
        fibers = numpy.array([fibers.fiber_array_r,fibers.fiber_array_a,fibers.fiber_array_s])
        register.moving = fibers

        register.maxfun = test_maxfun
        register.sigma = sigma
        register.compute()
        objective_list.append(register.objective_function_values[test_objective_index])
        all_rotations.append(r)
        all_sigmas.append(sigma)

rotations = numpy.array(all_rotations)
sigmas = numpy.array(all_sigmas)

objectives = numpy.array(objective_list)
sigma_5 = numpy.nonzero(sigmas == 5)[0]
sigma_10 = numpy.nonzero(sigmas == 10)[0]
sigma_20 = numpy.nonzero(sigmas == 20)[0]
sigma_30 = numpy.nonzero(sigmas == 30)[0]

plt.rcParams['font.size'] = 19

plt.figure() 
plt.plot(rotations[sigma_5], objectives[sigma_5], linewidth=2)
plt.plot(rotations[sigma_10], objectives[sigma_10], linewidth=2)
plt.plot(rotations[sigma_20], objectives[sigma_20], linewidth=2)
plt.plot(rotations[sigma_30], objectives[sigma_30], linewidth=2)

plt.savefig('objective_function_shear_x_prob_one_step.pdf')

plt.figure() 
plt.plot(rotations[sigma_5], objectives[sigma_5],'k', linewidth=2)
plt.plot(rotations[sigma_10], objectives[sigma_10],'k', linewidth=2)
plt.plot(rotations[sigma_20], objectives[sigma_20],'k', linewidth=2)
plt.plot(rotations[sigma_30], objectives[sigma_30],'k', linewidth=2)

plt.savefig('objective_function_shear_x_prob_bw_one_step.pdf')

plt.figure() 
plt.plot(rotations[sigma_10], objectives[sigma_10], linewidth=2)
plt.plot(rotations[sigma_20], objectives[sigma_20], linewidth=2)
plt.plot(rotations[sigma_30], objectives[sigma_30], linewidth=2)

plt.savefig('objective_function_shear_x_prob_sig_10_30_one_step.pdf')

plt.figure() 
plt.plot(rotations[sigma_10], objectives[sigma_10],'k', linewidth=2)
plt.plot(rotations[sigma_20], objectives[sigma_20],'k', linewidth=2)
plt.plot(rotations[sigma_30], objectives[sigma_30],'k', linewidth=2)

plt.savefig('objective_function_shear_x_prob_bw_sig_10_30_one_step.pdf')

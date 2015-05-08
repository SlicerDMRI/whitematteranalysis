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

#pd_fixed = wma.io.read_polydata('test_objective_pd_fixed.vtp')
#pd_moving = wma.io.read_polydata('test_objective_pd_moving.vtp')
# same as in MICCAI paper.
pd_fixed = wma.io.read_polydata('test_data_objective/test_objective_pd_fixed.vtp')
#pd_moving = wma.io.read_polydata('test_data_objective/test_objective_pd_moving.vtp')
pd_moving = wma.io.read_polydata('test_data_objective/test_objective_pd_fixed.vtp')

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
        register = wma.congeal.CongealTractography()
        register.parallel_jobs = 3
        register.fiber_sample_size = 2000
        register.add_subject(pd_fixed)
        register.add_subject(transformer.GetOutput())
        register.maxfun = 1
        register.sigma = sigma
        register.compute()
        objective_list.append(register.objective_function_values[0])
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

plt.savefig('objective_function_trans_x_prob.pdf')

plt.figure() 
plt.plot(translations[sigma_5], objectives[sigma_5],'k', linewidth=2)
plt.plot(translations[sigma_10], objectives[sigma_10],'k', linewidth=2)
plt.plot(translations[sigma_20], objectives[sigma_20],'k', linewidth=2)
plt.plot(translations[sigma_30], objectives[sigma_30],'k', linewidth=2)

plt.savefig('objective_function_trans_x_prob_bw.pdf')

plt.figure() 
plt.plot(translations[sigma_10], objectives[sigma_10], linewidth=2)
plt.plot(translations[sigma_20], objectives[sigma_20], linewidth=2)
plt.plot(translations[sigma_30], objectives[sigma_30], linewidth=2)

plt.savefig('objective_function_trans_x_prob_sig_10_30.pdf')

plt.figure() 
plt.plot(translations[sigma_10], objectives[sigma_10],'k', linewidth=2)
plt.plot(translations[sigma_20], objectives[sigma_20],'k', linewidth=2)
plt.plot(translations[sigma_30], objectives[sigma_30],'k', linewidth=2)

plt.savefig('objective_function_trans_x_prob_bw_sig_10_30.pdf')


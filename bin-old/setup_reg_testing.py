#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/ipython

import os
import glob
import matplotlib.pyplot as plt
import numpy

import vtk

import whitematteranalysis as wma


print 'Read and preprocess'

indir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/test_scale_tracts/GAUSSIAN/sigma20_weight_greater_2'
outdir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/test_scale_tracts/GAUSSIAN/test_reg_tues_jan31'

minimum_length = 60
number_of_fibers = 150
number_of_fibers_step_one = 50
number_of_fibers_step_two = 75
number_of_fibers_step_three = 100

import multiprocessing
number_of_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', number_of_jobs

inputMask = "{0}/*.vtp".format(indir)
inputPolyDatas = glob.glob(inputMask)
input_pds = list()

#inputPolyDatas = inputPolyDatas[0:5]
inputPolyDatas = inputPolyDatas[0:8]
print inputPolyDatas

for fname in inputPolyDatas:
    print fname
    pd = wma.io.read_polydata(fname)
    pd3 = wma.filter.downsample(pd, number_of_fibers)
    input_pds.append(pd3)
 
# view input data
ren = wma.registration_functions.view_polydatas(input_pds)

register = wma.congeal.CongealTractography()
register.parallel_jobs = number_of_jobs
register.threshold = 0
register.points_per_fiber = 5
register.fiber_sample_size = number_of_fibers_step_one

# inputs are fixed, moving
for pd in input_pds:
    register.add_subject(pd)

model_pds=list()
for subj in register._subjects:
    model_pds.append(subj._original_fibers.convert_to_polydata())

#ren = wma.registration_functions.view_polydatas(model_pds)

# RUN registration (initial)
#register.sigma = 10
register.sigma = 20
register.maxfun = 200
register.translate_only()
register.compute()

register.rotate_only()
register.compute()

register.scale_only()
register.compute()

# view output data
output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
ren = wma.registration_functions.view_polydatas(output_pds)

plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
plt.savefig('objective_function.pdf')

if 0:
    # narrower rho parameters, solution should move less
    inc_rot = (3.0 / 180.0) * numpy.pi
    inc_trans = 5.0
    inc_scale = .2 
    register.set_rhobeg(inc_rot, inc_trans, inc_scale)
    inc_rot = (1.0 / 180.0) * numpy.pi
    inc_trans = 2.0
    inc_scale = .01
    register.set_rhoend(inc_rot, inc_trans, inc_scale)

# RUN registration (second)
#register.translate_only()
#register.maxfun = 250
#register.compute()


import numpy

import vtk
from joblib import Parallel, delayed

import whitematteranalysis as wma

print 'Read and preprocess'

number_of_fibers = 500
import multiprocessing
number_of_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', number_of_jobs
#number_of_jobs = 12

pd = wma.io.read_polydata('/Users/lauren/Data/TBI/Tracts/01030-dwi-filt-Ed-DTI-tract.vtp')
pd2 = wma.filter.preprocess(pd,80)
pd2 = wma.filter.downsample(pd2,number_of_fibers)
pdA = wma.filter.mask(pd2,numpy.ones(number_of_fibers),numpy.ones(number_of_fibers))


test = wma.congeal.RegistrationInformation()
test.points_per_fiber = 5
test.fiber_sample_size = 50
test.initialize(pd2)
test.initialize_fiber_sample()

test.transform[1]=1
#test.transform[8]=10

test.apply_transform()
threshold = 5
sigmasq = 30*30

similarity = Parallel(
    n_jobs=number_of_jobs, verbose=1)(
        delayed(wma.similarity.total_similarity)(
            test._moving_fibers.get_fiber(lidx),
            test._moving_fibers,
            threshold,
            sigmasq)
            for lidx in test._moving_fiber_sample)

print similarity[0]


vtktrans = test.convert_transform_to_vtk()
transformer = vtk.vtkTransformPolyDataFilter()
transformer.SetInput(pd2)
transformer.SetTransform(vtktrans)
transformer.Update()
pd3 = transformer.GetOutput()
test2 = wma.congeal.RegistrationInformation()
test2.points_per_fiber = 5
test2.fiber_sample_size = 50
test2.initialize(pd3)
test2.initialize_fiber_sample()
test2._moving_fiber_sample = test._moving_fiber_sample
test2.apply_transform()
threshold = 5
sigmasq = 30*30

similarity2 = Parallel(
    n_jobs=number_of_jobs, verbose=1)(
        delayed(wma.similarity.total_similarity)(
            test2._moving_fibers.get_fiber(lidx),
            test2._moving_fibers,
            threshold,
            sigmasq)
            for lidx in test2._moving_fiber_sample)

print similarity2[0]


print vtktrans.GetMatrix()

print test._original_fibers.fiber_array_r[0][0], \
    test._original_fibers.fiber_array_a[0][0], \
    test._original_fibers.fiber_array_s[0][0]

print test._moving_fibers.fiber_array_r[0][0], \
    test._moving_fibers.fiber_array_a[0][0], \
    test._moving_fibers.fiber_array_s[0][0]


print test2._moving_fibers.fiber_array_r[0][0], \
    test2._moving_fibers.fiber_array_a[0][0], \
    test2._moving_fibers.fiber_array_s[0][0]

error_r = test2._moving_fibers.fiber_array_r - test._moving_fibers.fiber_array_r
error_a = test2._moving_fibers.fiber_array_a - test._moving_fibers.fiber_array_a
error_s = test2._moving_fibers.fiber_array_s - test._moving_fibers.fiber_array_s

print 'Error r', numpy.max(error_r)
print 'Error a', numpy.max(error_a)
print 'Error s', numpy.max(error_s)

# test transform 1
print test._original_fibers.fiber_array_r[0][0]*(-numpy.sin(test.transform[1]) ) +numpy.cos(test.transform[1])*test._original_fibers.fiber_array_s[0][0]

tmp = test._original_fibers.fiber_array_r*(-numpy.sin(test.transform[1]) ) +numpy.cos(test.transform[1])*test._original_fibers.fiber_array_s
print tmp[0][0]

#!/usr/bin/env python

import whitematteranalysis as wma
import vtk
import numpy
print 'Read and preprocess'

minimum_length = 60

number_of_fibers = 1000
#number_of_fibers = 2000
import multiprocessing
number_of_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', number_of_jobs
#number_of_jobs = 12

pd = wma.io.read_polydata('/Users/lauren/Data/TBI/Tracts/01030-dwi-filt-Ed-DTI-tract.vtp')
pd2 = wma.filter.preprocess(pd, minimum_length)
pd2 = wma.filter.downsample(pd2, number_of_fibers)
pdA = wma.filter.mask(pd2, numpy.ones(number_of_fibers), numpy.ones(number_of_fibers))

ren1 = wma.render.render(pdA, axes=True)


align = wma.midsagalign.MidsagittalAlignment()

align.parallel_jobs = number_of_jobs
#align.parallel_jobs = 1
align.threshold = 0
#align.sigma = 10
align.sigma = 5
align.sigma = 3
align.points_per_fiber = 5
#align.fiber_sample_size = 100
#align.fiber_sample_size = 500
align.fiber_sample_size = 300

wm_align, transform = align.compute(pdA)

#transform = align.convert_transform_to_vtk()

#transformer = vtk.vtkTransformPolyDataFilter()
#transformer.SetInput(pdA)
#transformer.SetTransform(transform)
#transformer.Update()
#pdB = transformer.GetOutput()

ren = wma.render.render(wm_align, axes=True)
ren.save_views()








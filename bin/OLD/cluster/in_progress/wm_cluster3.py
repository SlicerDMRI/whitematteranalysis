#!/usr/bin/env python

import whitematteranalysis as wma
import vtk
import numpy
print 'Read and preprocess'

number_of_fibers = 2000
import multiprocessing
number_of_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', number_of_jobs
#number_of_jobs = 12

pd = wma.io.read_polydata('/Users/lauren/Data/TBI/Tracts/01030-dwi-filt-Ed-DTI-tract.vtp')
pd2 = wma.filter.preprocess(pd,80)
pd2 = wma.filter.downsample(pd2,number_of_fibers)
pdA = wma.filter.mask(pd2,numpy.ones(number_of_fibers),numpy.ones(number_of_fibers))

pd = wma.io.read_polydata('/Users/lauren/Data/TBI/Tracts/01035-dwi-filt-Ed-DTI-tract.vtp')
pd2 = wma.filter.preprocess(pd,80)
pd2 = wma.filter.downsample(pd2,number_of_fibers)
pdB = wma.filter.mask(pd2,numpy.ones(number_of_fibers),numpy.multiply(numpy.ones(number_of_fibers),2))

appender = vtk.vtkAppendPolyData()
appender.AddInput(pdA)
appender.AddInput(pdB)
appender.Update()
pd3 = appender.GetOutput()
wma.render.render(pd3)

register = wma.register.RegisterTractography()
register.parallelJobs = number_of_jobs
register.threshold = 5
register.sigma = 10
#register.sigma = 5
register.fiber_sample_size = 100
# RUN registration (initial)
# inputs are fixed, moving
register.initialize(pdA,pdB)
register.translate_only()
#register.maxiter = 10
register.maxiter = 100
register.maxfun = 50
register.compute()

register.fiber_sample_size = 400
register.maxiter = 200
register.maxfun = 200
register.sigma = 5
register.compute()

#register.translate_and_rotate=1
register.rotate_only()
register.compute()

register.translate_and_rotate()
register.fiber_sample_size = 1000
register.maxiter = 300
register.maxfun = 500
register.compute()

register.translate_and_rotate_and_scale()
register.fiber_sample_size = 1000
register.maxiter = 300
register.maxfun = 1000
register.threshold = 3
register.sigma = 5
register.compute()

transform = register.convert_transform_to_vtk()
transformer = vtk.vtkTransformPolyDataFilter()
transformer.SetInput(pdB)
transformer.SetTransform(transform)
transformer.Update()
pdB_moved = transformer.GetOutput()

pdB_moved = wma.filter.mask(pdB_moved,numpy.ones(number_of_fibers),numpy.multiply(numpy.ones(number_of_fibers),3))
appender = vtk.vtkAppendPolyData()
appender.AddInput(pdA)
#appender.AddInput(pdB)
appender.AddInput(pdB_moved)
appender.Update()
pd3 = appender.GetOutput()
#ren = wma.render.render(pd3, opacity=0.3, depth_peeling=True)
ren = wma.render.render(pd3)
ren.save_views()

# now only cluster the registered one
del pd3
del appender
appender = vtk.vtkAppendPolyData()
appender.AddInput(pdA)
appender.AddInput(pdB_moved)
appender.Update()
pd3 = appender.GetOutput()

# run the clustering
#output_polydata_h, cluster_idx_h = wma.cluster.hierarchical(pd3)
output_polydata_h, cluster_idx_h = wma.cluster.hierarchical(pd3,300,2,0.95,number_of_jobs)
wma.render.render(output_polydata_h)
import matplotlib.pyplot as plt
print 'View results'
plt.figure()
not_used = plt.hist(cluster_idx_h, max(cluster_idx_h))
#ren = wma.render.render(output_polydata_h)
wma.cluster.view_cluster_number(output_polydata_h,1)



output_polydata_s, cluster_numbers_s, color, embed = wma.cluster.spectral(pd3,number_of_jobs=number_of_jobs)

# view the whole thing
ren = wma.render.render(output_polydata_s)

# view it cluster by cluster (run line 2 several times)
cidx = 1
#fiber_mask = cluster_numbers_s == cidx; pd4 = wma.filter.mask(output_polydata_s, fiber_mask, color); wma.render.render(pd4); cidx = cidx+1

#wma.cluster.view_cluster_number(output_polydata_s,1,cluster_numbers_s)

wma.cluster.view_cluster_number(output_polydata_s, cidx)
cidx = cidx+1; wma.cluster.view_cluster_number(output_polydata_s, cidx)
# evaluate

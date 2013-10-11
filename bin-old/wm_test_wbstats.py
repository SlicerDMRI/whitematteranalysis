#!/usr/bin/env python

import whitematteranalysis as wma
import vtk
import numpy
print 'Read and preprocess'

#number_of_fibers = 5000
#number_of_fibers = 3000
number_of_fibers = 8000
fiber_length  = 30
import multiprocessing
number_of_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', number_of_jobs
#number_of_jobs = 12

print 'reading'
#pd = wma.io.read_polydata('/Users/lauren/Data/1T-fb/1T-01045.vtk')
pd = wma.io.read_polydata('/Users/odonnell/Data/TBI/TenMatchedTracts/01030-dwi-filt-Ed-DTI-tract.vtp')
print 'processing'
pd2 = wma.filter.preprocess(pd,fiber_length)
print 'random downsampling'
pd2 = wma.filter.downsample(pd2,number_of_fibers)
print 'putting group number into data'
pdA = wma.filter.mask(pd2,numpy.ones(number_of_fibers),numpy.ones(number_of_fibers))

print 'reading'
#pd = wma.io.read_polydata('/Users/lauren/Data/2T-fb/2T-wb-01045_100_0030_015.vtk')
pd = wma.io.read_polydata('/Users/odonnell/Data/TBI/2T-Tracts/2T-01030-state.vtk')
print 'processing'
pd2 = wma.filter.preprocess(pd,fiber_length)
print 'random downsampling'
pd2 = wma.filter.downsample(pd2,number_of_fibers)
print 'putting group number into data'
pdB = wma.filter.mask(pd2,numpy.ones(number_of_fibers),numpy.multiply(numpy.ones(number_of_fibers),2))

appender = vtk.vtkAppendPolyData()
appender.AddInput(pdA)
appender.AddInput(pdB)
appender.Update()
pd3 = appender.GetOutput()
pd3.GetCellData().GetArray(0).SetName('GROUP')
wma.render.render(pd3)

input_polydata = pd3
threshold = 0
number_of_jobs = 8
distance_matrix = \
    wma.cluster._pairwise_distance_matrix(input_polydata, threshold,
                                          number_of_jobs)


# first thing: how many from each group are within n mm?
sample_threshold = 4.0
n_fibers = distance_matrix.shape[0]

distance_thresholded = distance_matrix <= sample_threshold

g_1 = distance_thresholded[:, 0:number_of_fibers]
g_2 = distance_thresholded[:, number_of_fibers:2*number_of_fibers]

avg_fibers = list()
fiber_array = wma.fibers.FiberArray()
fiber_array.points_per_fiber = 15
fiber_array.convert_from_polydata(input_polydata)

for f_idx in range(n_fibers):
    print f_idx
    # find which fibers made the cutoff
    indices = numpy.nonzero(distance_thresholded[f_idx])[0]
    # average these fibers
    f_avg = fiber_array.get_fiber(indices[0])
    f_count = 1
    for f_avg_idx in indices[1:]:
        f_avg += fiber_array.get_fiber(f_avg_idx)
        f_count += 1
    f_avg /= f_count
    # output the average in a new pd
    avg_fibers.append(f_avg)

stats_group_1 = list()
stats_group_2 = list()
for f_idx in range(n_fibers):
    stats_group_1.append(numpy.sum(g_1[f_idx]))
    stats_group_2.append(numpy.sum(g_2[f_idx]))

# make a pd out of it with those as masks
current_fiber_array = wma.fibers.FiberArray()    
current_fiber_array.number_of_fibers = len(avg_fibers)
current_fiber_array.points_per_fiber = fiber_array.points_per_fiber
dims = [current_fiber_array.number_of_fibers, current_fiber_array.points_per_fiber]
# fiber data
current_fiber_array.fiber_array_r = numpy.zeros(dims)
current_fiber_array.fiber_array_a = numpy.zeros(dims)
current_fiber_array.fiber_array_s = numpy.zeros(dims)

curr_fidx = 0
for curr_fib in avg_fibers:
    current_fiber_array.fiber_array_r[curr_fidx] = curr_fib.r
    current_fiber_array.fiber_array_a[curr_fidx] = curr_fib.a
    current_fiber_array.fiber_array_s[curr_fidx] = curr_fib.s
    curr_fidx += 1

# output as pd
outpd = current_fiber_array.convert_to_polydata()
out_g1 = vtk.vtkFloatArray()
out_g2 = vtk.vtkFloatArray()
out_percent_g2  = vtk.vtkFloatArray()
for f_idx in range(n_fibers):
    out_g1.InsertNextTuple1(stats_group_1[f_idx])
    out_g2.InsertNextTuple1(stats_group_2[f_idx])
    out_percent_g2.InsertNextTuple1(100.0 * stats_group_2[f_idx] / (stats_group_1[f_idx] + stats_group_2[f_idx]))

out_g1.SetName('GROUP_1')
out_g2.SetName('GROUP_2')
out_percent_g2.SetName('PERCENT_GROUP_2')
outpd.GetCellData().AddArray(out_g1)
outpd.GetCellData().AddArray(out_g2)
outpd.GetCellData().AddArray(out_percent_g2)
outpd.GetCellData().SetActiveScalars('PERCENT_GROUP_2')

ren = wma.render.render(outpd, scalar_bar=True)
ren.save_views()

existence = ((numpy.array(stats_group_1) > 0) * (numpy.array(stats_group_2) > 0)) 
#outpd_exist = wma.filter.mask(outpd,numpy.ones(existence.shape),existence)
# in both
#outpd_exist = wma.filter.mask(outpd,existence,numpy.array(stats_group_1))
# in only one
outpd_exist = wma.filter.mask(outpd, ~existence, numpy.array(stats_group_1)>0)
ren2 = wma.render.render(outpd_exist, scalar_bar=True)

# run the clustering
#output_polydata_h, cluster_idx_h = wma.cluster.hierarchical(pd3)
#output_polydata_h, cluster_idx_h = wma.cluster.hierarchical(pd3,300,2,0.95,number_of_jobs)
#wma.render.render(output_polydata_h)
#import matplotlib.pyplot as plt
#print 'View results'
#plt.figure()
#not_used = plt.hist(cluster_idx_h, max(cluster_idx_h))

#output_polydata_s, cluster_numbers_s, color, embed = wma.cluster.spectral(pd3,number_of_jobs=number_of_jobs)

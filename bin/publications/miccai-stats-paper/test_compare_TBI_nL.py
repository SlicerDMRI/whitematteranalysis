#!/usr/bin/env python

# get rid of 1250

import os
import glob
import matplotlib.pyplot as plt
import numpy
import scipy.stats

import vtk

import whitematteranalysis as wma

import multiprocessing
parallel_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', parallel_jobs
parallel_jobs = 15
print 'Using N jobs:', parallel_jobs

input_directory = '/Users/odonnell/Dropbox/Coding/OUTPUTS/MICCAI2012/test_TBI_reg_feb11b/iteration_4'
input_mask = "{0}/*.vtk".format(input_directory)
input_poly_datas = glob.glob(input_mask)
#number_of_fibers_per_subject = 500
number_of_fibers_per_subject = 1000
percent_of_fibers_per_subject = 10
percent_of_fibers_per_subject = 50
fiber_length  = 60

# the first half are controls and the second half are patients
group = numpy.zeros(len(input_poly_datas))
group[0:13] = 1

print input_poly_datas
idx = 0

input_pds = list()
input_fiber_numbers = list()

for fname in input_poly_datas:
    print 'reading:', fname
    pd = wma.io.read_polydata(fname)
    print 'processing'
    pd2 = wma.filter.preprocess(pd, fiber_length)
    print 'random downsampling'
    n_fib = pd2.GetNumberOfLines() * percent_of_fibers_per_subject /100
    pd3 = wma.filter.downsample(pd2, n_fib)
    input_fiber_numbers.append(n_fib)
    #pd3 = wma.filter.downsample(pd2, number_of_fibers_per_subject)
    print 'putting group number into data'
    pdA = wma.filter.mask(pd3,numpy.ones(n_fib),numpy.ones(n_fib)*group[idx])
    print 'appending to list'
    input_pds.append(pdA)
    idx += 1

print 'Appending inputs into one polydata'
appender = vtk.vtkAppendPolyData()
for pd in input_pds:
    appender.AddInput(pd)

appender.Update()
pd_all = appender.GetOutput()
pd_all.GetCellData().GetArray(0).SetName('GROUP')
#ren = wma.render.render(pd_all)

# pick a subset of fibers
fiber_sample_size = 1000
points_per_fiber = 5
total_number_of_fibers = pd_all.GetNumberOfLines()
fiber_sample = numpy.random.random_integers(
            0, total_number_of_fibers - 1,
            fiber_sample_size)

# convert to array representation
fiber_array = wma.fibers.FiberArray()
fiber_array.convert_from_polydata(pd_all, points_per_fiber)

avg_fibers = list()
avg_fiber_array = wma.fibers.FiberArray()
avg_fiber_array.points_per_fiber = 15
print "converting polydata to array representation..."
avg_fiber_array.convert_from_polydata(pd_all)
print "done converting polydata to array representation..."

# find the sample's distances to all other fibers
distances = numpy.zeros([fiber_sample_size, total_number_of_fibers])

for idx in range(fiber_sample_size):
    print idx, '/', fiber_sample_size
    fiber = fiber_array.get_fiber(fiber_sample[idx])
    distances[idx,:] = wma.similarity.fiber_distance(fiber, fiber_array, threshold=0, distance_method='Hausdorff')
    
# figure out indices for each neighborhood
#sample_threshold = 30.0
sample_threshold = 20.0
sample_threshold = 10.0
distances_thresholded = distances <= sample_threshold
#group1 = fiber_sample < total_number_of_fibers/2
#group2 = fiber_sample > total_number_of_fibers/2
#group1 = range(total_number_of_fibers/2)
#group2 = range(total_number_of_fibers/2, total_number_of_fibers)
input_fiber_numbers = numpy.array(input_fiber_numbers)
ng1 = numpy.sum(input_fiber_numbers[0:13])
ng2 = numpy.sum(input_fiber_numbers[13:])
group1 = range(ng1)
group2 = range(ng1, ng1+ng2)
print len(group1)
print len(group2)
score_group_1 = numpy.sum(distances_thresholded[:,group1],1)
score_group_2 = numpy.sum(distances_thresholded[:,group2],1)
print len(score_group_1)
print len(score_group_2)

# figure out fiber counts per subject
number_of_subjects = len(input_poly_datas)
#f_per_subject = total_number_of_fibers/number_of_subjects
#assert number_of_fibers_per_subject == f_per_subject
subject_counts = numpy.zeros([fiber_sample_size, number_of_subjects])
for s_idx in range(number_of_subjects):
    #f_idxs = range(s_idx*number_of_fibers_per_subject, (s_idx + 1)*number_of_fibers_per_subject)
    f_idxs = range(numpy.sum(input_fiber_numbers[0:s_idx]),numpy.sum(input_fiber_numbers[0:s_idx+1]))
    subject_counts[:,s_idx] = numpy.sum(distances_thresholded[:, f_idxs],1)

subject_counts_group_1 = subject_counts[:,range(number_of_subjects/2)]
subject_counts_group_2 = subject_counts[:,range(number_of_subjects/2, number_of_subjects)]

t_stats = numpy.zeros(fiber_sample_size)
p_values = numpy.zeros(fiber_sample_size)
for f_idx in range(fiber_sample_size):
    t_stats[f_idx], p_values[f_idx] = \
        scipy.stats.ttest_ind(subject_counts_group_1[f_idx], subject_counts_group_2[f_idx])


for f_idx in range(fiber_sample_size):
    if (f_idx % 100) == 0:
        print f_idx
    # find which fibers made the cutoff
    indices = numpy.nonzero(distances_thresholded[f_idx])[0]
    # average these fibers (guaranteed to be at least 1, self)
    f_avg = avg_fiber_array.get_fiber(indices[0])
    f_count = 1
    for f_avg_idx in indices[1:]:
        f_avg += avg_fiber_array.get_fiber(f_avg_idx)
        f_count += 1
    f_avg /= f_count
    # output the average in a new pd
    avg_fibers.append(f_avg)

# make a pd out of it with those as masks
current_fiber_array = wma.fibers.FiberArray()    
current_fiber_array.number_of_fibers = len(avg_fibers)
current_fiber_array.points_per_fiber = avg_fiber_array.points_per_fiber
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
out_t = vtk.vtkFloatArray()
out_p = vtk.vtkFloatArray()

for f_idx in range(fiber_sample_size):
    out_g1.InsertNextTuple1(score_group_1[f_idx])
    out_g2.InsertNextTuple1(score_group_2[f_idx])
    out_percent_g2.InsertNextTuple1(100.0 * score_group_2[f_idx] / (score_group_1[f_idx] + score_group_2[f_idx]))
    out_t.InsertNextTuple1(t_stats[f_idx])
    #out_p.InsertNextTuple1((p_values[f_idx]<0.05)*(p_values[f_idx]*100))
    out_p.InsertNextTuple1(p_values[f_idx])
    
out_g1.SetName('GROUP_1')
out_g2.SetName('GROUP_2')
out_percent_g2.SetName('PERCENT_GROUP_2')
out_t.SetName('t-statistic')
out_p.SetName('p_value')

outpd.GetCellData().AddArray(out_g1)
outpd.GetCellData().AddArray(out_g2)
outpd.GetCellData().AddArray(out_percent_g2)
outpd.GetCellData().AddArray(out_t)
outpd.GetCellData().AddArray(out_p)
outpd.GetCellData().SetActiveScalars('PERCENT_GROUP_2')
outpd.GetCellData().SetActiveScalars('p_value')
outpd.GetCellData().SetActiveScalars('GROUP_2')

ren = wma.render.render(outpd, scalar_bar=True)
ren.save_views()

existence = ((numpy.array(score_group_1) > 0) * (numpy.array(score_group_2) > 0)) 
#outpd_exist = wma.filter.mask(outpd,numpy.ones(existence.shape),existence)
# in both
#outpd_exist = wma.filter.mask(outpd,existence,numpy.array(score_group_1))
# in only one
outpd_exist = wma.filter.mask(outpd, ~existence, numpy.array(score_group_1)>0)
ren2 = wma.render.render(outpd_exist, scalar_bar=True, scalar_range=[0,1])

f = open('test.txt','w')
for lidx in range(len(p_values)):
    f.write(str(p_values[lidx]))
    f.write('\n')
f.close()

print 'load test.txt into matlab to test fdr for now'

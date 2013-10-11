import os
import glob
import matplotlib.pyplot as plt
import numpy
import scipy.stats

import vtk

import whitematteranalysis as wma

import multiprocessing


# create some polydata objects to view the results
def fiber_list_to_fiber_array(fiber_list):
    fiber_array = wma.fibers.FiberArray()    
    fiber_array.number_of_fibers = len(fiber_list)
    fiber_array.points_per_fiber = len(fiber_list[0].r)
    dims = [fiber_array.number_of_fibers, fiber_array.points_per_fiber]
    # fiber data
    fiber_array.fiber_array_r = numpy.zeros(dims)
    fiber_array.fiber_array_a = numpy.zeros(dims)
    fiber_array.fiber_array_s = numpy.zeros(dims)
    curr_fidx = 0
    for curr_fib in fiber_list:
        fiber_array.fiber_array_r[curr_fidx] = curr_fib.r
        fiber_array.fiber_array_a[curr_fidx] = curr_fib.a
        fiber_array.fiber_array_s[curr_fidx] = curr_fib.s
        curr_fidx += 1
    return fiber_array


def add_array_to_polydata(pd, array, array_name='Test', array_type='Cell'):
    out_array = vtk.vtkFloatArray()
    for idx in range(len(array)):
        out_array.InsertNextTuple1(array[idx])
    out_array.SetName(array_name)
    ret = pd.GetCellData().AddArray(out_array)
    print ret
    pd.GetCellData().SetActiveScalars(array_name)
    return(pd)


parallel_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', parallel_jobs
#parallel_jobs *= 3
#parallel_jobs = 101
#parallel_jobs = 15
parallel_jobs = 10
print 'Using N jobs:', parallel_jobs

group_indices = [1, 0, 1, 0, 0, 1, 1, 0]
# 1 T, 2 C, 3 T, 4 C, 5 C, 6 T, 7 T, 8 C

indir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/MICCAI2012/mice_with_scalars'


input_mask = "{0}/*.vtk".format(indir)
input_poly_datas = glob.glob(input_mask)

print input_poly_datas

input_pds = list()
input_pds_downsampled = list()

number_of_fibers_per_subject = 3000
number_of_fiber_centroids = 1000
number_of_subjects = len(input_poly_datas)
points_per_fiber = 30
neighborhood_threshold = 15

#input_poly_datas = input_poly_datas[0:1]
# read in ones with scalars already
for fname in input_poly_datas:
    print fname
    pd = wma.io.read_polydata(fname)
    input_pds.append(pd)

# grab scalars of interest
input_mean_mds_per_subject = list()
input_pds_downsampled = list()
for pd in input_pds:
    pd2, fiber_indices = wma.filter.downsample(pd, number_of_fibers_per_subject,return_indices=True)
    # get MD only at fibers of interest
    md = pd.GetCellData().GetArray('mean_MD')
    md_subj = list()
    for idx in fiber_indices:
        md_subj.append(md.GetTuple1(idx))
    #md_subj = numpy.array(md_subj)
    #md_subj = numpy.array(md_avg_list)[fiber_indices]
    input_mean_mds_per_subject.append(md_subj)    
    input_pds_downsampled.append(pd2)

# assign to flat lists with subject idx
subj_idx = 0
input_data_per_fiber = list()
input_subject_idx_per_fiber = list()
input_group_idx_per_fiber = list()
for md_subj in input_mean_mds_per_subject:
    input_data_per_fiber += md_subj
    for data_point in md_subj:
        input_subject_idx_per_fiber.append(subj_idx)
        input_group_idx_per_fiber.append(group_indices[subj_idx])
    subj_idx +=1
    
# convert to arrays for dist and averaging
# use entire appended polydata (perhaps in future compute per-subject)
print 'Appending inputs into one polydata'
appender = vtk.vtkAppendPolyData()
for pd in input_pds_downsampled:
    appender.AddInput(pd)

appender.Update()

# convert to array representation
print 'Converting fibers to array representation for dist and averaging'
fiber_array = wma.fibers.FiberArray()
fiber_array.convert_from_polydata(appender.GetOutput(), points_per_fiber)
print 'Done converting fibers to array representation for dist and averaging'

# try to do some statistics
# random sample of fibers for stats
total_number_of_fibers = number_of_fibers_per_subject*number_of_subjects
fiber_sample = numpy.random.permutation(total_number_of_fibers - 1)
fiber_sample = fiber_sample[0:number_of_fiber_centroids]

# compute dists
# find the sample's distances to all other fibers
distances = numpy.zeros([number_of_fiber_centroids, total_number_of_fibers])

for idx in range(number_of_fiber_centroids):
    print idx, '/', number_of_fiber_centroids
    fiber = fiber_array.get_fiber(fiber_sample[idx])
    distances[idx,:] = wma.similarity.fiber_distance(fiber, fiber_array, threshold=0, distance_method='Hausdorff')
    
# according to neighborhood definition:
# compute avg fibers and MD stats in neighborhoods
average_fibers_list = list()
statistic_list = list()
significance_list = list()
data_0_list = list()
data_1_list = list()

for idx in range(number_of_fiber_centroids):
    print idx, '/', number_of_fiber_centroids
    neighborhood_indices = numpy.nonzero(distances[idx,:] < neighborhood_threshold)[0]
    hood_count = len(neighborhood_indices)
    # find average fiber in neighborhood
    avg_fiber = fiber_array.get_fiber(neighborhood_indices[0])
    for hood_idx in neighborhood_indices[1:]:
        avg_fiber += fiber_array.get_fiber(hood_idx)
    avg_fiber /= hood_count
    average_fibers_list.append(avg_fiber)
    # compute statistic(s) in neighborhood
    data_list = list()
    group_list = list()
    for hood_idx in neighborhood_indices:
        data_list.append(input_data_per_fiber[hood_idx])
        group_list.append(input_group_idx_per_fiber[hood_idx])
    # statistic: MD in group 1 vs MD in group 2
    data_list = numpy.array(data_list)
    group_list = numpy.array(group_list)
    g0 = numpy.nonzero(group_list == 0)[0]
    g1 = numpy.nonzero(group_list == 1)[0]
    g0_data = data_list[g0]
    g1_data = data_list[g1]
    if len(g0_data) and len(g1_data):
        t, p = scipy.stats.ttest_ind(g0_data, g1_data)
        data_0_list.append(numpy.mean(numpy.array(g0_data)))
        data_1_list.append(numpy.mean(numpy.array(g1_data)))
    else:
        t = p = numpy.nan
        data_0_list.append(numpy.nan)
        data_1_list.append(numpy.nan)
    statistic_list.append(t)
    significance_list.append(p)


f = open('pvals.txt','w')
for lidx in range(len(significance_list)):
    f.write(str(significance_list[lidx]))
    f.write('\n')
f.close()

print 'load pvals.txt into matlab to test fdr for now'
# output as pd
outpd = fiber_list_to_fiber_array(average_fibers_list).convert_to_polydata()
outpd = add_array_to_polydata(outpd, significance_list)
outpd = add_array_to_polydata(outpd, data_0_list, array_name='Data0')
outpd = add_array_to_polydata(outpd, data_1_list, array_name='Data1')

mask = ~numpy.isnan(data_1_list)
mask_idx = numpy.nonzero(mask)[0]
mask_pd = wma.filter.mask(outpd, mask)
mask_pd = add_array_to_polydata(mask_pd, numpy.array(data_1_list)[mask_idx], array_name='Data1')
mask_pd = add_array_to_polydata(mask_pd, numpy.array(data_0_list)[mask_idx], array_name='Data0')
mask_pd = add_array_to_polydata(mask_pd, numpy.array(data_0_list)[mask_idx]-numpy.array(data_1_list)[mask_idx], array_name='Difference')
mask_pd = add_array_to_polydata(mask_pd, numpy.array(significance_list)[mask_idx], array_name='P')
mask_pd.GetCellData().SetActiveScalars('Data0')
mask_pd.GetCellData().SetActiveScalars('Data1')
mask_pd.GetCellData().SetActiveScalars('Difference')
mask_pd.GetCellData().SetActiveScalars('P')

ren = wma.render.render(outpd, scalar_bar=True, scalar_range=[0,1])

# from matlab fdr
thresh = 0.0083
# .10
thresh = 0.0199
ren2 = wma.render.render(mask_pd, scalar_bar=True, scalar_range=[0,thresh])

# nhood 15mm, 10% fdr
thresh = 0.0382

# nhood 20mm, 5% fdr
thresh = 0.0242

# nhood 15mm, % fdr
thresh = 0.0162

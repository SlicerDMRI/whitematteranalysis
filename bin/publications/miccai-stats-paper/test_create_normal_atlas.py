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


#indir = '/Users/odonnell/Dropbox/Work/Publications/Results/2012_MICCAI_reg_results/test_nL_reg_N26/iteration_4/'
indir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/MICCAI2012/controls_with_scalars'

input_mask = "{0}/*.vtk".format(indir)
input_poly_datas = glob.glob(input_mask)

#input_poly_datas = input_poly_datas[0:2]
print input_poly_datas

input_pds = list()
input_pds_downsampled = list()

#number_of_fibers_per_subject = 1000
#number_of_fiber_centroids = 400

number_of_fibers_per_subject = 3000
number_of_fiber_centroids = 1500
#number_of_fiber_centroids = 1000
# this is about 2.4 GB of memory for the distances...
#number_of_fibers_per_subject = 6000
#number_of_fiber_centroids = 2000
number_of_subjects = len(input_poly_datas)
points_per_fiber = 30
fiber_length = 30

# read in ones with scalars already
for fname in input_poly_datas:
    print fname
    pd = wma.io.read_polydata(fname)
    input_pds.append(pd)

# downsample pds
input_pds_downsampled = list()
downsample_indices = list()
for pd in input_pds:
    pd2, fiber_indices1 = wma.filter.preprocess(pd, fiber_length, return_indices=True)
    pd2, fiber_indices2 = wma.filter.downsample(pd2, number_of_fibers_per_subject,return_indices=True)
    fiber_indices = fiber_indices1[fiber_indices2]
    #pd2, fiber_indices = wma.filter.downsample(pd, number_of_fibers_per_subject,return_indices=True)
    downsample_indices.append(fiber_indices)
    input_pds_downsampled.append(pd2)
    
# grab scalars of interest
input_mean_perp_per_subject = list()
input_mean_para_per_subject = list()
input_mean_md_per_subject = list()
input_mean_fa_per_subject = list()
pidx = 0
for pd in input_pds:
    print pidx, '----------------------------------------------------------'
    mean_perp = pd.GetCellData().GetArray('mean_perpendicular_diffusivity')
    mean_para = pd.GetCellData().GetArray('mean_parallel_diffusivity')
    mean_md = pd.GetCellData().GetArray('mean_MD')
    mean_fa = pd.GetCellData().GetArray('mean_FA')
    mean_perp_subj = list()
    mean_para_subj = list()
    mean_md_subj = list()
    mean_fa_subj = list()
    fiber_indices = downsample_indices[pidx]
    for idx in fiber_indices:
        mean_perp_subj.append(mean_perp.GetTuple1(idx))
        mean_para_subj.append(mean_para.GetTuple1(idx))
        mean_md_subj.append(mean_md.GetTuple1(idx))
        mean_fa_subj.append(mean_fa.GetTuple1(idx))
    input_mean_perp_per_subject.append(mean_perp_subj)    
    input_mean_para_per_subject.append(mean_para_subj)    
    input_mean_md_per_subject.append(mean_md_subj)    
    input_mean_fa_per_subject.append(mean_fa_subj)    
    pidx += 1

# convert to arrays for dist and averaging
# use entire appended polydata (perhaps in future compute per-subject)
print 'Appending inputs into one polydata'
appender = vtk.vtkAppendPolyData()
for pd in input_pds_downsampled:
    appender.AddInput(pd)

appender.Update()
print 'Done appending inputs into one polydata'

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


# ------------------------------------------------------
# Can re-run the below after changing neighborhood_threshold
# Or after changing group membership
# All slow processing happens above this line
# ------------------------------------------------------
#neighborhood_threshold = 30.0
#neighborhood_threshold = 15.0
neighborhood_threshold = 20.0

data_to_analyze = input_mean_fa_per_subject
    
# assign to flat lists with subject idx
subj_idx = 0
input_data_per_fiber = list()
input_subject_idx_per_fiber = list()
for fa_subj in data_to_analyze:
    input_data_per_fiber += fa_subj
    for data_point in fa_subj:
        input_subject_idx_per_fiber.append(subj_idx)
    subj_idx +=1
    
# according to neighborhood definition:
# compute avg fibers and FA stats in neighborhoods
average_fibers_list = list()
average_data_list = list()
hood_density_list = list()

sigma = neighborhood_threshold/2.0
sigma_sq = sigma*sigma

# TO LOOK AT AVERAGE BRAIN, AND GROUP STATS
for idx in range(number_of_fiber_centroids):
    print idx, '/', number_of_fiber_centroids
    neighborhood_indices = numpy.nonzero(distances[idx,:] < neighborhood_threshold)[0]
    d_sq = numpy.multiply(distances[idx,neighborhood_indices],distances[idx,neighborhood_indices])
    neighborhood_weights = numpy.exp(numpy.divide(-d_sq, sigma_sq))
    
    hood_count = len(neighborhood_indices)
    # find average fiber in neighborhood
    avg_fiber = fiber_array.get_fiber(neighborhood_indices[0])
    avg_fiber *= neighborhood_weights[0]
    data_list = list()
    for idx in range(1, hood_count):
        curr_fiber = fiber_array.get_fiber(neighborhood_indices[idx])
        curr_fiber *= neighborhood_weights[idx]
        avg_fiber += curr_fiber
        data_list.append(input_data_per_fiber[neighborhood_indices[idx]] * neighborhood_weights[idx])
        
    #avg_fiber /= hood_count
    total_weight = numpy.sum(neighborhood_weights)
    avg_fiber /= total_weight
    data_list = numpy.array(data_list)
    avg_data = numpy.sum(data_list)
    avg_data /= total_weight
    average_fibers_list.append(avg_fiber)
    average_data_list.append(avg_data)
    hood_density_list.append(total_weight)
    
# output as pd
outpd = fiber_list_to_fiber_array(average_fibers_list).convert_to_polydata()
outpd = add_array_to_polydata(outpd, average_data_list, array_name='AverageData')
outpd = add_array_to_polydata(outpd, hood_density_list, array_name='FiberDensity')
outpd.GetCellData().SetActiveScalars('AverageData')
wma.io.write_polydata(outpd,'atlas_info.vtp')
# mean FA
ren = wma.render.render(outpd, scalar_bar=True, scalar_range=[0.35,0.65])
# max FA
#ren = wma.render.render(outpd, scalar_bar=True, scalar_range=[0.5,0.95])
# lambda parallel
#ren = wma.render.render(outpd, scalar_bar=True, scalar_range=[0.0003,0.0011])

mask = numpy.array(hood_density_list) > 3.0
mask_idx = numpy.nonzero(mask)[0]
mask_pd = wma.filter.mask(outpd, mask)
mask_pd = add_array_to_polydata(mask_pd, numpy.array(average_data_list)[mask_idx], array_name='AverageData')
mask_pd.GetCellData().SetActiveScalars('AverageData')
ren2 = wma.render.render(mask_pd, scalar_bar=True, scalar_range=[0.35,0.65])


import os
import glob
import matplotlib.pyplot as plt
import numpy
import scipy.stats

import vtk

import whitematteranalysis as wma

import multiprocessing

READ_DATA = 1
# intial experiment
#number_of_fibers_per_subject = 3000
#number_of_fiber_centroids = 1000
#number_of_fibers_per_subject = 5000
#number_of_fiber_centroids = 2000
# this is about 2.4 GB of memory for the distances...
#number_of_fibers_per_subject = 6000
#number_of_fiber_centroids = 2000

number_of_fibers_per_subject = 10000
number_of_fiber_centroids = 2000

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
print 'Using N jobs:', parallel_jobs

#execfile('/Users/odonnell/Dropbox/Coding/Python/WhiteMatterAnalysis/bin/miccai-stats-paper/test_compute_FA.py')

indir = '/Users/odonnell/Data/tbi_with_scalars'

input_mask = "{0}/*.vtk".format(indir)
input_poly_datas = glob.glob(input_mask)

print input_poly_datas

number_of_subjects = len(input_poly_datas)
points_per_fiber = 30
#points_per_fiber = 20

#input_poly_datas = input_poly_datas[0:1]

# this is to produce the files with scalars
#if 0:
#    for fname in input_poly_datas:
#        print fname
#        pd = wma.io.read_polydata(fname)
#        pd, fa_lines_list, fa_avg_list = compute_scalar_measures(pd)
#        fname2 =  'scalars_' + os.path.basename(fname)
#        wma.io.write_polydata(pd, fname2)

# read in ones with scalars already
if READ_DATA:
    input_pds = list()
    for fname in input_poly_datas:
        print fname
        pd = wma.io.read_polydata(fname)
        input_pds.append(pd)

# grab scalars of interest
input_mean_fas_per_subject = list()
input_pds_downsampled = list()
downsample_indices = list()
for pd in input_pds:
    pd2, fiber_indices = wma.filter.downsample(pd, number_of_fibers_per_subject,return_indices=True)
    # get FA only at fibers of interest
    fa = pd.GetCellData().GetArray('mean_FA')
    fa_subj = list()
    for idx in fiber_indices:
        fa_subj.append(fa.GetTuple1(idx))
    #fa_subj = numpy.array(fa_subj)
    input_mean_fas_per_subject.append(fa_subj)    
    input_pds_downsampled.append(pd2)
    downsample_indices.append(fiber_indices)

# convert to arrays for dist and averaging
# use entire appended polydata (perhaps in future compute per-subject)
print 'Appending inputs into one polydata'
appender = vtk.vtkAppendPolyData()
for pd in input_pds_downsampled:
    appender.AddInputData(pd)

appender.Update()
print 'Done appending inputs into one polydata'

# convert to array representation
print 'Converting fibers to array representation for dist and averaging'
fiber_array = wma.fibers.FiberArray()
fiber_array.convert_from_polydata(appender.GetOutput(), points_per_fiber)
print 'Done converting fibers to array representation for dist and averaging'

# try to do some statistics
# random sample of fibers for stats
# this sample must be from the "atlas" or the control data!
# this is the first half of the subjects.
total_number_of_fibers = number_of_fibers_per_subject*number_of_subjects
fiber_sample = numpy.random.permutation(total_number_of_fibers/2 - 1)
fiber_sample = fiber_sample[0:number_of_fiber_centroids]

# compute dists
# find the sample's distances to all other fibers
distances = numpy.zeros([number_of_fiber_centroids, total_number_of_fibers])

for idx in range(number_of_fiber_centroids):
    print idx, '/', number_of_fiber_centroids
    fiber = fiber_array.get_fiber(fiber_sample[idx])
    distances[idx,:] = wma.similarity.fiber_distance(fiber, fiber_array, threshold=0, distance_method='Hausdorff')

print 'Done computing distances'


import os
import glob
import multiprocessing

import numpy
import vtk
import whitematteranalysis as wma
from joblib import Parallel, delayed

print 'LAUREN can similarity avoid explicit matrix replication'
print 'LAUREN use hausdorff distance and percentage of subjects represented. Can be interpreted!'

print 'Read and preprocess'

minimum_length = 60

#outdir = os.path.join('.', 'test_outputs')
#os.mkdir(outdir)

number_of_fibers = 300
number_of_fibers = 1000

number_of_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', number_of_jobs
#number_of_jobs = 12
#number_of_jobs = 16

input_dir = '/Users/lauren/Dropbox/Coding/OUTPUTS/group_register_test_13s_brage_moredatathreads'
input_mask = "{0}/*.vtk".format(input_dir)
input_pd_fnames = glob.glob(input_mask)
input_pds = list()

#input_pd_fnames = input_pd_fnames[0:8]

print input_pd_fnames

for fname in input_pd_fnames:
    print fname
    pd = wma.io.read_polydata(fname)
    pd2 = wma.filter.preprocess(pd, minimum_length)
    pd3 = wma.filter.downsample(pd2, number_of_fibers)
    input_pds.append(pd3)

# test clustering the data
appender = vtk.vtkAppendPolyData()
for pd in input_pds:
    appender.AddInput(pd)

appender.Update()
pd_all_registered = appender.GetOutput()

# figure out which fibers are from which subject
subject_idx = list()
sidx = 0
for pd in input_pds:
    num_lines = pd.GetNumberOfLines()
    print sidx, num_lines
    for line in range(0, num_lines):
        subject_idx.append(sidx)
    sidx = sidx + 1

total_number_of_fibers = len(subject_idx)

# USE REAL DISTANCES
#distances = wma.cluster._pairwise_distance_matrix(pd_all_registered, 0, number_of_jobs)


print 'Convert to array'
fiber_array = wma.fibers.FiberArray()
threshold = 0
distance_method = 'Hausdorff'
#fiber_array.convert_from_polydata(pd_all_registered, points_per_fiber=5)
fiber_array.convert_from_polydata(pd_all_registered, points_per_fiber=7)

# pairwise distance matrix
all_fibers = range(0, fiber_array.number_of_fibers)

print 'Compute fiber distances'
distances = Parallel(n_jobs=number_of_jobs,
                     verbose=1)(
                 delayed(wma.similarity.fiber_distance)(
                     fiber_array.get_fiber(lidx),
                     fiber_array,
                     threshold, distance_method)
                     for lidx in all_fibers)

distances = numpy.array(distances)

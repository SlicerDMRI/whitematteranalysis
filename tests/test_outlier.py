from __future__ import print_function
import numpy
import pickle
import whitematteranalysis as wma

fname = 'test_data_small/brain_0001.vtk'
expected_results_file = 'test_outlier.pkl'

print("Reading input file:", fname)

pd = wma.io.read_polydata(fname)

#print pd
min_fiber_distance = 15.0
n_jobs = 0
distance_method ='Mean'
#distance_method ='Hausdorff'
keep, mask, reject = wma.filter.remove_outliers(pd, min_fiber_distance, n_jobs=n_jobs, distance_method = distance_method)

wma.io.write_polydata(keep, 'keep.vtp')
wma.io.write_polydata(reject, 'reject.vtp')


# to store any new output mask
#pickle.dump(mask, open(expected_results_file, 'wb'))

# Open saved output mask, check for differences
expected_mask = pickle.load(open(expected_results_file, 'rb'))

print("Testing for differences in saved mask output and current outlier mask:")
#print numpy.array(expected_mask) - numpy.array(mask)
print("Overall difference?:", numpy.max(numpy.array(expected_mask) - numpy.array(mask)))


print('Done')


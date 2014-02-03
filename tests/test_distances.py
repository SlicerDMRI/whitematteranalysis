import numpy
import pickle
import whitematteranalysis as wma

expected_results_file = 'test_distances.pkl'
test_methods = ['Mean', 'Hausdorff', 'MeanSquared', 'StrictSimilarity']
sigma = 30
sigmasq = sigma*sigma

# test data
fname = 'test_data_small/brain_0001.vtk'
print fname
pd = wma.io.read_polydata(fname)

# test indices
indices = [1, 100, 300]

# get one fiber
mask = numpy.zeros(pd.GetNumberOfLines())
#idx = numpy.random.randint(pd.GetNumberOfLines())
#print idx
idx = indices[0]
mask[idx] = 1
pd_test = wma.filter.mask(pd, mask)
fiber_array_self = wma.fibers.FiberArray()
fiber_array_self.convert_from_polydata(pd_test, points_per_fiber=20)

# get a different fiber
mask = numpy.zeros(pd.GetNumberOfLines())
#idx = numpy.random.randint(pd.GetNumberOfLines())
#print idx
idx = indices[1]
mask[idx] = 1
pd_test = wma.filter.mask(pd, mask)
fiber_array_other = wma.fibers.FiberArray()
fiber_array_other.convert_from_polydata(pd_test, points_per_fiber=20)


# get another different fiber
mask = numpy.zeros(pd.GetNumberOfLines())
idx = indices[1]
mask[idx] = 2
pd_test = wma.filter.mask(pd, mask)
fiber_array_other2 = wma.fibers.FiberArray()
fiber_array_other2.convert_from_polydata(pd_test, points_per_fiber=20)

# make sure distances are correct for fiber to itself.
out_dist = []
for distance_method in test_methods:
    print "====="
    print "Test method", distance_method
    # distance to self
    distance1 = wma.similarity.fiber_distance(fiber_array_self.get_fiber(0), fiber_array_self, threshold=0, distance_method=distance_method)[0]
    # distance to other
    distance2 = wma.similarity.fiber_distance(fiber_array_self.get_fiber(0), fiber_array_other, threshold=0, distance_method=distance_method)[0]
    # distance to other
    distance3 = wma.similarity.fiber_distance(fiber_array_self.get_fiber(0), fiber_array_other2, threshold=0, distance_method=distance_method)[0]

    out_dist.append(distance1)
    out_dist.append(distance2)
    out_dist.append(distance3)
    
    if distance_method == 'StrictSimilarity':
        print "Self similarity calculated as:", distance1
        print "Other similarity calculated as:", distance2
        print "Other similarity 2 calculated as:", distance3
    else:
        print "Self distance calculated as:", distance1, "similarity:", wma.similarity.distance_to_similarity(distance1, sigmasq)
        print "Other distance calculated as:", distance2, "similarity:", wma.similarity.distance_to_similarity(distance2, sigmasq)
        print "Other distance 2 calculated as:", distance3, "similarity:", wma.similarity.distance_to_similarity(distance3, sigmasq)

# to store any new output distances
####pickle.dump(out_dist, open(expected_results_file, 'wb'))

# Open saved output file, check for differences
expected_dist = pickle.load(open(expected_results_file, 'rb'))

print "Testing for differences in saved distance output and current:"
print numpy.array(expected_dist) - numpy.array(out_dist)

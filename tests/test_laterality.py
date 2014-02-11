import numpy
import pickle
import whitematteranalysis as wma

expected_results_file = 'test_laterality.pkl'

sigma = 30.0
threshold = 0.0
equal_fibers = False
points_per_fiber = 20

# parameter ranges to test
test_sigmas = [10.0, 30.0, 50.0]
test_thresholds = [0.0, 2.0]
test_equal_fibers = [True, False]
test_points_per_fiber = [5, 20]

# test data
fname = 'test_data_small/brain_0001.vtk'
print fname
pd = wma.io.read_polydata(fname)

pd_symm = wma.filter.symmetrize(pd)


laterality = wma.laterality.WhiteMatterLaterality()
laterality.sigma = sigma
laterality.threshold = threshold
laterality.equal_fiber_num = equal_fibers
laterality.points_per_fiber = points_per_fiber
                
print ""
print "QUANTITATIVE RESULTS (CODE CHANGE) TESTING"
print "------------------------------------------"
expected_results_file = 'test_laterality.pkl'
laterality_results = laterality.compute(pd)
# to store any new output distances
#####pickle.dump(laterality_results.laterality_index, open(expected_results_file, 'wb'))
# Open saved output file, check for differences
expected_laterality = pickle.load(open(expected_results_file, 'rb'))

print "Testing for differences in saved distance output and current:"
error = numpy.array(expected_laterality) - laterality_results.laterality_index
print "difference range (should be 0.0):", numpy.min(error), numpy.max(error)



print ""
print "SYNTHETIC DATA TESTING"
print "----------------------------"
laterality_results = laterality.compute(pd_symm)
# this must be all 0.0 since it is perfectly symmetric
print "These should be 0 for a perfectly symmetric brain:"
#print len(numpy.nonzero(laterality_results.laterality_index)[0])
print "LI range:", numpy.min(laterality_results.laterality_index), numpy.max(laterality_results.laterality_index)

# asymmetric test data
mask_right = laterality_results.hemisphere == 1
mask_left = laterality_results.hemisphere == -1
mask_hem = laterality_results.hemisphere != 0

pd_right = wma.filter.mask(pd_symm, mask_right)
pd_left = wma.filter.mask(pd_symm, mask_left)

laterality_results = laterality.compute(pd_right)
print "These should be all +1 for a brain with only right hemisphere fibers"
print "LI range:", numpy.min(laterality_results.laterality_index), numpy.max(laterality_results.laterality_index)

laterality_results = laterality.compute(pd_left)
print "These should be all -1 for a brain with only left hemisphere fibers"
print "LI range:", numpy.min(laterality_results.laterality_index), numpy.max(laterality_results.laterality_index)

# this only works if there are fibers in each hemisphere
equal_fibers = True
laterality.equal_fiber_num = equal_fibers
# random test data
for rand_idx in range(3):
    number_of_rejected_fibers = 10
    mask_rand = numpy.random.permutation(pd_symm.GetNumberOfLines()) > number_of_rejected_fibers - 1
    pd_rand = wma.filter.mask(pd_symm, mask_rand)

    hem = mask_hem[mask_rand]
    laterality_results = laterality.compute(pd_rand)
    print "These should be from -1 to +1 for a brain with", number_of_rejected_fibers, "randomly asymmetric fibers"
    print "LI range (mean, median):", numpy.min(laterality_results.laterality_index), numpy.max(laterality_results.laterality_index), "(", numpy.mean(laterality_results.laterality_index[hem]), ",", numpy.median(laterality_results.laterality_index[hem]), ")"


print ""
print "REAL DATA PARAMETER TESTING"
print "----------------------------"
print "Testing initial input data, not synthetically symmetrized or asymetrized."
for equal_fibers in test_equal_fibers:
    for threshold in test_thresholds:
        for sigma in test_sigmas:
            for points_per_fiber in test_points_per_fiber:
                laterality.sigma = sigma
                laterality.threshold = threshold
                laterality.equal_fiber_num = equal_fibers
                laterality.points_per_fiber = points_per_fiber
                laterality_results = laterality.compute(pd)
                mask_hem = laterality_results.hemisphere != 0
                print "[sigma/threshold/equal-hem/ppf]", "[",sigma,"/",threshold,"/",equal_fibers,",", points_per_fiber, "]:", "(",numpy.min(laterality_results.laterality_index), ",", numpy.max(laterality_results.laterality_index), ")", "(", numpy.mean(laterality_results.laterality_index[mask_hem]), ",", numpy.median(laterality_results.laterality_index[mask_hem]), ")"





from __future__ import print_function
import pickle
import numpy
import os

print("First run test_wm_cluster_atlas.sh to produce output files for testing.")

outdir = 'output/test_cluster_atlas'
test_files = ['test_cluster_atlas_numbers.pkl',
              'test_cluster_atlas_colors.pkl',
              'test_cluster_atlas_embeddings.pkl']

for file_name in test_files:
    
    expected_results_file = file_name
    test_results_file = os.path.join(outdir, file_name)

    expected_value = pickle.load(open(expected_results_file, 'rb'))
    test_value = pickle.load(open(test_results_file, 'rb'))

    print("Testing for differences in", test_results_file, "vs", expected_results_file, ":")
    testval = numpy.array(expected_value) - numpy.array(test_value)
    print(testval)
    print(numpy.mean(numpy.abs(testval)))
    print(numpy.min(testval))
    print(numpy.max(testval))
    




import whitematteranalysis as wma
import numpy
import scipy.stats
#---------------------------
# Usage of tract_measurement
#---------------------------

print '====================='
print 'test 1: Column Comma'
measurement_file = "./test_data_measurement/column_comma.txt"
hierarchy = "Column"
separator = "Comma"
tm = wma.tract_measurement.load_measurement(measurement_file, hierarchy, separator)

print ' <tm> Extracted measurements:\n', tm.measurement_header
print ' <tm> Number of clusters:', tm.cluster_number
print ' <tm> Cluster ploydata path:', tm.cluster_path
print ' <tm> Measurement matrix: shape', tm.measurement_matrix.shape
print '      Each row gives the different measurements from each cluster.'
print '      Each column gives certain measurement from all clusters.'
print '      Examples:'
print '      # Measurements from Cluster 0', tm.cluster_path[0], 'are (first 4):', tm.measurement_matrix[0, 0:4]
print '      # Measurements of', tm.measurement_header[0], 'from all clusters are', tm.get_measurements_by_index(0)
print '      # Mean value of', tm.measurement_header[2], 'of all clusters is', numpy.mean(tm.get_measurements_by_index(2))


##
print '====================='
print 'test 2: Column Tab'
measurement_file = "./test_data_measurement/column_tab.txt"
hierarchy = "Column"
separator = "Tab"
tm = wma.tract_measurement.load_measurement(measurement_file, hierarchy, separator)

print ' <tm> Extracted measurements:\n', tm.measurement_header
print ' <tm> Number of clusters:', tm.cluster_number
print ' <tm> Cluster ploydata path:', tm.cluster_path
print ' <tm> Measurement matrix: shape', tm.measurement_matrix.shape
print '      Each row gives the different measurements from each cluster.'
print '      Each column gives certain measurement from all clusters.'
print '      Examples:'
print '      # Measurements from Cluster 0', tm.cluster_path[0], 'are (first 4)', tm.measurement_matrix[0, 0:4]
print '      # Measurements of', tm.measurement_header[0], 'from all clusters are', tm.get_measurements_by_index(0)
print '      # Mean value of', tm.measurement_header[2], 'of all clusters is', numpy.mean(tm.get_measurements_by_index(2))


##
print '====================='
print 'test 3: multi subjects'
measurement_folder = '/Users/fan/Desktop/tmp' # TODO: change data path, containing multiple measurement files
hierarchy = "Column"
separator = "Comma"

measurement_list = wma.tract_measurement.load_measurement_in_folder(measurement_folder, hierarchy, separator)

print ' <tm> Extract FA distribution of cluster 0 across all subjects'
vec_FA = []
for tm in measurement_list:
    vec_FA.append(tm.get_measurements_by_name('tensors.FractionalAnisotropyMean')[0]) # TODO: change the name according to 'measurement_header'

print ' <tm> FA vector length:', len(vec_FA)
print ' <tm> FA vector (first 5):', vec_FA[0:5]
print ' <tm> Mean FA of cluster 0 from all subjects:', numpy.mean(vec_FA)
print ' <tm> Perform t-test between the first 10 subjects and the second 10 subjects'
t, p = scipy.stats.ttest_ind(vec_FA[0:9], vec_FA[10:19], equal_var=False) # TODO: provide indices of different groups (e.g., patient V.S. control)
print ' <tm> t =', t, ', p =', p
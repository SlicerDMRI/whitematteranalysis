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
print '      # Measurements from Cluster 0', tm.cluster_path[0], 'are (first 4):', tm.measurement_matrix[0, 0:3]
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
print '      # Measurements from Cluster 0', tm.cluster_path[0], 'are (first 4)', tm.measurement_matrix[0, 0:3]
print '      # Measurements of', tm.measurement_header[0], 'from all clusters are', tm.get_measurements_by_index(0)
print '      # Mean value of', tm.measurement_header[2], 'of all clusters is', numpy.mean(tm.get_measurements_by_index(2))


##
print '====================='
print 'test 3: multi subjects'
measurement_folder = './test_data_measurement/group_test'
hierarchy = "Column"
separator = "Tab"

measurement_list = wma.tract_measurement.load_measurement_in_folder(measurement_folder, hierarchy, separator)

print ' <tm> Number of cases to be analyzed:', len(measurement_list)
idx = 0
for tm in measurement_list:
    print ' <tm> Case', idx, 'ID:', tm.case_id
    idx = idx + 1

print ' <tm> Extract FA distribution of cluster 0 across all subjects'
vec_FA = []
for tm in measurement_list:
    vec_FA.append(tm.get_measurements_by_name('tensor1.FractionalAnisotropy')[0])

print ' <tm> FA vector length:', len(vec_FA)
print ' <tm> FA vector:', vec_FA
print ' <tm> Mean FA of cluster 0 from all subjects:', numpy.mean(vec_FA)
print ' <tm> Perform t-test between the first 3 subjects and the second 3 subjects'
t, p = scipy.stats.ttest_ind(vec_FA[0:2], vec_FA[3:5], equal_var=False)
print ' <tm> t =', t, ', p =', p

##
print '====================='

print 'test 4: demographics'

dg = wma.tract_measurement.load_demographics('./test_data_measurement/demographics.xlsx')
print '<tm> Loading demographics file:', dg.demographics_file
print '<tm> Number of subject:', len(dg.case_id_list)
print '<tm> Number of fields:', len(dg.demographics_header)
print '<tm> Fields are: ', dg.demographics_header

# data access way 1
case_id_list = dg.case_id_list
group_id_list = dg.group_id_list
age_list = dg.get_demographics_by_header('Age')
age_list = map(float, age_list)
print '<tm> First 5 case ID:', case_id_list[:5]
print '<tm> First 5 case group:', group_id_list[0:5]
print '<tm> First 5 case age:', age_list[:5]

# data access way 2
case_id_list = dg.get_demographics_by_index(0)
group_id_list = dg.get_demographics_by_index(1)
age_list = dg.get_demographics_by_index(2)
age_list = map(float, age_list)
print '<tm>  Last 5 case ID:', case_id_list[-5:]
print '<tm>  Last 5 case group:', group_id_list[-5:]
print '<tm>  Last 5 case age:', age_list[-5:]
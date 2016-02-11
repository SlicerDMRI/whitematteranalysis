import whitematteranalysis as wma

##
print 'test 1: Column Comma'
measurement_file = "./test_data_measurement/column_comma.txt"
hierarchy = "Column"
separator = "Comma"
matrix = wma.tract_measurement.load_measurement(measurement_file, hierarchy, separator)
print matrix

##
print 'test 2: Column Tab'
measurement_file = "./test_data_measurement/column_tab.txt"
hierarchy = "Column"
separator = "Tab"
matrix = wma.tract_measurement.load_measurement(measurement_file, hierarchy, separator)
print matrix

##
print 'test 3: Column Space'
measurement_file = "./test_data_measurement/column_space.txt"
hierarchy = "Column"
separator = "Space"
matrix = wma.tract_measurement.load_measurement(measurement_file, hierarchy, separator)
print matrix

##
print 'test 4: multi subjects'
measurement_folder = '/Users/fan/Desktop/tmp' # change path for testing
hierarchy = "Column"
separator = "Comma"
measurement_list = wma.tract_measurement.load_measurement_in_folder(measurement_folder, hierarchy, separator)
print len(measurement_list)
print measurement_list[0]

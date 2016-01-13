import whitematteranalysis as wma

# measurement_file = "./test_data_measurement/column_comma.txt"
# hierarchy = "Column"
# separator = "Comma"
# matrix = wma.tract_measurement.load_measurement(measurement_file, hierarchy, separator)

# print 'c_c'
# print matrix

# measurement_file = "./test_data_measurement/column_tab.txt"
# hierarchy = "Column"
# separator = "Tab"
# matrix = wma.tract_measurement.load_measurement(measurement_file, hierarchy, separator)

# print 'c_t'
# print matrix

# measurement_file = "./test_data_measurement/column_space.txt"
# hierarchy = "Column"
# separator = "Space"
# matrix = wma.tract_measurement.load_measurement(measurement_file, hierarchy, separator)

# print 'c_s'
# print matrix

measurement_file = "./test_data_measurement/row_comma.txt"
hierarchy = "Row"
separator = "Comma"
matrix = wma.tract_measurement.load_measurement(measurement_file, hierarchy, separator)

print 'r_c'
print matrix
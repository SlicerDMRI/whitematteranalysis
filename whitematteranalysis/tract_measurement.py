import csv
import os
import glob
import numpy

class TractMeasurement:
    """Fiber tract scalar measurement obtained from Slicer module FiberTractScalarMeasurement."""

    def __init__(self):
        self.measurement_file = None
        self.cluster_number = None
        self.cluster_path = None
        self.measurement_header = None
        self.measurement_matrix = None
        self.hierarchy = 'Column'
        self.separator = 'Tab'


    def load(self):
        if not os.path.isfile(self.measurement_file):
            print "<tract_measurement.py> Error: Input file", self.measurement_file , "does not exist."
            raise AssertionError

        separator_list = ['Comma', 'Tab', 'Space'] 
        if not any(self.separator in s for s in separator_list):
            print "<tract_measurement.py> Error: Separator should be one of Comma, Tab or Space. "
            raise AssertionError
        if self.separator == 'Comma':
            separator_char = ','
        elif self.separator == 'Tab':
            separator_char = '\t'
        elif self.separator == 'Space':
            separator_char = ' '

        hierarchy_list = ['Row', 'Column'] 
        if not any(self.hierarchy in s for s in hierarchy_list):
            print "<tract_measurement.py> Error: Hierarchy should be one of Row or Column. "
            raise AssertionError

        txt_matrix = []
        with open(self.measurement_file, 'r') as txtfile:
            reader = csv.reader(txtfile, delimiter=separator_char, skipinitialspace=True,  quoting=csv.QUOTE_NONE)
            for row in reader:
                row = map(str.strip, row)
                txt_matrix.append(row)

        if self.hierarchy == 'Row':
            # TODO: transfer Row output into list
            print "<tract_measurement.py> Error: Only support Column currently. "
            raise AssertionError

        tmp_matrix = numpy.array(txt_matrix)

        self.cluster_path = tmp_matrix[1:, 0]
        self.measurement_header = tmp_matrix[0, 1:]
        self.measurement_matrix = tmp_matrix[1:, 1:].astype(numpy.float)
        self.cluster_number = self.measurement_matrix.shape[0]

    def check(self):
        # Simple check if the first three fields are Name, Num_Point and Num_Fiber
        header = self.measurement_header
        if header[0] != 'Num_Points' or header[1] != 'Num_Fibers':
            print "<tract_measurement.py> Error: Measurement loading failed. First three fields extracted are: \n 1. ", measures[0], "\n 2. ", measures[1], "\n 3. ", measures[2], "\nwhich are expected to be \n 1. Name \n 2. Num_Points \n 3. Num_Fibers. "
            raise AssertionError

    def get_measurements_by_name(self, query_header_name):
        if not numpy.sum(self.measurement_header == query_header_name) == 1:
            print " Error: Header", query_header_name, "cannot be found. Select from:\n", self.measurement_header
            return None

        header_index = numpy.where(self.measurement_header == query_header_name)[0][0]
        return self.measurement_matrix[:, header_index]

    def get_measurements_by_index(self, query_index):
        if query_index >= self.measurement_header.shape[0]:
            print " Error: Index", query_index, "should range from 0 to", self.measurement_header.shape[0]

        return self.measurement_matrix[:, query_index]


def load_measurement(measurement_file, hierarchy = 'Column', separator = 'Tab'):
    """ Load measurement for one subject
    """

    tm = TractMeasurement()
    tm.measurement_file = measurement_file
    tm.hierarchy = hierarchy
    tm.separator = separator  
    
    tm.load()
    tm.check()

    return tm

def load_measurement_in_folder(measurement_folder, hierarchy = 'Column', separator = 'Tab'):
    """ Load measurements for multiple subjects
    """

    # txt of csv files will be handled
    input_mask = "{0}/*.txt".format(measurement_folder)
    input_mask2 = "{0}/*.csv".format(measurement_folder)

    measurement_files = glob.glob(input_mask) + glob.glob(input_mask2)

    measurement_list = list()
    for m in measurement_files:
        measurement_list.append(load_measurement(m, hierarchy, separator))

    return measurement_list

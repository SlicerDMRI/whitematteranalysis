import csv
import os
import glob

class TractMeasurement:
    """Fiber tract scalar measurement obtained from Slicer module FiberTractScalarMeasurement."""

    def __init__(self):
        self.measurement_file = None
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

        self.measurement_matrix = [];
        with open(self.measurement_file, 'r') as txtfile:
            reader = csv.reader(txtfile, delimiter=separator_char, skipinitialspace=True,  quoting=csv.QUOTE_NONE)
            for row in reader:
                row = map(str.strip, row)
                self.measurement_matrix.append(row)

        if self.hierarchy == 'Row':
            # TODO: transfer Row output into list
            print "<tract_measurement.py> Error: Only support Column currently. "
            raise AssertionError

    def check(self):
        # Simple check if the first three fields are Name, Num_Point and Num_Fiber
        measures = self.measurement_matrix[0]
        if measures[0] != 'Name' or measures[1] != 'Num_Points' or measures[2] != 'Num_Fibers':
            print "<tract_measurement.py> Error: Measurement loading failed. First three fields extracted are: \n 1. ", measures[0], "\n 2. ", measures[1], "\n 3. ", measures[2], "\nwhich are expected to be \n 1. Name \n 2. Num_Points \n 3. Num_Fibers. "
            raise AssertionError

        
def load_measurement(measurement_file, hierarchy = 'Column', separator = 'Tab'):
    """ Load measurement for one subject
    """

    tm = TractMeasurement()
    tm.measurement_file = measurement_file
    tm.hierarchy = hierarchy
    tm.separator = separator  
    
    tm.load()
    tm.check()

    return tm.measurement_matrix

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

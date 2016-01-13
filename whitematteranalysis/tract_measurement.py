import csv
import os

class TractMeasurement:
    
    def __init__(self):
        self.measurement_file = None
        self.measurement_matrix = None
        self.hierarchy = 'Column'
        self.separator = ','

    def load(self):

        if not os.path.isfile(self.measurement_file):
            print "<tract_measurement.py> Error: Input file", self.measurement_file , "does not exist."
            raise AssertionError

        separator_list = ['Comma', 'Tab', 'Space'] 
        if not any(self.separator in s for s in separator_list):
            print "<tract_measurement.py> Error: Separator shold be one of Comma, Tab or Space. "
            raise AssertionError
        if self.separator == 'Comma':
            self.separator = ','
        elif self.separator == 'Tab':
            self.separator = '\t'
        elif self.separator == 'Space':
            self.separator = ' '

        hierarchy_list = ['Row', 'Column'] 
        if not any(self.hierarchy in s for s in hierarchy_list):
            print "<tract_measurement.py> Error: Hierarchy shold be one of Row or Column. "
            raise AssertionError

        self.measurement_matrix = [];
        with open(self.measurement_file, 'r') as txtfile:
            reader = csv.reader(txtfile, delimiter=self.separator, skipinitialspace=True,  quoting=csv.QUOTE_NONE)
            for row in reader:
                row = map(str.strip, row)
                self.measurement_matrix.append(row)

        if self.hierarchy == 'Row':
            self.measurement_matrix = map(list, zip(*self.measurement_matrix))

    def check(self):
        # Simple check if the first three fields are Name, Num_Point and Num_Fiber
        measures = self.measurement_matrix[0]
        if measures[0] != 'Name' or measures[1] != 'Num_Points' or measures[2] != 'Num_Fibers':
            print "<tract_measurement.py> Error: Measurement loading failed. First three fileds extracted are: \n 1. ", measures[0], "\n 2. ", measures[1], "\n 3. ", measures[2], "\nwhich are expected to be \n 1. Name \n 2. Num_Points \n 3. Num_Fibers. " 
            return False
        return True
        
def load_measurement(measurement_file, hierarchy = 'Column', separator = 'Comma'):
    tm = TractMeasurement()
    
    tm.measurement_file = measurement_file
    tm.hierarchy = hierarchy
    tm.separator = separator  
    tm.load()

    if not tm.check():
        raise AssertionError

    return tm.measurement_matrix




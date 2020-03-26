import csv
import os
import glob
import numpy
import xlrd

class TractMeasurement:
    """Fiber tract scalar measurement obtained from Slicer module FiberTractScalarMeasurement."""

    def __init__(self):
        self.case_id = None
        self.measurement_file = None
        self.cluster_number = None
        # For Fibers_File_Folder, cluster_path gives the cluster file (vtk/vtp) path
        # For Fibers_Hierarchy extraction, cluster_path gives the hierarchy name provided in the mrml file
        # TODO: This variable name is confusing if using Fibers_Hierarchy but still kept to avoid breaking wm_quality_control_cluster_measurements.py
        self.cluster_path = None
        self.measurement_header = None
        self.measurement_matrix = None
        self.hierarchy = 'Column'
        self.separator = 'Tab'


    def load(self):
        if not os.path.isfile(self.measurement_file):
            print("<tract_measurement.py> Error: Input file", self.measurement_file , "does not exist.")
            raise AssertionError

        separator_list = ['Comma', 'Tab', 'Space'] 
        if not any(self.separator in s for s in separator_list):
            print("<tract_measurement.py> Error: Separator should be one of Comma, Tab or Space. ")
            raise AssertionError
        if self.separator == 'Comma':
            separator_char = ','
        elif self.separator == 'Tab':
            separator_char = '\t'
        elif self.separator == 'Space':
            separator_char = ' '

        hierarchy_list = ['Row', 'Column'] 
        if not any(self.hierarchy in s for s in hierarchy_list):
            print("<tract_measurement.py> Error: Hierarchy should be one of Row or Column. ")
            raise AssertionError

        txt_matrix = []
        with open(self.measurement_file, 'r') as txtfile:
            reader = csv.reader(txtfile, delimiter=separator_char, skipinitialspace=True,  quoting=csv.QUOTE_NONE)
            for row in reader:
                row = list(map(str.strip, row))
                row = [('NAN' if (len(r)==0) else r) for r in row] # Replace '' to NAN
                txt_matrix.append(row)

        if self.hierarchy == 'Row':
            # TODO: transfer Row output into list
            print("<tract_measurement.py> Error: Only support Column currently. ")
            raise AssertionError

        tmp_matrix = numpy.array(txt_matrix)

        self.case_id = os.path.splitext(os.path.split(self.measurement_file)[1])[0]
        self.cluster_path = tmp_matrix[1:, 0]
        self.measurement_header = tmp_matrix[0, 1:]
        self.measurement_matrix = tmp_matrix[1:, 1:].astype(numpy.float)
        self.cluster_number = self.measurement_matrix.shape[0]

    def check(self):
        # Simple check if the first two fields are Num_Point and Num_Fiber
        header = self.measurement_header
        if not (header[0] == 'Num_Fibers' or header[1] == 'Num_Fibers'):
            print("<tract_measurement.py> Error: Measurement loading failed. One of three first fields should contain Num_Fibers. ")
            raise AssertionError

    def get_measurements_by_name(self, query_header_name):
        if not numpy.sum(self.measurement_header == query_header_name) == 1:
            print(" Error: Header", query_header_name, "cannot be found. Select from:\n", self.measurement_header)
            return None

        header_index = numpy.where(self.measurement_header == query_header_name)[0][0]
        return self.measurement_matrix[:, header_index]

    def get_measurements_by_index(self, query_index):
        if query_index >= self.measurement_header.shape[0]:
            print(" Error: Index", query_index, "should range from 0 to", self.measurement_header.shape[0])

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
    measurement_files = sorted(measurement_files)

    measurement_list = list()
    for m in measurement_files:
        measurement_list.append(load_measurement(m, hierarchy, separator))

    return measurement_list


class Demographics:
    """Group demographics."""

    def __init__(self):
        self.demographics_file = None
        self.demographics_header = None
        self.case_id_list = None
        self.group_id_list = None
        self.demographics = None

    def load(self):
        if not os.path.isfile(self.demographics_file):
            print("<tract_measurement.py> Error: Input file", self.demographics_file , "does not exist.")
            raise AssertionError

        if os.path.splitext(self.demographics_file)[1] != '.xls' and os.path.splitext(self.demographics_file)[1] != '.xlsx':
            print("<tract_measurement.py> Error: Either .xls or .xlsx file format is required.")
            raise AssertionError

        try:
            wb = xlrd.open_workbook(self.demographics_file)
            sh = wb.sheet_by_index(0)
            self.demographics_header = list(map(str, sh.row_values(0)))
        except:
            print("<tract_measurement.py> Error: Fail to load:", self.demographics_file)
            print("                       Please make sure the file has the demographics in the first sheet.")
            raise AssertionError

        self.demographics = list()
        for h in range(len(self.demographics_header)):
            col = list(map(str, sh.col_values(h)[1:]))
            self.demographics.append(col)

        self.case_id_list = self.demographics[0]
        self.group_id_list = self.demographics[1]

    def check(self):
        # Simple check if the first two fields are subjectID and groupID
        header = self.demographics_header
        if header[0] != 'subjectID' or header[1] != 'groupID':
            print("<tract_measurement.py> Error: Demographics loading failed. \n" \
                  "                       First two fields extracted are [", header[0], "] and [", header[1], "], which are expected to be [ subjectID ] and [ groupID ].")
            raise AssertionError

    def get_demographics_by_index(self, query_index):
        if query_index >= len(self.demographics_header):
            print(" Error: Index", query_index, "should range from 0 to", len(self.demographics_header))

        return self.demographics[query_index]

    def get_demographics_by_header(self, query_header_name):
        for h in range(len(self.demographics_header)):
            if query_header_name == self.demographics_header[h]:
                return self.demographics[h]

        print(" Error: Header", query_header_name, "cannot be found. Select from:\n", self.demographics_header)
        return None

def load_demographics(xlsx):
    """ Load load_demographics file
        Each row is one case, including subjectID, groupID, Age, etc.
    """

    dg = Demographics()
    dg.demographics_file = xlsx
    dg.load()
    dg.check()

    return dg

import csv
import os
import argparse

parser = argparse.ArgumentParser(
    description="Load the scalar measurement calculated by FiberTractMeasurements in 3D Slicer.",
    epilog="Written by Fan Zhang (fzhang@bwh.harvard.edu) and Lauren O\'Donnell. Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputMeasurementFile',
    help='A txt file obtained by FiberTractMeasurements.')
parser.add_argument(
    'outputDataFile',
    help='The output of the load measurement.')
parser.add_argument(
    '-separator', action="store", dest="separator",
    help='Comma (default), Tab or Space')
parser.add_argument(
    '-hierarchy', action="store", dest="hierarchy", 
    help='Row (default) or Column')

args = parser.parse_args()

if not os.path.isfile(args.inputMeasurementFile):
    print "<wm_load_measurement.py> Error: Input file", args.inputMeasurementFile, "does not exist."
    exit()

measure_file = args.inputMeasurementFile

out_data_file = args.outputDataFile
if os.path.exists(out_data_file):
    print "<wm_load_measurement.py> Output directory", out_data_file, "exists, deleting it."
    os.remove(out_data_file)

separator = ','
if args.separator is not None:
    separator_list = ['Comma', 'Tab', 'Space'] 
    if not any(args.separator in s for s in separator_list):
        print "<wm_load_measurement.py> Error: Separator shold be one of Comma, Tab or Space. "
        exit()
    if args.separator == 'Tab':
        separator = '\t'
    elif args.separator == 'Space':
        separator = ' '

hierarchy = 'Row'
if args.hierarchy is not None:
    hierarchy_list = ['Row', 'Column'] 
    if not any(args.hierarchy in s for s in hierarchy_list):
        print "<wm_load_measurement.py> Error: Hierarchy shold be one of Row or Column. "
        exit()
    hierarchy = args.hierarchy

print "<wm_load_measurement.py> Starting loading."
print ""
print "=====input measurements ======\n", measure_file
print "=====output Python data =====\n", out_data_file
print "=========================="
print ""

measure_matrix = []
with open(measure_file, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=separator, quoting=csv.QUOTE_NONE)
    for row in reader:
        row = map(str.strip, row)
        measure_matrix.append(row)

if hierarchy == 'Column':
    measure_matrix = map(list, zip(*measure_matrix))

measures = measure_matrix[0]
if measures[0] != 'Name' or measures[1] != 'Num_Points' or measures[2] != 'Num_Fibers':
    print "<wm_load_measurement.py> Error: Data loading failed. First three measures extracted are: \n 1. ", measures[0], "\n 2. ", measures[1], "\n 3. ", measures[2], "\nwhich should be \n 1. Name \n 2. Num_Points \n 3. Num_Fibers. " 
    exit()


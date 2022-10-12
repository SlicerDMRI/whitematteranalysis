#!/usr/bin/env python
import argparse
import os
import numpy

try:
    import whitematteranalysis as wma
except:
    print("<wm_label_from_atlas.py> Error importing white matter analysis package\n")
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Measures mean FA, etc. in tractography clusters in a directory.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputDirectory',
    help='A directory of tractography as vtkPolyData (.vtk or .vtp).')
parser.add_argument(
    'outputFile',
    help='Output measurements will be recorded here.')

args = parser.parse_args()

if not os.path.isdir(args.inputDirectory):
    print("Error: Input directory", args.inputDirectory, "does not exist.")
    exit()

print("<wm_measure_all_clusters.py> Starting computation.")
print("")
print("=====input directory ======\n", args.inputDirectory)
print("=====output file =====\n", args.outputFile)
print("==========================")
print("")


# =======================================================================
# Above this line is argument parsing. Below this line is the pipeline.
# =======================================================================

def compute_point_data_stats(pd, array_name):
    point_array = pd.GetPointData().GetArray(array_name)
    if point_array is None:
        return None
    # make sure this is a one-component scalar
    if point_array.GetNumberOfComponents() > 1:
        print("Error in compute_point_data_stats: Array", array_name, "has more than one component", point_array.GetNumberOfComponents(), ".")
        return None
    print(point_array)
    num_points = pd.GetNumberOfPoints()
    points_copy = numpy.zeros(num_points)
    
    for pidx in range(0, num_points):
        if (pidx % 1000) == 0:
            print("Point", pidx, '/', num_points)
        # this assumes we have scalars here
        points_copy[pidx] = point_array.GetTuple(pidx)[0]

    points_mean = numpy.mean(points_copy)
    #points_std = numpy.std(points_copy)
    #points_median = numpy.median(points_copy)

    print("Mean ", array_name, ":", points_mean)

    return points_mean



input_polydatas = wma.io.list_vtk_files(args.inputDirectory)
number_of_clusters = len(input_polydatas)

input_polydatas = input_polydatas[0:10]

print("<wm_measure_all_clusters.py> Input number of vtk/vtp files: ", number_of_clusters)

scalars = ['FA', 'Trace', 'FA1', 'FA2', 'Trace1', 'Trace2']

output_rows = list()

# read in data
input_pds = list()
for fname in input_polydatas:
    print(fname)
    # read data
    print("<wm_cluster_atlas.py> Reading input file:", fname)
    pd = wma.io.read_polydata(fname)
    # preprocessing step: minimum length
    print("<wm_cluster_atlas.py> Computing stats for input file:", fname)
    output_row = list()
    output_row.append(fname)
    for sc in scalars:
        stat = compute_point_data_stats(pd, sc)
        if stat is not None:
            output_row.append(sc)
            output_row.append(stat)
    output_rows.append(output_row)

# output a csv file
f = open(args.outputFile, 'w')

for row in output_rows:
    outstr = ''
    for item in row:
        if outstr != '':
            outstr = outstr + '\t'
        outstr = outstr + str(item)
    outstr = outstr + '\n'
    f.write(outstr)

f.close()

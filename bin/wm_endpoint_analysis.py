#!/usr/bin/env python
import argparse
import os

try:
    import whitematteranalysis as wma
except:
    print "<wm_measure_all_subjects> Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Endpoint analysis of fiber clusters ",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputDirectory',
    help='Directory of fiber clustering results obtained by <wm_cluster_from_altas.py> for one subject. ')
parser.add_argument(
    'inputLabelMap',
    help='Path of the label map file.')
parser.add_argument(
    'outputDirectory',
    help='Directory of output CSV files of fiber scalar measurement (computed using Slicer FiberEndPointFromLabelMap module).')
parser.add_argument(
    'modulePath',
    help='Path of the FiberEndPointFromLabelMap module.')

args = parser.parse_args()

if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

if not os.path.exists(args.inputLabelMap):
    print "Error: label map file", args.inputLabelMap, "does not exist."
    exit()

if not os.path.exists(args.modulePath):
    print "Error: FiberEndPointFromLabelMap", args.Slicer, "does not exist."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

inputdir = args.inputDirectory
if inputdir.endswith('/'):
    inputdir = inputdir[:-1]

sub_folder_name = os.path.split(inputdir)[1]

print "<wm_endpoint_analysis>. Starting processing."
print ""
print "=====input directory======\n", args.inputDirectory
print "=====input label map======\n", args.inputLabelMap
print "=====output directory=====\n", args.outputDirectory
print "=====module path====\n", args.modulePath

pds = wma.io.list_vtk_files(args.inputDirectory)

print "<wm_measure_all_subjects> found", len(pds), "vtk/vtp files."

os.system(args.modulePath + ' ' +  args.inputLabelMap + ' ' + args.inputDirectory + ' ' + \
          os.path.join(args.outputDirectory, sub_folder_name+'_endpoint.txt'))


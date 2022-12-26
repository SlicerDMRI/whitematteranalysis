#!/usr/bin/env python
import argparse
import os
import glob
from joblib import Parallel, delayed

try:
    import whitematteranalysis as wma
except:
    print("<wm_measure_endpoint_overlap> Error importing white matter analysis package\n")
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Measure overlaps of fiber clusters with cortical parcellation or fMRI functional areas. This is based on the 3D Slicer module FiberEndPointFromLabelMap.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputTractDirectory',
    help='Directory of fiber clustering results obtained by <wm_cluster_from_altas.py> of multiple subjects. Make sure only the fiber clustering results are stored in this folder, making one subdirectory corresponding to one subject.')
parser.add_argument(
    'inputLabelMapDirectory',
    help='Contains the parcellation or functional areas as label map files. Make sure that the input tract files and the label map files match each other in alphabetical order.')
parser.add_argument(
    'outputDirectory',
    help='Directory of output CSV files that shows the percentage of fiber clusters\' endpoint connecting to certain regions in the label map.')
parser.add_argument(
    'modulePath',
    help='Path of the 3D Slicer FiberEndPointFromLabelMap module.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')

args = parser.parse_args()

if not os.path.isdir(args.inputTractDirectory):
    print("Error: Input directory", args.inputTractDirectory, "does not exist.")
    exit()

if not os.path.isdir(args.inputLabelMapDirectory):
    print("Error: Input label map directory", args.inputLabelMapDirectory, "does not exist.")
    exit()

if not os.path.exists(args.modulePath):
    print("Error: FiberEndPointFromLabelMap", args.modulePath, "does not exist.")
    exit()

if args.numberOfJobs is not None:
    number_of_jobs = args.numberOfJobs
else:
    number_of_jobs = 1

if not os.path.exists(args.outputDirectory):
    print("Output directory", args.outputDirectory, "does not exist, creating it.")
    os.makedirs(args.outputDirectory)

print("<wm_endpoint_analysis>. Starting processing.")
print("")
print("=====input fiber cluster directory======\n", args.inputTractDirectory)
print("=====input label map directory======\n", args.inputLabelMapDirectory)
print("=====output directory=====\n", args.outputDirectory)
print("=====module path====\n", args.modulePath)
print('=====using N jobs:', number_of_jobs, "====\n")

tract_dir_list = os.listdir(args.inputTractDirectory)
tract_dir_list = sorted(tract_dir_list)

print("<wm_endpoint_analysis> found", len(tract_dir_list), "subjects.")

def list_label_map_files(input_dir):
    # Find input files
    input_mask = f"{input_dir}/*.nrrd"
    input_mask2 = f"{input_dir}/*.nhdr"
    input_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_fnames = sorted(input_fnames)
    return(input_fnames)

label_map_file_list = list_label_map_files(args.inputLabelMapDirectory)

print("<wm_endpoint_analysis> found", len(label_map_file_list), "label maps. \n")

if len(tract_dir_list) != len(label_map_file_list):
    print("Error: The number of subjects", len(tract_dir_list), "should be equal to the number of label maps", len(label_map_file_list))
    exit()

def extract_endpoint(tract_dir, lalel_map_file, args):

    pds = wma.io.list_vtk_files(os.path.join(args.inputTractDirectory, tract_dir))
    print("<wm_endpoint_analysis> Computing:", os.path.join(args.inputTractDirectory, tract_dir))
    print("                            using", lalel_map_file)
    print("                            with:", len(pds), "vtk/vtp files.")

    sub_name = os.path.split(tract_dir)[1]
    os.system(args.modulePath + ' ' +  lalel_map_file + ' ' + os.path.join(args.inputTractDirectory, tract_dir) + ' ' + \
              os.path.join(args.outputDirectory, sub_name+'_endpoint.txt')+ \
              ' > ' + os.path.join(args.outputDirectory, 'log'+sub_name))

Parallel(n_jobs=number_of_jobs, verbose=1)(
    delayed(extract_endpoint)(tract_dir, label_map_file, args)
    for tract_dir, label_map_file in zip(tract_dir_list, label_map_file_list))

def list_txt_files(input_dir):
    # Find input files
    input_mask = f"{input_dir}/*.txt"
    input_fnames = glob.glob(input_mask)
    input_fnames = sorted(input_fnames)
    
    return input_fnames

endpoint_txt_list = list_txt_files(args.outputDirectory)
print("<wm_endpoint_analysis> Endpoint analysis were measured for", len(endpoint_txt_list), "subjects.")

if len(tract_dir_list) != len(endpoint_txt_list):
    print("Error: The numbers of inputs and outputs are different. Check the log file of each subeject.")
else:
    os.system("rm -rf "+os.path.join(args.outputDirectory, 'log*'))

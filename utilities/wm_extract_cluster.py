#!/usr/bin/env python
print("wm_extract_cluster.py importing packages")

import glob
import argparse
import os
import shutil
import numpy
import whitematteranalysis as wma

parser = argparse.ArgumentParser(
    description="Grab one cluster from within all subject directories and rename to include subject ID. Output all clusters into the output directory",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'cluster', type=int, 
    help='The cluster number you would like to copy for all subjects.')
parser.add_argument(
    'inputDirectory',
    help='A directory containing subdirectories for all clustered subjects.')
parser.add_argument(
    'outputDirectory',
    help='The output directory will be created if it does not exist.')



args = parser.parse_args()


if not os.path.isdir(args.inputDirectory):
    print("<wm_cluster_atlas.py> Error: Input directory", args.inputDirectory, "does not exist or is not a directory.")
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print("<wm_cluster_atlas.py> Output directory", outdir, "does not exist, creating it.")
    os.makedirs(outdir)

# cluster filename we want
fname_c = f'cluster_{args.cluster:05d}.vtp'
# output MRML filename
fname_mrml = f'all_clusters_{args.cluster:05d}.mrml'
fname_mrml = os.path.join(outdir, fname_mrml)

# Find input directories that may contain clusters
input_mask = f"{args.inputDirectory}/*"
input_directories = sorted(glob.glob(input_mask))

# Loop over inputs and try to find clusters
fname_list = list()
for dir in input_directories:
    if os.path.isdir(dir):
        subject_id = os.path.basename(dir)
        print(dir)
        fname = os.path.join(dir, fname_c)
        if os.path.exists(fname):
            fname1 = subject_id+'_'+fname_c
            fname_list.append(fname1)
            fname2 = os.path.join(outdir, fname1)
            print(fname, "===>>>>", fname2)
            shutil.copy(fname, fname2)

# also output a MRML file to load them all into Slicer
# for now, random colors
colors = numpy.random.random((len(fname_list),3))*255.0
wma.mrml.write(fname_list, colors, fname_mrml, ratio = 1.0)



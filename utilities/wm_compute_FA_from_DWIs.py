#!/usr/bin/env python
import argparse
import os
import numpy
import glob
try:
    import whitematteranalysis as wma
except:
    print("<wm_label_from_atlas.py> Error importing white matter analysis package\n")
    raise

def list_nhdr_files(input_dir):
    # Find input files
    input_mask = f"{input_dir}/*.nhdr"
    input_mask2 = f"{input_dir}/*.nrrd"
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)

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
    'inputDirectoryDWI',
    help='A directory of DWIs (nhdr or nrrd).')
parser.add_argument(
    'inputDirectoryMask',
    help='A directory of masks (nhdr or nrrd).')
parser.add_argument(
    'outputDirectory',
    help='Output data will be saved here.')

args = parser.parse_args()

if not os.path.isdir(args.inputDirectoryDWI):
    print("Error: Input directory", args.inputDirectory, "does not exist.")
    exit()

if not os.path.isdir(args.inputDirectoryMask):
    print("Error: Input directory", args.inputDirectory, "does not exist.")
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print("<register> Output directory", outdir, "does not exist, creating it.")
    os.makedirs(outdir)

# get inputs
dwi_list = list_nhdr_files(args.inputDirectoryDWI)
mask_list = list_nhdr_files(args.inputDirectoryMask)

# create FA images
for (dwi, mask) in zip(dwi_list, mask_list):
    subject_id = os.path.splitext(os.path.basename(dwi))[0]
    fname_out_dti = os.path.join(args.outputDirectory, subject_id + '_DTI.nhdr')
    fname_out_b0 = os.path.join(args.outputDirectory, subject_id + '_B0.nhdr')
    fname_out_fa = os.path.join(args.outputDirectory, subject_id + '_FA.nhdr')
    print("/Applications/Slicer.app/Contents/lib/Slicer-4.4/cli-modules/DWIToDTIEstimation -m ", mask, dwi, fname_out_dti, fname_out_b0)
    print("/Applications/Slicer.app/Contents/lib/Slicer-4.4/cli-modules/DiffusionTensorScalarMeasurements", fname_out_dti, fname_out_fa, "-e FractionalAnisotropy &")

    

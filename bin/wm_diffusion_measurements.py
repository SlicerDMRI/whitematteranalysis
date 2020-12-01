#!/usr/bin/env python
import argparse
import os

try:
    import whitematteranalysis as wma
except:
    print("<wm_diffusion_measurements> Error importing white matter analysis package\n")
    raise

def main():
    #-----------------
    # Parse arguments
    #-----------------
    parser = argparse.ArgumentParser(
        description="Compute diffusion scalar measurements (such as FA, MD, etc). This script reports the mean statistics of each fiber cluster (or fiber tract) within the input folder.",
        epilog="Written by Fan Zhang (fzhang@bwh.harvard.edu)")
    
    parser.add_argument("-v", "--version",
        action="version", default=argparse.SUPPRESS,
        version='1.0',
        help="Show program's version number and exit")
    parser.add_argument(
        'inputDirectory',
        help='Directory of fiber clustering results obtained by <wm_cluster_from_altas.py> of multiple subjects. Make sure only the fiber clustering results are stored in this folder, making one subdirectory corresponding to one subject.')
    parser.add_argument(
        'outputCSV',
        help='Directory of output CSV files of fiber scalar measurement (computed using Slicer FiberTractMeasurements module).')
    parser.add_argument(
        'Slicer',
        help='Path of 3D Slicer.')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.inputDirectory):
        print("Error: Input directory", args.inputDirectory, "does not exist.")
        exit()
    
    if not os.path.exists(args.Slicer.split()[0]):
        print("Error: 3D Slicer", args.Slicer, "does not exist.")
        exit()
    
    module_FTSM = args.Slicer + ' '
    
    outdir = os.path.split(args.outputCSV)[0]
    if not os.path.exists(outdir):
        print("Output directory", outdir, "does not exist, creating it.")
        os.makedirs(outdir)
    
    print("<wm_diffusion_measurements>. Starting scalar measurement extraction.")
    print("")
    print("=====input directory======\n", args.inputDirectory)
    print("=====output directory=====\n", outdir)
    print("=====3D Slicer====\n", args.Slicer)
    print("==========================")
    
    os.system(module_FTSM + \
              ' --inputtype Fibers_File_Folder --format Column_Hierarchy --separator Comma ' + \
              ' --inputdirectory ' + args.inputDirectory + \
              ' --outputfile ' + args.outputCSV)
    
    print("<wm_diffusion_measurements> Measurements done at:", args.outputCSV)

if __name__ == '__main__':
    main()

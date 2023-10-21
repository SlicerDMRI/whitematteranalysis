#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import multiprocessing
import os

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Converts all vtp files in input directory to vtk files, which are saved in output directory.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu")
    parser.add_argument(
        'inputDirectory',
        help='Contains input tractography as vtkPolyData vtp format file(s).')
    parser.add_argument(
        'outputDirectory',
        help='The output directory should be a new empty directory. It will be created if needed. The actual output will be placed into a directory within this output directory, to preserve the informative name of the input directory that may contain subject ID, tract name, etc.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputDirectory):
        print(f"Error: Input directory {args.inputDirectory} does not exist.")
        exit()
    
    outdir = args.outputDirectory
    if not os.path.exists(args.outputDirectory):
        print(f"Output directory {outdir} does not exist, creating it.")
        os.makedirs(args.outputDirectory)
    
    # Preserve input directory name which may contain provenance like subject ID or structure, e.g UF
    subject_dir = os.path.splitext(os.path.basename(args.inputDirectory))[0]
    print(subject_dir)
    outdir_subject = os.path.join(args.outputDirectory, subject_dir+"_vtk")
    if not os.path.exists(outdir_subject):
        print(f"Output directory {outdir_subject} does not exist, creating it.")
        os.makedirs(outdir_subject)
    
    print("")
    print("wm_vtp2vtk.py: Convert all vtp files in input directory to vtk files in output directory.")
    print(f"=====input directory======\n {args.inputDirectory}")
    print(f"=====top-level output directory requested by user=====\n {args.outputDirectory}")
    print(f"=====final output directory=====\n {outdir_subject}")
    print("==========================")

    def list_vtp_files(input_dir):
        # Find input files (JUST vtp)
        #input_mask = "{0}/*.vtk".format(input_dir)
        input_mask2 = f"{input_dir}/*.vtp"
        #input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
        input_pd_fnames = glob.glob(input_mask2)
        input_pd_fnames = sorted(input_pd_fnames)
        return(input_pd_fnames)
    
    # Loop over input vtps
    inputPolyDatas = list_vtp_files(args.inputDirectory)
    
    print(f'Input number of vtp files found: {len(inputPolyDatas)}\n')
    
    if len(inputPolyDatas) < 1:
        print("Warning: no input vtp files found in input directory.")
        print("")
        exit()
        
    # for testing
    #inputPolyDatas = inputPolyDatas[0:2]
    
    for pd_fname in inputPolyDatas:
        
        # get subject/cluster identifier from unique input filename
        # -------------------
        subjectID = os.path.splitext(os.path.basename(pd_fname))[0]
    
        # read input vtk data
        # -------------------
        print(f"**Reading input: {pd_fname} ({subjectID})")
    
        wm = wma.io.read_polydata(pd_fname)
            
        # outputs
        # -------------------
        fname = os.path.join(outdir_subject, subjectID+'.vtk')
        try:
            print(f"====> Writing output polydata {fname}...")
            wma.io.write_polydata(wm, fname)
        except:
            print("Unknown exception in IO")
            raise
        del wm
    
    print("")
    print("Finished converting all vtp files to vtk files.")
    print("")
    
if __name__ == '__main__':
    main()

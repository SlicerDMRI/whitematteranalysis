#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import os
import shutil
import warnings

# If the cluster was created in bilateral clustering and is not commissural, separate into separate directories
# output should be a right_hem and a left_hem directory
# copy MRML files from toplevel directory into this one
# perhaps make a commissural directory also. why not.
# note that it could be cleaner to have a cell data array in the polydata and the option to 
# measure left and right hem parts of a cluster separately.
# Note: if this is performed at the atlas clustering stage, it can be used to separate clusters into groups,
# and this can be learned. At that point all data are midsagitally aligned, which this requires.
# For running this per-subject, the alignment should be performed to handle the tracts near the midline better.
# That should be added as an option.
import numpy as np
import vtk

import whitematteranalysis as wma
from whitematteranalysis.utils.utils import hemisphere_loc_name_typo_warn_msg


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Separate each cluster into left/right/commissural tracts based on the fiber location information computed by <wm_assess_cluster_location_by_hemisphere.py>. "
                    "The output is three directories of fiber bundles according to left hemisphere, right hemisphere, and commissural tracts. ",
        epilog="Written by Fan Zhang (fzhang@bwh.harvard.edu) and Lauren O\'Donnell (odonnell@bwh.harvard.edu). Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
    parser.add_argument(
        'inputDirectory',
        help='A directory of clustered whole-brain tractography as vtkPolyData (.vtk or .vtp).')
    parser.add_argument(
        'outputDirectory',
        help='The output directory will be created if it does not exist.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputDirectory):
        print(f"<{os.path.basename(__file__)}> Error: Input directory {args.inputDirectory} does not exist or is not a directory.")
        exit()
    
    outdir = args.outputDirectory
    if not os.path.exists(outdir):
        print(f"<{os.path.basename(__file__)}> Output directory {outdir} does not exist, creating it.")
        os.makedirs(outdir)
    
    outdir_right = os.path.join(outdir, 'tracts_right_hemisphere')
    if not os.path.exists(outdir_right):
        os.makedirs(outdir_right)
    outdir_left = os.path.join(outdir, 'tracts_left_hemisphere')
    if not os.path.exists(outdir_left):
        os.makedirs(outdir_left)
    outdir_commissure = os.path.join(outdir, 'tracts_commissural')
    if not os.path.exists(outdir_commissure):
        os.makedirs(outdir_commissure)
    
    
    def list_cluster_files(input_dir):
        # Find input files
        input_mask = f"{input_dir}/cluster_*.vtk"
        input_mask2 = f"{input_dir}/cluster_*.vtp"
        input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
        input_pd_fnames = sorted(input_pd_fnames)
        return(input_pd_fnames)
    
    def write_mask_location_to_vtk(inpd, mask_location):
    
        inpointdata = inpd.GetPointData()
        if inpointdata.GetNumberOfArrays() > 0:
            point_data_array_indices = list(range(inpointdata.GetNumberOfArrays()))
            for idx in point_data_array_indices:
                array = inpointdata.GetArray(idx)
                if array.GetName() == 'HemisphereLocataion':
                    warnings.warn(
                        hemisphere_loc_name_typo_warn_msg,
                        PendingDeprecationWarning)
                    print('  -- HemisphereLocataion is in the input data: skip updating the vtk file.')
                    return inpd
                elif array.GetName() == 'HemisphereLocation':
                    print('  -- HemisphereLocation is in the input data: skip updating the vtk file.')
                    return inpd
    
        vtk_array = vtk.vtkDoubleArray()
        vtk_array.SetName('HemisphereLocation')
    
        inpd.GetLines().InitTraversal()
        for lidx in range(0, inpd.GetNumberOfLines()):
            ptids = vtk.vtkIdList()
            inpd.GetLines().GetNextCell(ptids)
            for pidx in range(0, ptids.GetNumberOfIds()):
                vtk_array.InsertNextTuple1(mask_location[lidx])
    
        inpd.GetPointData().AddArray(vtk_array)
        inpd.Update()
    
        return inpd
    
    def read_mask_location_from_vtk(inpd):

        def _read_location(_inpd):
            _flag_location = True

            _inpd.GetLines().InitTraversal()
            for lidx in range(0, _inpd.GetNumberOfLines()):
                ptids = vtk.vtkIdList()
                _inpd.GetLines().GetNextCell(ptids)
                for pidx in range(0, ptids.GetNumberOfIds()):
                    mask_location[lidx] = array.GetTuple(ptids.GetId(pidx))[0]

            return _flag_location

        mask_location = np.zeros(pd.GetNumberOfLines())
    
        inpointdata = inpd.GetPointData()
        flag_location = False
        if inpointdata.GetNumberOfArrays() > 0:
            point_data_array_indices = list(range(inpointdata.GetNumberOfArrays()))
            for idx in point_data_array_indices:
                array = inpointdata.GetArray(idx)
                if array.GetName() == 'HemisphereLocataion':
                    warnings.warn(
                        hemisphere_loc_name_typo_warn_msg,
                        PendingDeprecationWarning)
                    flag_location = _read_location(inpd)
                    break
                elif array.GetName() == 'HemisphereLocation':
                    flag_location = _read_location(inpd)
                    break
    
        return flag_location, mask_location
    
    print(f"<{os.path.basename(__file__)}> Starting computation.")
    print("")
    print(f"=====input directory ======\n {args.inputDirectory}")
    print(f"=====output directory =====\n {args.outputDirectory}")
    print("==========================")
    print("")
    
    input_polydatas = list_cluster_files(args.inputDirectory)
    
    number_of_clusters = len(input_polydatas)
    
    print(f"<{os.path.basename(__file__)}> Input number of vtk/vtp files: {number_of_clusters}")
    
    # midsagittal alignment step to get transform to apply to each polydata
    # alternatively could use the transform into atlas space that is found
    # during the labeling. This should be an option at the labeling step.
    # this transform is the one to use.
    
    # read in data
    for fname, c_idx in zip(input_polydatas, list(range(len(input_polydatas)))):
    
        # figure out filename and extension
        fname_base = os.path.basename(fname)
    
        # read data
        print(f"<{os.path.basename(__file__)}> Separating input file: {fname}")
        pd = wma.io.read_polydata(fname)
    
        flag_location, mask_location = read_mask_location_from_vtk(pd)
           
         # If HemisphereLocation is not defined in the input vtk file, the location of each fiber in the cluster is decided.
        if not flag_location:
            print(f"Error: {fname} has no hemisphere location information")
            exit()
    
        # for sanity check
        if len(np.where(mask_location ==0)[0]) > 1:
            print(f"Error: Not all fibers in {fname} is labeled with hemisphere location information.")
            exit()
    
        # output separated clusters
        mask_right = np.zeros(pd.GetNumberOfLines())
        mask_left = np.zeros(pd.GetNumberOfLines())
        mask_commissure = np.zeros(pd.GetNumberOfLines())
    
        mask_left[np.where(mask_location==1)[0]] = 1
        mask_right[np.where(mask_location==2)[0]] = 1
        mask_commissure[np.where(mask_location==3)[0]] = 1
    
        pd_right = wma.filter.mask(pd, mask_right, preserve_point_data=True, preserve_cell_data=True, verbose=False)
        pd_left = wma.filter.mask(pd, mask_left, preserve_point_data=True, preserve_cell_data=True, verbose=False)
        pd_commissure = wma.filter.mask(pd, mask_commissure, preserve_point_data=True, preserve_cell_data=True, verbose=False)
    
        fname_output = os.path.join(outdir_right, fname_base)
        wma.io.write_polydata(pd_right, fname_output)
        fname_output = os.path.join(outdir_left, fname_base)
        wma.io.write_polydata(pd_left, fname_output)
        fname_output = os.path.join(outdir_commissure, fname_base)
        wma.io.write_polydata(pd_commissure, fname_output)
    
    print("")
    print(f"<{os.path.basename(__file__)}> Done!!!")

if __name__ == '__main__':
    main()

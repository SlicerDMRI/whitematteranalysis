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
        description="Assess if a fiber within the clusters belongs to left hemispheric, right hemispheric, or commissural tracts. "
                    "This code needs to be run in the ATLAS space, where the brain is clustered at the RAS origin. ",
        epilog="Written by Fan Zhang (fzhang@bwh.harvard.edu) and Lauren O\'Donnell (odonnell@bwh.harvard.edu). Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
    parser.add_argument(
        'inputDirectory',
        help='A directory of clustered whole-brain tractography as vtkPolyData (.vtk or .vtp).')
    parser.add_argument(
        '-pthresh', action="store", dest="hemispherePercentThreshold", type=float,
        help='The percent of a fiber that has to be in one hemisphere to consider the fiber as part of that hemisphere (rather than a commissural fiber). '
             'Default number is 0.6, where a higher number will tend to label fewer fibers as hemispheric and more fibers as commissural (not strictly in one hemisphere or the other), '
             'while a lower number will be stricter about what is classified as commissural. '
             'This parameter is not applicable if -clusterLocationFile was provided. '
             'Note: performing hemisphere separation using -pthresh needs the input clusters be in the ATLAS space.')
    parser.add_argument(
        '-clusterLocationFile', action="store", dest="clusterLocationFile",
        help='A csv file defining the location of each cluster, i.e., hemispheric or commissural. '
             'The hemispherePercentThreshold will be varied across the clusters if their locations are already known.'
             'Note: performing hemisphere separation using -clusterLocationFile needs the input clusters be in the ATLAS space.')
    parser.add_argument(
        '-outputDirectory',
        help='If this is given, separated clusters will be output under this folder. The output directory will be created if it does not exist.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

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
        # inpd.Update()
    
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

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.isdir(args.inputDirectory):
        print(f"<{os.path.basename(__file__)}> Error: Input directory {args.inputDirectory} does not exist or is not a directory.")
        exit()
    
    clusterLocationFile = args.clusterLocationFile
    location_data = None
    if clusterLocationFile is not None:
        print(clusterLocationFile)
        if not os.path.exists(clusterLocationFile):
            print(f"<{os.path.basename(__file__)}> Cluster location file is assigned but the file", clusterLocationFile, "does not exist.")
            exit()
        else:
            location_data = np.loadtxt(open(clusterLocationFile, "rb"),
                                          dtype={'names': ('Cluster Index', 'Location Label'), 'formats': ('S17', 'S1')},
                                          delimiter="\t", skiprows=1)

            # change CPC to be hemispheric
            location_data[145][1] = b'h'
            location_data[159][1] = b'h'
            location_data[557][1] = b'h'
            location_data[677][1] = b'h'
            location_data[770][1] = b'h'
    
    # default to be changed if user input is there
    hemisphere_percent_threshold = 0.6
    if args.hemispherePercentThreshold is not None:
        if (args.hemispherePercentThreshold > 0.5) & (args.hemispherePercentThreshold <= 1.0):
            hemisphere_percent_threshold = args.hemispherePercentThreshold
        else:
            print(f"<{os.path.basename(__file__)}> Hemisphere fiber percent threshold", args.hemispherePercentThreshold, "must be between 0.5 and 1. (0.6 is recommended).")
            exit()
    
    outdir = args.outputDirectory
    if outdir is not None:
        print(f"<{os.path.basename(__file__)}> Separated clusters will be ouput.")
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
    
    
    print(f"<{os.path.basename(__file__)}> Starting computation.")
    print("")
    print("=====input directory ======\n", args.inputDirectory)
    print("==========================")
    print("")
    
    if location_data is None:
        print(f"<{os.path.basename(__file__)}> Hemisphere fiber percent threshold", hemisphere_percent_threshold)
    else:
        print(f"<{os.path.basename(__file__)}> Separating clusters using location file:", clusterLocationFile)
    print("")
    
    # relatively high number of points for accuracy
    points_per_fiber = 40
    
    input_polydatas = list_cluster_files(args.inputDirectory)
    
    number_of_clusters = len(input_polydatas)
    
    if location_data is not None and number_of_clusters != len(location_data):
        print(f"<{os.path.basename(__file__)}> Error: Number of clusters (%d) is not the same as in the location file (%d)." % (number_of_clusters, len(location_data)))
        exit()
    
    print(f"<{os.path.basename(__file__)}> Input number of vtk/vtp files: ", number_of_clusters)
    
    # midsagittal alignment step to get transform to apply to each polydata
    # alternatively could use the transform into atlas space that is found
    # during the labeling. This should be an option at the labeling step.
    # this transform is the one to use.
    
    # read in data
    for fname, c_idx in zip(input_polydatas, list(range(len(input_polydatas)))):
    
        # figure out filename and extension
        fname_base = os.path.basename(fname)
    
        # read data
        print(f"<{os.path.basename(__file__)}> Separating input file:", fname)
        pd = wma.io.read_polydata(fname)
    
        flag_location, mask_location = read_mask_location_from_vtk(pd)
           
         # If HemisphereLocation is not defined in the input vtk file, the location of each fiber in the cluster is decided.
        if not flag_location:
            if clusterLocationFile is None:
                # internal representation for fast similarity computation
                # this also detects which hemisphere fibers are in
                fibers = wma.fibers.FiberArray()
                fibers.points_per_fiber = points_per_fiber
                fibers.hemisphere_percent_threshold = hemisphere_percent_threshold
                # must request hemisphere computation from object
                fibers.hemispheres = True
                # Now convert to array with points and hemispheres as above
                fibers.convert_from_polydata(pd)
    
                # separate into right and left hemispheres
                # note: this assumes RAS coordinates as far as right and left labels
                # LPS would be switched
                # -------------------------
    
                mask_location[fibers.index_left_hem] = 1
                mask_location[fibers.index_right_hem] = 2
                mask_location[fibers.index_commissure] = 3
            else:
                if location_data[c_idx][1] == b'c' or location_data[c_idx][1] == b'ng':
    
                    mask_location[:] = 3
    
                elif location_data[c_idx][1] == b'h':
                    hemisphere_percent_threshold = 0.5001
    
                    # internal representation for fast similarity computation
                    # this also detects which hemisphere fibers are in
                    fibers = wma.fibers.FiberArray()
                    fibers.points_per_fiber = points_per_fiber
                    fibers.hemisphere_percent_threshold = hemisphere_percent_threshold
                    # must request hemisphere computation from object
                    fibers.hemispheres = True
                    # Now convert to array with points and hemispheres as above
                    fibers.convert_from_polydata(pd)
    
                    # separate into right and left hemispheres
                    # note: this assumes RAS coordinates as far as right and left labels
                    # LPS would be switched
                    # -------------------------
                  
                    mask_location[fibers.index_left_hem] = 1
                    mask_location[fibers.index_right_hem] = 2
    
                    if len(fibers.index_commissure) > 0:
                        if len(fibers.index_left_hem) <= len(fibers.index_right_hem):
                            mask_location[fibers.index_commissure] = 1
                        else:
                            mask_location[fibers.index_commissure] = 2
        
             # Update the input vtk file
            pd = write_mask_location_to_vtk(pd, mask_location)
            wma.io.write_polydata(pd, fname)
    
        # for sanity check 
        if len(np.where(mask_location ==0)[0]) > 1:
            print("Error: Not all fibers in", fname, "is labeled with hemisphere location infromation.")
            exit()
    
        if outdir is not None:
            
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
    
    f = open(os.path.join(args.inputDirectory, 'cluster_location_by_hemisphere.log'), 'a')
    f.write('<wm_assess_cluster_location_by_hemisphere.py> Done!!!')
    f.close()
    
    print("")
    print(f"<{os.path.basename(__file__)}> Done!!!")

if __name__ == '__main__':
    main()

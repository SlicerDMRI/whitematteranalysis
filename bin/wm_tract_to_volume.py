#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import nibabel as nib
import numpy as np
import vtk
from nibabel.affines import apply_affine

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Convert a fiber tract or cluster (vtk) to a voxel-wise fiber density image (nii.gz). ",
        epilog="Written by Fan Zhang")
    parser.add_argument(
        'inputVTK',
        help='Input VTK/VTP file that is going to be converted.')
    parser.add_argument(
        'refvolume',
        help='A volume image that the cluster will be converted on.')
    parser.add_argument(
        'outputVol',
        help='Output volume image, where the value of each voxel represents the number of fibers passing though the voxel.')
    parser.add_argument(
        '-m', action="store", dest="measure", type=str,
        help="diffusion measure; if not provided, the output will be density map")

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.exists(args.inputVTK):
        print(f"Error: Input directory {args.inputVTK} does not exist.")
        exit()

    def convert_cluster_to_volume_with_sz(inpd, volume, sampling_size=0.5):
    
        volume_shape = volume.get_fdata().shape
        new_voxel_data = np.zeros(volume_shape)
    
        resampler = vtk.vtkPolyDataPointSampler()
        resampler.GenerateEdgePointsOn()
        resampler.GenerateVertexPointsOff()
        resampler.GenerateInteriorPointsOff()
        resampler.GenerateVerticesOff()
        resampler.SetDistance(sampling_size)
    
        inpoints = inpd.GetPoints()
    
        inpd.GetLines().InitTraversal()
        for lidx in range(0, inpd.GetNumberOfLines()):
    
            ptids = vtk.vtkIdList()
            inpd.GetLines().GetNextCell(ptids)
    
            tmpPd = vtk.vtkPolyData()
            tmpPoints = vtk.vtkPoints()
            tmpCellPtIds = vtk.vtkIdList()
            tmpLines =  vtk.vtkCellArray()
            
            for pidx in range(0, ptids.GetNumberOfIds()):
                point = inpoints.GetPoint(ptids.GetId(pidx))
                idx_ = tmpPoints.InsertNextPoint(point)
                tmpCellPtIds.InsertNextId(idx_)
    
            tmpLines.InsertNextCell(tmpCellPtIds)
    
            tmpPd.SetLines(tmpLines)
            tmpPd.SetPoints(tmpPoints)
    
            if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                resampler.SetInputData(tmpPd)
            else:
                resampler.SetInput(tmpPd)
    
            resampler.Update()
    
            sampledCellPts = resampler.GetOutput().GetPoints()
            sampledNpts = resampler.GetOutput().GetNumberOfPoints()
    
            line_tmp_voxel_data = np.zeros(volume_shape)
            for pidx in range(0, sampledNpts):
                point = sampledCellPts.GetPoint(pidx)
    
                point_ijk = apply_affine(np.linalg.inv(volume.affine), point)
                point_ijk = np.rint(point_ijk).astype(np.int32)
    
                line_tmp_voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])] += 1
    
            new_voxel_data = new_voxel_data + line_tmp_voxel_data
    
        return new_voxel_data
    
    
    def convert_cluster_to_volume(inpd, volume, measure=None):
    
        volume_shape = volume.get_fdata().shape
        new_voxel_data = np.zeros(volume_shape)
        new_voxel_measure = np.zeros(volume_shape)
    
        inpoints = inpd.GetPoints()
        if measure is not None:
            label_array = inpd.GetPointData().GetArray(measure)
            if label_array is None:
                print("Error: Check if the input measure (-m) is stored in the vtk file.")
                exit()

        inpd.GetLines().InitTraversal()
        for lidx in range(0, inpd.GetNumberOfLines()):
    
            ptids = vtk.vtkIdList()
            inpd.GetLines().GetNextCell(ptids)
            # print 'Line', lidx, ' - ', ptids.GetNumberOfIds()
            for pidx in range(0, ptids.GetNumberOfIds()):
                point = inpoints.GetPoint(ptids.GetId(pidx))
    
                point_ijk = apply_affine(np.linalg.inv(volume.affine), point)
                point_ijk = np.rint(point_ijk).astype(np.int32)

                if measure is not None:
                    count = new_voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])]
                    pre_val = new_voxel_measure[(point_ijk[0], point_ijk[1], point_ijk[2])]
                    val = label_array.GetTuple(ptids.GetId(pidx))[0]
                    new_voxel_measure[(point_ijk[0], point_ijk[1], point_ijk[2])] = (pre_val * count + val) / (count + 1)

                new_voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])] += 1

        if measure is not None:
            new_voxel_data = new_voxel_measure
            
        return new_voxel_data

    volume = nib.load(args.refvolume)
    print(f'<{os.path.basename(__file__)}> {args.refvolume} input volume shape: f{volume.get_fdata().shape}')

    inpd = wma.io.read_polydata(args.inputVTK)
    
    new_voxel_data = convert_cluster_to_volume(inpd, volume, measure=args.measure)
    
    volume_new = nib.Nifti1Image(new_voxel_data, volume.affine, volume.header)
    
    nib.save(volume_new, args.outputVol)
    
    print(f'Done: save tract map to {args.outputVol}')

if __name__ == '__main__':
    main()

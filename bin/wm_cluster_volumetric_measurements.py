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
        description="Compute measurements of a cluster based on a provided volumetric scalar map (e.g. an FA image).",
        epilog="Written by Fan Zhang")
    parser.add_argument(
        'inputDirectory',
        help='Directory of input VTK/VTP files.')
    parser.add_argument(
        'inputVolumetricMap',
        help='Path to input volumetric map, e.g. an FA image. Note: this input image needs to be a nifti (nii or nii.gz) file. To read and write this format, NiBabel package is needed, by running: pip install nibabel ')
    parser.add_argument(
        'outputDirectory',
        help='Directory of output statistics.')
    parser.add_argument(
        'outputStatFile', action="store", type=str,
        help="File name of the output statistics.")
    parser.add_argument(
        '-sampleSize', action="store", type=float,
        help='Fiber sample size')
    parser.add_argument(
        '-outputLabelmap', action='store_true',
        help='Generate a label map of each input cluster if given.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    inputdir = os.path.abspath(args.inputDirectory)
    if not os.path.isdir(args.inputDirectory):
        print(f"Error: Input directory {args.inputDirectory} does not exist.")
        exit()
    
    inputvol = os.path.abspath(args.inputVolumetricMap)
    if not os.path.exists(inputvol):
        print(f"Error: Input volumetric map {inputvol} does not exist.")
        exit()
    
    outdir = os.path.abspath(args.outputDirectory)
    if not os.path.exists(args.outputDirectory):
        print(f"Output directory {args.outputDirectory} does not exist, creating it.")
        os.makedirs(outdir)

    input_volume = nib.load(inputvol)
    print(f'<{os.path.basename(__file__)}> Input volume shape: {input_volume.get_data().shape}')

    input_vtk_list = wma.io.list_vtk_files(inputdir)
    print(f'<{os.path.basename(__file__)}> Number of input clusters: {len(input_vtk_list)}')
    
    output_stats_file = os.path.join(outdir, args.outputStatFile)
    
    # compute cluster statistics based on the input volumetric map
    def compute_stat(inpd, volume, sampling_size=None):
    
        volume_data = volume.get_data()
    
        volume_shape = volume.get_data().shape 
        voxel_size = volume.header.get_zooms()
    
        new_voxel_data = np.zeros(volume_shape)
    
        if sampling_size is not None:
            resampler = vtk.vtkPolyDataPointSampler()
            resampler.GenerateEdgePointsOn()
            resampler.GenerateVertexPointsOff()
            resampler.GenerateInteriorPointsOff()
            resampler.GenerateVerticesOff()
            resampler.SetDistance(sampling_size)
    
        value_list = []
        inpoints = inpd.GetPoints()
        inpd.GetLines().InitTraversal()
        for lidx in range(0, inpd.GetNumberOfLines()):
            ptids = vtk.vtkIdList()
            inpd.GetLines().GetNextCell(ptids)
    
            if sampling_size is not None:
                # print '-- sampling with size', sampling_size
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
            else:
                # print '-- No sampling'
                sampledCellPts = inpoints
                sampledNpts = inpd.GetNumberOfPoints()
    
            for pidx in range(0, sampledNpts):
                point = sampledCellPts.GetPoint(pidx)
    
                point_ijk = apply_affine(np.linalg.inv(volume.affine), point)
                point_ijk = np.rint(point_ijk).astype(np.int32)
    
                if point_ijk[0] > new_voxel_data.shape[0] or point_ijk[1]> new_voxel_data.shape[1] or point_ijk[2] > new_voxel_data.shape[2]:
                    print(f'Warning: point {point_ijk} is outside the input volumetric map.')
                    continue
    
                if new_voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])] == 0:
                    new_voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])] = 1
                    value_list.append(volume_data[(point_ijk[0], point_ijk[1], point_ijk[2])])
    
        if len(value_list) > 0:
            mean_v = np.mean(value_list)
            var_v = np.var(value_list)
            max_v = np.max(value_list)
            min_v = np.min(value_list)
            median_v = np.median(value_list)
        else:
            mean_v = np.nan
            var_v = np.nan
            max_v = np.nan
            min_v = np.nan
            median_v = np.nan
         
        num_voxels = len(value_list)
        volume_size = voxel_size[0] * voxel_size[1] * voxel_size[2] * num_voxels
    
        return new_voxel_data, mean_v, var_v, max_v, min_v, median_v, num_voxels, volume_size
    
    str_head = 'Name,Num_Voxels,Volume,Mean,Median,Min,Max,Variance'
    str_out = str_head
    for input_vtk_path in input_vtk_list:
    
        input_vtk = wma.io.read_polydata(input_vtk_path)
    
        vtk_file_name = os.path.split(input_vtk_path)[1][:-4]
        print(f"<{os.path.basename(__file__)}> Working on {vtk_file_name}")
    
        new_voxel_data, mean_v, var_v, max_v, min_v, median_v, num_voxels, volume_size = compute_stat(input_vtk, input_volume, args.sampleSize)
    
        str_line = f"{vtk_file_name},{str(num_voxels)},{str(volume_size)},{str(mean_v)},{str(median_v)},{str(min_v)},{str(max_v)},{str(var_v)}"
        print(f' - {str_head}')
        print(f' - {str_line}')
    
        str_out = str_out + '\n' + str_line
    
        if args.outputLabelmap:
            volume_new = nib.Nifti1Image(new_voxel_data, input_volume.affine, input_volume.header)
            output_labelmap = os.path.join(outdir, f'{vtk_file_name}.nii.gz')
            nib.save(volume_new, output_labelmap)
    
    output_file = open(output_stats_file, 'w')
    output_file.write(str_out)
    output_file.close()
    
    print('')
    print(f'<{os.path.basename(__file__)}> Done! Result is in {output_stats_file}')

if __name__ == '__main__':
    main()

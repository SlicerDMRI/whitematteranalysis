#!/usr/bin/env python

import argparse
import os
import nibabel
import vtk
import numpy
from nibabel.affines import apply_affine

try:
    import whitematteranalysis as wma
except:
    print("Error importing white matter analysis package\n")
    raise

def main():
    #-----------------
    # Parse arguments
    #-----------------
    parser = argparse.ArgumentParser(
        description="Compute measurements of a cluster based on a provided volumetric scalar map (e.g. an FA image).",
        epilog="Written by Fan Zhang")
    
    parser.add_argument("-v", "--version",
        action="version", default=argparse.SUPPRESS,
        version='1.0',
        help="Show program's version number and exit")
    
    parser.add_argument(
        'inputDirectory',
        help='Directory of input VTK/VTP files.')
    
    parser.add_argument(
        'inputVolumetricMap',
        help='Path to input volumetric map, e.g. an FA image. Note: this input image needs to be a nifti (nii or nii.gz) file. To read and writ this format, NiBabel package is needed, by runing: pip install nibabel ')
    
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
    
    args = parser.parse_args()
    
    inputdir = os.path.abspath(args.inputDirectory)
    if not os.path.isdir(args.inputDirectory):
        print("Error: Input directory", args.inputDirectory, "does not exist.")
        exit()
    
    inputvol = os.path.abspath(args.inputVolumetricMap)
    if not os.path.exists(inputvol):
        print("Error: Input volumetri map", inputvol, "does not exist.")
        exit()
    
    outdir = os.path.abspath(args.outputDirectory)
    if not os.path.exists(args.outputDirectory):
        print("Output directory", args.outputDirectory, "does not exist, creating it.")
        os.makedirs(outdir)
    
    input_volume = nibabel.load(inputvol)
    print('<wm_cluster_volumetric_measurements.py> Input volume shape:', input_volume.get_data().shape)
    
    input_vtk_list = wma.io.list_vtk_files(inputdir)
    print('<wm_cluster_volumetric_measurements.py> Number of input clusters:', len(input_vtk_list))
    
    output_stats_file = os.path.join(outdir, args.outputStatFile)
    
    # compute cluster statistics based on the input volumetric map
    def compute_stat(inpd, volume, sampling_size=None):
    
        volume_data = volume.get_data()
    
        volume_shape = volume.get_data().shape 
        voxel_size = volume.header.get_zooms()
    
        new_voxel_data = numpy.zeros(volume_shape)
    
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
    
                point_ijk = apply_affine(numpy.linalg.inv(volume.affine), point)
                point_ijk = numpy.rint(point_ijk).astype(numpy.int32)
    
                if point_ijk[0] > new_voxel_data.shape[0] or point_ijk[1]> new_voxel_data.shape[1] or point_ijk[2] > new_voxel_data.shape[2]:
                    print('Warning: point', point_ijk, 'is outside the input volumetric map.')
                    continue
    
                if new_voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])] == 0:
                    new_voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])] = 1
                    value_list.append(volume_data[(point_ijk[0], point_ijk[1], point_ijk[2])])
    
        if len(value_list) > 0:
            mean_v = numpy.mean(value_list)
            var_v = numpy.var(value_list)
            max_v = numpy.max(value_list)
            min_v = numpy.min(value_list)
            median_v = numpy.median(value_list)
        else:
            mean_v = numpy.nan
            var_v = numpy.nan
            max_v = numpy.nan
            min_v = numpy.nan
            median_v = numpy.nan
         
        num_voxels = len(value_list)
        volume_size = voxel_size[0] * voxel_size[1] * voxel_size[2] * num_voxels
    
        return new_voxel_data, mean_v, var_v, max_v, min_v, median_v, num_voxels, volume_size
    
    str_head = 'Name,Num_Voxels,Volume,Mean,Median,Min,Max,Variance'
    str_out = str_head
    for input_vtk_path in input_vtk_list:
    
        input_vtk = wma.io.read_polydata(input_vtk_path)
    
        vtk_file_name = os.path.split(input_vtk_path)[1][:-4]
        print("<wm_cluster_volumetric_measurements.py> Working on", vtk_file_name)
    
        new_voxel_data, mean_v, var_v, max_v, min_v, median_v, num_voxels, volume_size = compute_stat(input_vtk, input_volume, args.sampleSize)
    
        str_line = vtk_file_name+','+str(num_voxels)+','+str(volume_size)+','+str(mean_v)+','+str(median_v)+','+str(min_v)+','+str(max_v)+','+str(var_v)
        print(' -', str_head)
        print(' -', str_line)
    
        str_out = str_out + '\n' + str_line
    
        if args.outputLabelmap:
            volume_new = nibabel.Nifti1Image(new_voxel_data, input_volume.affine, input_volume.header)
            output_labelmap = os.path.join(outdir, vtk_file_name+'.nii.gz')
            nibabel.save(volume_new, output_labelmap)
    
    output_file = open(output_stats_file, 'w')
    output_file.write(str_out)
    output_file.close()
    
    print('')
    print('<wm_cluster_volumetric_measurements.py> Done! Result is in', output_stats_file)

if __name__ == '__main__':
    main()

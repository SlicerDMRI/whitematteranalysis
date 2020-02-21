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
        description="Convert a fiber tract or cluster (vtk) to a voxel-wise fiber density image (nii.gz). ",
        epilog="Written by Fan Zhang")
    
    parser.add_argument("-v", "--version",
        action="version", default=argparse.SUPPRESS,
        version='1.0',
        help="Show program's version number and exit")
    
    parser.add_argument(
        'inputVTK',
        help='Input VTK/VTP file that is going to be converted.')
    parser.add_argument(
        'refvolume',
        help='A volume image that the cluster will be converted on.')
    parser.add_argument(
        'outputVol',
        help='Output volume image, where the value of each voxel represents the number of fibers passing though the voxel.')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.inputVTK):
        print("Error: Input directory", args.inputVTK, "does not exist.")
        exit()
    
    def convert_cluster_to_volume_with_sz(inpd, volume, sampling_size=0.5):
    
        volume_shape = volume.get_data().shape
        new_voxel_data = numpy.zeros(volume_shape)
    
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
    
            line_tmp_voxel_data = numpy.zeros(volume_shape)
            for pidx in range(0, sampledNpts):
                point = sampledCellPts.GetPoint(pidx)
    
                point_ijk = apply_affine(numpy.linalg.inv(volume.affine), point)
                point_ijk = numpy.rint(point_ijk).astype(numpy.int32)
    
                line_tmp_voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])] += 1
    
            new_voxel_data = new_voxel_data + line_tmp_voxel_data
    
        return new_voxel_data
    
    
    def convert_cluster_to_volume(inpd, volume):
    
        volume_shape = volume.get_data().shape
        new_voxel_data = numpy.zeros(volume_shape)
    
        inpoints = inpd.GetPoints()
    
        inpd.GetLines().InitTraversal()
        for lidx in range(0, inpd.GetNumberOfLines()):
    
            ptids = vtk.vtkIdList()
            inpd.GetLines().GetNextCell(ptids)
            # print 'Line', lidx, ' - ', ptids.GetNumberOfIds()
            for pidx in range(0, ptids.GetNumberOfIds()):
                point = inpoints.GetPoint(ptids.GetId(pidx))
    
                point_ijk = apply_affine(numpy.linalg.inv(volume.affine), point)
                point_ijk = numpy.rint(point_ijk).astype(numpy.int32)
    
                new_voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])] += 1
    
        return new_voxel_data
    
    volume = nibabel.load(args.refvolume)
    print('<wm_tract_to_volume>', args.refvolume, ', input volume shape: ', volume.get_data().shape)
    
    inpd = wma.io.read_polydata(args.inputVTK)
    
    new_voxel_data = convert_cluster_to_volume_with_sz(inpd, volume)
    
    volume_new = nibabel.Nifti1Image(new_voxel_data, volume.affine, volume.header)
    
    nibabel.save(volume_new, args.outputVol)
    
    print('Done: save tract map to', args.outputVol)

if __main__ == '__main__':
    main()

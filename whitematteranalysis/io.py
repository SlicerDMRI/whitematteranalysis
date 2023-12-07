# -*- coding: utf-8 -*-

""" io.py

This module provides input of vtk polydata tractography files (vtk/vtp).
It also provides a class for I/O of laterality results.

read_polydata

Function to read vtkPolyData in .vtk or .vtp form

write_laterality_results

Function to write laterality indices, histograms, polydata to summarize
laterality output

read_laterality_results

This function reads in the laterality data for further analysis.

"""

import glob
import os
import pickle
import time

import numpy as np
import vtk
from joblib import Parallel, delayed

from . import filter, render

VERBOSE = 0


def read_polydata(filename):
    """Read whole-brain tractography as vtkPolyData format."""

    if VERBOSE:
        print(f"Reading in data from {filename}...")

    basename, extension = os.path.splitext(filename)

    if   (extension == '.vtk'):
        reader = vtk.vtkPolyDataReader()
    elif (extension == '.vtp'):
        reader = vtk.vtkXMLPolyDataReader()
    else:
        print('Cannot recognize model file format')
        return None

    reader.SetFileName(filename)
    reader.Update()
    outpd = reader.GetOutput()
    del reader
    if VERBOSE:
        print(f"Done reading in data from {filename}")
        print(f"Number of lines found: {outpd.GetNumberOfLines()}")

    return outpd

def list_vtk_files(input_dir):
    # Find input files
    input_mask = f"{input_dir}/*.vtk"
    input_mask2 = f"{input_dir}/*.vtp"
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)

def list_transform_files(input_dir):
    # Find input files
    input_mask = f"{input_dir}/*.tfm"
    input_tf_fnames = glob.glob(input_mask)
    input_tf_fnames = sorted(input_tf_fnames)
    return (input_tf_fnames)

def read_and_preprocess_polydata_directory(input_dir, fiber_length, number_of_fibers, random_seed=None, fiber_length_max=None):
    """ Find and read all .vtk and .vtp files in the given directory
    input_dir. Preprocess with fiber length threshold and downsample
    to desired number of fibers."""
    
    input_pd_fnames = list_vtk_files(input_dir)
    num_pd = len(input_pd_fnames)
    
    print(f"<{os.path.basename(__file__)}> =======================================")
    print(f"<{os.path.basename(__file__)}> Reading vtk and vtp files from directory: {input_dir}")
    print(f"<{os.path.basename(__file__)}> Total number of files found: {num_pd}")
    print(f"<{os.path.basename(__file__)}> =======================================")

    input_pds = list()
    subject_ids = list()
    sidx = 0

    for fname in input_pd_fnames:
        subject_id = os.path.splitext(os.path.basename(fname))[0]
        subject_ids.append(subject_id)
        print(f"<{os.path.basename(__file__)}> {sidx + 1} / {num_pd} {subject_id} Reading {fname}...")
        pd = read_polydata(fname)
        print(f"<{os.path.basename(__file__)}> {sidx + 1} / {num_pd} {subject_id} Input number of fibers: {pd.GetNumberOfLines()}")
        pd2 = filter.preprocess(pd, min_length_mm=fiber_length, verbose=False, max_length_mm=fiber_length_max)
        print(f"<{os.path.basename(__file__)}> {sidx + 1} / {num_pd} {subject_id} Length threshold {fiber_length} mm. Number of fibers retained: {pd2.GetNumberOfLines()}")
        pd3 = filter.downsample(pd2, number_of_fibers, verbose=False, random_seed=random_seed)
        print(f"<{os.path.basename(__file__)}> {sidx + 1} / num_pd {subject_id} Downsample to {number_of_fibers} fibers. Number of fibers retained: {pd3.GetNumberOfLines()}")
        input_pds.append(pd3)
        sidx += 1
        print(f"<{os.path.basename(__file__)}> =======================================")

    print(f"<{os.path.basename(__file__)}> =======================================")
    print(f"<{os.path.basename(__file__)}> Done reading vtk and vtp files from directory: {input_dir}")
    print(f"<{os.path.basename(__file__)}> Total number of files read: {len(input_pds)}")
    print(f"<{os.path.basename(__file__)}> =======================================")
        
    return input_pds, subject_ids

    
                            
def write_polydata(polydata, filename):
    """Write polydata as vtkPolyData format, according to extension."""

    if VERBOSE:
        print("Writing ", filename, "...")

    basename, extension = os.path.splitext(filename)

    if   (extension == '.vtk'):
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileTypeToBinary()
    elif (extension == '.vtp'):
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetDataModeToBinary()
    else:
        print('Cannot recognize model file format')
        return None

    writer.SetFileName(filename)
    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        writer.SetInputData(polydata)
    else:
        writer.SetInput(polydata)
    writer.Update()

    del writer

    if VERBOSE:
        print(f"Done writing {filename}")

def transform_polydata_from_disk(in_filename, transform_filename, out_filename):
    # Read it in.
    print(f"<{os.path.basename(__file__)}> Transforming {in_filename} -> {out_filename}...")

    # Read the transform from disk because we cannot pickle it
    (root, ext) = os.path.splitext(transform_filename)
    print(root, ext)
    if ext == '.xfm':
        reader = vtk.vtkMNITransformReader()
        reader.SetFileName(transform_filename)
        reader.Update()
        transform = reader.GetTransform()
    elif ext == '.img':
        reader = vtk.vtkImageReader()
        reader.SetFileName(transform_filename)
        reader.Update()
        coeffs = reader.GetOutput()
        transform = vtk.vtkBSplineTransform()
        transform.SetCoefficients(coeffs)
        print(coeffs)
        print(transform)
    else:
        f = open(transform_filename)
        transform = vtk.vtkTransform()
        matrix = vtk.vtkMatrix4x4()
        for i in range(0,4):
            for j in range(0,4):
                matrix_val = float(f.readline())
                matrix.SetElement(i,j, matrix_val)
        transform.SetMatrix(matrix)
        del matrix

    start_time = time.time()
    pd = read_polydata(in_filename)
    elapsed_time = time.time() - start_time
    print(f"READ: {elapsed_time}")
    # Transform it.
    start_time = time.time()
    transformer = vtk.vtkTransformPolyDataFilter()
    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        transformer.SetInputData(pd)
    else:
        transformer.SetInput(pd)
    transformer.SetTransform(transform)
    transformer.Update()
    elapsed_time = time.time() - start_time
    print(f"TXFORM: {elapsed_time}")

    # Write it out.
    start_time = time.time()
    pd2 = transformer.GetOutput()
    write_polydata(pd2, out_filename)
    elapsed_time = time.time() - start_time
    print(f"WRITE: {elapsed_time}")

    # Clean up.
    del transformer
    del pd2
    del pd
    del transform

def transform_polydata_from_disk_using_transform_object(in_filename, transform, out_filename):
    # Read it in.
    print(f"<{os.path.basename(__file__)}> Transforming {in_filename} -> {out_filename}...")
    start_time = time.time()
    pd = read_polydata(in_filename)
    elapsed_time = time.time() - start_time
    print(f"READ: {elapsed_time}")
    # Transform it.
    start_time = time.time()
    transformer = vtk.vtkTransformPolyDataFilter()
    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        transformer.SetInputData(pd)
    else:
        transformer.SetInput(pd)
    transformer.SetTransform(transform)
    transformer.Update()
    elapsed_time = time.time() - start_time
    print(f"TXFORM: {elapsed_time}")

    # Write it out.
    start_time = time.time()
    pd2 = transformer.GetOutput()
    write_polydata(pd2, out_filename)
    elapsed_time = time.time() - start_time
    print(f"WRITE: {elapsed_time}")

    # Clean up.
    del transformer
    del pd2
    del pd

def transform_polydatas_from_disk(input_dir, transforms, output_dir):
    """Loop over all input polydata files and apply the vtk transforms from the

    input transforms list. Save transformed polydata files in the output
    directory. As long as files were read in using list_vtk_files
    originally, they will be in the same order as the transforms now.
    """

    # Find input files
    input_pd_fnames = list_vtk_files(input_dir)
    num_pd = len(input_pd_fnames)
    print(f"<{os.path.basename(__file__)}> =======================================")
    print(f"<{os.path.basename(__file__)}> Transforming vtk and vtp files from directory: {input_dir}")
    print(f"<{os.path.basename(__file__)}> Total number of files found: {num_pd}")
    print(f"<{os.path.basename(__file__)}> Writing output to directory: {output_dir}")
    print(f"<{os.path.basename(__file__)}> =======================================")

    if not os.path.exists(output_dir):
        print(f"<{os.path.basename(__file__)}> ERROR: Output directory does not exist.")
        return
    if not os.path.exists(input_dir):
        print(f"<{os.path.basename(__file__)}> ERROR: Output directory does not exist.")
        return

    # Transform the files
    for idx in range(0, len(input_pd_fnames)):
        in_filename = input_pd_fnames[idx]
        subject_id = os.path.splitext(os.path.basename(in_filename))[0]
        out_filename = os.path.join(output_dir, f'{subject_id}_reg.vtk')
        transform_polydata_from_disk_using_transform_object(in_filename, transforms[idx], out_filename)

# This function was faster but is not always safe. Could crash due to missing file if any write/disk access issue
def transform_polydatas_from_diskUNSAFE(input_dir, transforms, output_dir, parallel_jobs=3):
    """Loop over all input polydata files and apply the vtk transforms from the

    input transforms list. Save transformed polydata files in the output
    directory. As long as files were read in using list_vtk_files
    originally, they will be in the same order as the transforms now.
    """

    # Find input files
    input_pd_fnames = list_vtk_files(input_dir)
    num_pd = len(input_pd_fnames)
    print(f"<{os.path.basename(__file__)}> =======================================")
    print(f"<{os.path.basename(__file__)}> Transforming vtk and vtp files from directory: {input_dir}")
    print(f"<{os.path.basename(__file__)}> Total number of files found: {num_pd}")
    print(f"<{os.path.basename(__file__)}> Writing output to directory: {output_dir}")
    print(f"<{os.path.basename(__file__)}> =======================================")

    if not os.path.exists(output_dir):
        print(f"<{os.path.basename(__file__)}> ERROR: Output directory does not exist.")
        return
    if not os.path.exists(input_dir):
        print(f"<{os.path.basename(__file__)}> ERROR: Output directory does not exist.")
        return

    # Set up inputs for subprocesses
    fname_list = list()
    out_fname_list = list()
    message_list = list()
    transform_list = list()
    for idx in range(0, len(input_pd_fnames)):
        fname = input_pd_fnames[idx]
        tx = transforms[idx]
        subject_id = os.path.splitext(os.path.basename(fname))[0]
        out_fname = os.path.join(output_dir, f'{subject_id}_reg.vtk')
        fname_list.append(fname)
        out_fname_list.append(out_fname)
        # save the transform to disk because we cannot pickle it
        if tx.GetClassName() == 'vtkThinPlateSplineTransform':
            writer = vtk.vtkMNITransformWriter()
            writer.AddTransform(tx)
            fname = f'.tmp_vtk_txform_{str(subject_id)}.xfm'
            fname = os.path.join(output_dir, fname)
            writer.SetFileName(fname)
            writer.Write()
            del writer
        elif tx.GetClassName() == 'vtkBSplineTransform':
            # there is no vtk file format for bspline transforms, unfortunately.
            print("Saving bspline transformed polydata without writing transform to disk for multiprocessing")
            transform_polydata_from_disk_using_transform_object(fname, tx, out_fname)
            # this did not work, and also it would lose the grid size information
            ## image = tx.GetCoefficients()
            ## writer = vtk.vtkImageWriter()
            ## writer.SetFileDimensionality(1)
            ## fname = '.tmp_vtk_txform_' + str(subject_id) + '.img'
            ## fname = os.path.join(output_dir, fname)
            ## writer.SetFileName(fname)
            ## writer.SetInput(image)
            ## writer.Write()
            ## del writer
        else:
            fname = f'.tmp_vtk_txform_{str(subject_id)}.txt'
            f = open(fname, 'w')
            for i in range(0,4):
                for j in range(0,4):
                    f.write(str(tx.GetMatrix().GetElement(i,j))+'\n')
            f.close()
        transform_list.append(fname)

    # Run this in parallel
    if tx.GetClassName() != 'vtkBSplineTransform':
        ret = Parallel(
                n_jobs=parallel_jobs, verbose=0)(
                    delayed(transform_polydata_from_disk)(in_filename, transform_filename, out_filename)
                    for (in_filename, transform_filename, out_filename) in zip(fname_list, transform_list, out_fname_list))

        # remove the temporary transform files
        for fname in transform_list:
            os.remove(fname)

    #for idx in range(0, len(input_pd_fnames)):
        #print(f"<{os.path.basename(__file__)}> {idx + 1} / {num_pd} {subject_id} Transforming {in_filename} -> {out_filename}...")
        #transform_polydata_from_disk(in_filename, transform, out_filename)

def transform_polydatas_from_diskOLD(input_dir, transforms, output_dir):
    """Loop over all input polydata files and apply the vtk transforms from the

    input transforms list. Save transformed polydata files in the output
    directory. As long as files were read in using list_vtk_files
    originally, they will be in the same order as the transforms now.
    """

    # Find input files
    input_pd_fnames = list_vtk_files(input_dir)
    num_pd = len(input_pd_fnames)
    print(f"<{os.path.basename(__file__)}> =======================================")
    print(f"<{os.path.basename(__file__)}> Transforming vtk and vtp files from directory: {input_dir}")
    print(f"<{os.path.basename(__file__)}> Total number of files found: {num_pd}")
    print(f"<{os.path.basename(__file__)}> Writing output to directory: {output_dir}")
    print(f"<{os.path.basename(__file__)}> =======================================")

    if not os.path.exists(output_dir):
        print(f"<{os.path.basename(__file__)}> ERROR: Output directory does not exist.")
        return
    if not os.path.exists(input_dir):
        print(f"<{os.path.basename(__file__)}> ERROR: Output directory does not exist.")
        return

    for idx in range(0, len(input_pd_fnames)):
        # Read it in.
        fname = input_pd_fnames[idx]
        subject_id = os.path.splitext(os.path.basename(fname))[0]
        out_fname = os.path.join(output_dir, f'{subject_id}_reg.vtk')
        print(f"<{os.path.basename(__file__)}> {idx + 1} / {num_pd} {subject_id} Transforming {fname} -> {out_fname}...")
        pd = read_polydata(fname)
        # Transform it.
        transformer = vtk.vtkTransformPolyDataFilter()
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            transformer.SetInputData(pd)
        else:
            transformer.SetInput(pd)
        transformer.SetTransform(transforms[idx])
        transformer.Update()
        # Write it out.
        pd2 = transformer.GetOutput()
        write_polydata(pd2, out_fname)
        # Clean up.
        del transformer
        del pd2
        del pd

def write_transforms_to_itk_format(transform_list, outdir, subject_ids=None):
    """Write VTK affine or spline transforms to ITK 4 text file formats.

    Input transforms are in VTK RAS space and are forward transforms. Output
    transforms are in LPS space and are the corresponsing inverse
    transforms, according to the conventions for these file formats and for
    resampling images. The affine transform is straightforward. The spline
    transform file format is just a list of displacements that have to be in
    the same order as they are stored in ITK C code. This now outputs an ITK
    transform that works correctly to transform the tracts (or any volume in
    the same space) in Slicer. In the nonrigid case, we also output a vtk
    native spline transform file using MNI format.
    """

    idx = 0
    tx_fnames = list()
    for tx in transform_list:

        # save out the vtk transform to a text file as it is
        # The MNI transform reader/writer are available in vtk so use those:
        if tx.GetClassName() != 'vtkBSplineTransform':
            writer = vtk.vtkMNITransformWriter()
            writer.AddTransform(tx)
            if subject_ids is not None:
                fname = f'vtk_txform_{str(subject_ids[idx])}.xfm'
            else:
                fname = f'vtk_txform_{idx:05d}.xfm'
            writer.SetFileName(os.path.join(outdir, fname))
            writer.Write()

        # file name for itk transform written below
        if subject_ids is not None:
            fname = f'itk_txform_{str(subject_ids[idx])}.tfm'
        else:
            fname = f'itk_txform_{idx:05d}.tfm'
        fname = os.path.join(outdir, fname)
        tx_fnames.append(fname)

        # Save the itk transform as the inverse of this transform (resampling transform) and in LPS.
        # This will show the same transform in the slicer GUI as the vtk transform we internally computed
        # that is stored in the .xfm text file, above.
        # To apply our transform to resample a volume in LPS:
        # convert to RAS, use inverse of transform to resample, convert back to LPS
        if tx.GetClassName() == 'vtkThinPlateSplineTransform' or  tx.GetClassName() == 'vtkBSplineTransform':
            #print 'Saving nonrigid transform displacements in ITK format'

            # Deep copy to avoid modifying input transform that will be applied to polydata
            if tx.GetClassName() == 'vtkThinPlateSplineTransform':
                tps = vtk.vtkThinPlateSplineTransform()
            else:
                tps = vtk.vtkBSplineTransform()
            tps.DeepCopy(tx)

            #extent = tps.GetCoefficients().GetExtent()
            #origin = tps.GetCoefficients().GetOrigin()
            #spacing = tps.GetCoefficients().GetSpacing()
            #dims = tps.GetCoefficients().GetDimensions()
            #print "E:", extent
            #print "O:", origin
            #print "S:", spacing
            #print "D:", dims

            # invert to get the transform suitable for resampling an image
            tps.Inverse()

            # convert the inverse spline transform from RAS to LPS
            ras_2_lps = vtk.vtkTransform()
            ras_2_lps.Scale(-1, -1, 1)
            lps_2_ras = vtk.vtkTransform()
            lps_2_ras.Scale(-1, -1, 1)
            spline_inverse_lps = vtk.vtkGeneralTransform()
            spline_inverse_lps.Concatenate(lps_2_ras)
            spline_inverse_lps.Concatenate(tps)
            spline_inverse_lps.Concatenate(ras_2_lps)

            # Now, loop through LPS space. Find the effect of the
            # inverse transform on each point. This is essentially what
            # vtk.vtkTransformToGrid() does, but this puts things into
            # LPS.

            # This low-res grid produced small differences (order of 1-2mm) when transforming
            # polydatas inside Slicer vs. in this code. 
            #grid_size = [15, 15, 15]
            #grid_spacing = 10
            # This higher-res grid has fewer small numerical differences
            # grid_size = [50, 50, 50]
            # grid_spacing = 5
            # This higher-res grid has fewer small numerical differences, but files are larger
            #grid_size = [70, 70, 70]
            #grid_spacing = 3

            # This higher-res grid is sufficient to limit numerical
            # differences to under .1mm in tests.  However, files are
            # quite large (47M). As this is still much smaller than
            # the tractography files, and correctness is desired, we
            # will produce large transform files. A preferable
            # solution would be to store the forward transform we
            # compute at the grid points at which it is defined, but
            # there is no inverse flag available in the file
            # format. Therefore the inverse must be stored at high
            # resolution.
            grid_size = [105, 105, 105]
            grid_spacing = 2

            extent_0 = [-int((grid_size[0] - 1)/2), -int((grid_size[1] - 1)/2), -int((grid_size[2] - 1)/2)]
            extent_1 = [ int((grid_size[0] - 1)/2),  int((grid_size[1] - 1)/2),  int((grid_size[2] - 1)/2)]

            origin = -grid_spacing * (np.array(extent_1) - np.array(extent_0))/2.0

            grid_points_LPS = list()
            grid_points_RAS = list()

            # ordering of grid points must match itk-style array order for images
            for s in range(extent_0[0], extent_1[0]+1):
                for p in range(extent_0[1], extent_1[1]+1):
                    for l in range(extent_0[2], extent_1[2]+1):
                        grid_points_RAS.append([-l*grid_spacing, -p*grid_spacing, s*grid_spacing])
                        grid_points_LPS.append([l*grid_spacing, p*grid_spacing, s*grid_spacing])

            displacements_LPS = list()

            print(f"LPS grid for storing transform: {grid_points_LPS[0]} {grid_points_LPS[-1]} {grid_spacing}")

            lps_points = vtk.vtkPoints()
            lps_points2 = vtk.vtkPoints()
            for gp_lps in grid_points_LPS:
                lps_points.InsertNextPoint(gp_lps[0], gp_lps[1], gp_lps[2])

            spline_inverse_lps.TransformPoints(lps_points, lps_points2)
            pidx = 0
            for gp_lps in grid_points_LPS:
                pt = lps_points2.GetPoint(pidx)
                diff_lps = [pt[0] - gp_lps[0], pt[1] - gp_lps[1], pt[2] - gp_lps[2]]
                pidx += 1

                ## # this tested grid definition and origin were okay.
                ## diff_lps = [20,30,40]

                ## # this tested that the ordering of L,P,S is correct:
                ## diff_lps = [0, gp_lps[1], 0]
                ## diff_lps = [gp_lps[0], 0, 0]
                ## diff_lps = [0, 0, gp_lps[2]]

                ## # this tested that the ordering of grid points is correct
                ## # only the R>0, A>0, S<0 region shows a transform.
                ## if gp_lps[0] < 0 and gp_lps[1] < 0 and gp_lps[2] < 0:
                ##     diff_lps = [gp_lps[0]/2.0, 0, 0]
                ## else:
                ##     diff_lps = [0, 0, 0]

                displacements_LPS.append(diff_lps)

            # save the points and displacement vectors in ITK format.
            #print 'Saving in ITK transform format.'
            f = open(fname, 'w')
            f.write('#Insight Transform File V1.0\n')
            f.write('# Transform 0\n')
            # ITK version 3 that included an additive (!) affine transform
            #f.write('Transform: BSplineDeformableTransform_double_3_3\n')
            # ITK version 4 that does not include a second transform in the file
            f.write('Transform: BSplineTransform_double_3_3\n')
            f.write('Parameters: ')
            # "Here the data are: The bulk of the BSpline part are 3D
            # displacement vectors for each of the BSpline grid-nodes
            # in physical space, i.e. for each grid-node, there will
            # be three blocks of displacements defining dx,dy,dz for
            # all grid nodes."
            for block in [0, 1, 2]:
                for diff in displacements_LPS:
                    f.write(f'{diff[block]} ')

            #FixedParameters: size size size origin origin origin origin spacing spacing spacing (then direction cosines: 1 0 0 0 1 0 0 0 1)
            f.write('\nFixedParameters:')
            #f.write(' {0} {0} {0}'.format(2*sz+1))
            f.write(f' {grid_size[0]}')
            f.write(f' {grid_size[1]}')
            f.write(f' {grid_size[2]}')

            f.write(f' {origin[0]}')
            f.write(f' {origin[1]}')
            f.write(f' {origin[2]}')
            f.write(f' {grid_spacing} {grid_spacing} {grid_spacing}')
            f.write(' 1 0 0 0 1 0 0 0 1\n')

            f.close()
        else:
            tx_inverse = vtk.vtkTransform()
            tx_inverse.DeepCopy(tx)
            tx_inverse.Inverse()
            ras_2_lps = vtk.vtkTransform()
            ras_2_lps.Scale(-1, -1, 1)
            lps_2_ras = vtk.vtkTransform()
            lps_2_ras.Scale(-1, -1, 1)
            tx2 = vtk.vtkTransform()
            tx2.Concatenate(lps_2_ras)
            tx2.Concatenate(tx_inverse)
            tx2.Concatenate(ras_2_lps)

            three_by_three = list()
            translation = list()
            for i in range(0,3):
                for j in range(0,3):
                    three_by_three.append(tx2.GetMatrix().GetElement(i,j))
            translation.append(tx2.GetMatrix().GetElement(0,3))
            translation.append(tx2.GetMatrix().GetElement(1,3))
            translation.append(tx2.GetMatrix().GetElement(2,3))

            f = open(fname, 'w')
            f.write('#Insight Transform File V1.0\n')
            f.write('# Transform 0\n')
            f.write('Transform: AffineTransform_double_3_3\n')
            f.write('Parameters: ')
            for el in three_by_three:
                f.write(f'{el} ')
            for el in translation:
                f.write(f'{el} ')
            f.write('\nFixedParameters: 0 0 0\n')
            f.close()

        idx +=1
    return(tx_fnames)


class LateralityResults:

    """Results of laterality computation for a subject.

    This class defines the structure returned by ComputeWhiteMatterLaterality.
    I/O functions included in this class are read and write.

    """

    def __init__(self):
        # I/O parameters
        self.directory = ''
        # results data storage
        self.polydata = None
        self.laterality_index = None
        self.right_hem_distance = None
        self.left_hem_distance = None
        # computation parameters
        self.sigma = None
        self.points_per_fiber = None
        self.threshold = None
        self.left_hem_similarity = None
        self.right_hem_similarity = None
        self.hemisphere = None
        
    def write(self, dirname, savedist=False):
        """Write output laterality results for one subject."""
        print("a")
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        print("b")
        # output polydata
        writer = vtk.vtkPolyDataWriter()
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            writer.SetInputData(self.polydata)
        else:
            writer.SetInput(self.polydata)
        
        writer.SetFileName(os.path.join(dirname, 'tractography_with_LI.vtk'))
        writer.Write()
        print("c")
        # output LI and other data values to text file
        # first pickle everything for later python processing
        fid = open(os.path.join(dirname, 'pickle_laterality_index.txt'), 'w')
        pickle.dump(self.laterality_index, fid)
        fid.close()

        fid = open(os.path.join(dirname, 'pickle_left_hem_similarity.txt'), 'w')
        pickle.dump(self.left_hem_similarity, fid)
        fid.close()
        fid = open(os.path.join(dirname, 'pickle_right_hem_similarity.txt'), 'w')
        pickle.dump(self.right_hem_similarity, fid)
        fid.close()
        fid = open(os.path.join(dirname, 'pickle_hemisphere.txt'), 'w')
        pickle.dump(self.hemisphere, fid)
        fid.close()

        if savedist:
            fid = open(os.path.join(dirname, 'pickle_right_hem_distance.txt'), 'w')
            pickle.dump(self.right_hem_distance, fid)
            fid.close()
            fid = open(os.path.join(dirname, 'pickle_left_hem_distance.txt'), 'w')
            pickle.dump(self.left_hem_distance, fid)
            fid.close()
            print("d")
        # now output human-readable LI values
        fid = open(os.path.join(dirname, 'laterality_index_values.txt'), 'w')
        for idx in range(0, len(self.laterality_index)):
            fid.write(str(self.laterality_index[idx]))
            fid.write('\n')
        fid.close()
        print("e")
        # generate histogram (needs matplotlib)
        li_stats = self.laterality_index[np.nonzero(self.laterality_index)]
        ## if USE_MATPLOTLIB and 0:
        ##     try:
        ##         print "f"
        ##         matplotlib.pyplot.hist(li_stats, bins=25, color=[.6, .6, .6])
        ##         print "f1"
        ##         matplotlib.pyplot.savefig(
        ##             os.path.join(dirname, 'LI_histogram.pdf'))
        ##         print "f2"
        ##         matplotlib.pyplot.close()
        ##         print "g"
        ##     except Exception:
        ##         print(f"<{os.path.basename(__file__)}> matplotlib was unable to write histogram.")
        ##         raise
        ## print "1"
        # generate fiber visualization
        #try:
            # pd3_li = filter.maskFibers(pd3_out,abs(li)>max(abs(li))*.01, li)
            #ren = render.render(self.polydata)
            #ren.save_views(dirname)
        #except Exception:
        #    print(f"<{os.path.basename(__file__)}> vtk or rendering issue. Failed to save views.")
        #    print(f"<{os.path.basename(__file__)}> polydata was saved to disk so you can re-render.")
        #    raise

        #print "IMPLEMENT SAVING OF PARAMETERS TOO"

    def read(self, dirname, readpd=False, readdist=False):
        """Read output (class laterality.LateralityResults) for one subject."""

        if not os.path.isdir(dirname):
            print(f"<{os.path.basename(__file__)}> error: directory does not exist {dirname}")

        if readpd:
            # input polydata
            reader = vtk.vtkPolyDataReader()
            fname = os.path.join(dirname, 'tractography_with_LI.vtk')
            reader.SetFileName(fname)
            reader.Update()
            self.polydata = reader.GetOutput()

        # input LI and other data values using pickle
        fname = os.path.join(dirname, 'pickle_laterality_index.txt')
        fid = open(fname)
        self.laterality_index = pickle.load(fid)
        fid.close()

        fname = os.path.join(dirname, 'pickle_left_hem_similarity.txt')
        fid = open(fname)
        self.left_hem_similarity = pickle.load(fid)
        fid.close()
        fname = os.path.join(dirname, 'pickle_right_hem_similarity.txt')
        fid = open(fname)
        self.right_hem_similarity = pickle.load(fid)
        fid.close()
        fname = os.path.join(dirname, 'pickle_hemisphere.txt')
        fid = open(fname)
        self.hemisphere = pickle.load(fid)
        fid.close()

        if readdist:
            fid = open(os.path.join(dirname, 'pickle_right_hem_distance.txt'))
            self.right_hem_distance = pickle.load(fid)
            fid.close()
            fid = open(os.path.join(dirname, 'pickle_left_hem_distance.txt'))
            self.left_hem_distance = pickle.load(fid)
            fid.close()

        self.directory = dirname

        #print "IMPLEMENT READING OF PARAMETERS TOO"

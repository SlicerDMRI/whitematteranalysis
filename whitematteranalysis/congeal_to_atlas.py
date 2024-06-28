# -*- coding: utf-8 -*-

""" congeal_multisubject.py

re-implementation of fiber tractography registration (group)

class MultiSubjectRegistration


"""

import os

import numpy as np
import vtk

import whitematteranalysis as wma


class SubjectToAtlasRegistration:

    def __init__(self):
        # parameters that can be set by user
        self.output_directory = "."
        self.random_seed = 1000
        self.input_polydata_filename = None
        
        # performance options set by user
        self.verbose = False
        self.points_per_fiber = 15
        
        # optimizer parameters set by user
        self.maxfun = 300
        self.sigma = 20
        self.initial_step = 10
        self.final_step = 5
        self.mean_brain_size = 4000
        # A different sample of 1000 is taken every iteration, so we really
        # see more fibers than this over the course of the registration
        self.subject_brain_size = 1000
        # options are Affine, Rigid, Nonrigid
        self.mode = "Affine"

        self.fiber_length = 40
        self.fiber_length_max = 260

        # internal stuff
        self.subject_polydata = None
        self.atlas_polydata = None
        self.transform = None
        self.transform_as_array = None
        self.objectives = list()
        self.total_iterations = 0
        self.subject_id = None

        self.target_landmarks = list()
        self.nonrigid_grid_resolution = 6
                    
    def update_nonrigid_grid(self):
        """This updates the nonrigid grid. Subject data must be added first.
        
        """
        
        # create the bspline transform from this displacement field
        vtktrans = wma.register_two_subjects_nonrigid_bsplines.convert_transform_to_vtk(self.transform_as_array)
        # now compute a new displacement field from the transform, using the new grid resolution
        # This code always uses a grid of 200mm x 200mm x 200mm
        size_mm = 200.0
        dims = self.nonrigid_grid_resolution
        spacing = size_mm / (dims - 1)
        origin = -size_mm / 2.0
        grid = vtk.vtkTransformToGrid()
        grid.SetGridOrigin(origin, origin, origin)
        grid.SetGridSpacing(spacing, spacing, spacing)
        grid.SetGridExtent(0, dims-1, 0, dims-1, 0, dims-1)
        grid.SetGridScalarTypeToFloat()
        grid.SetInput(vtktrans)
        grid.Update()
        # convert the vtkImageData displacement field back to a numpy object
        displacement_field_vtk = grid.GetOutput()
        #print displacement_field_vtk.GetPointData().GetArray(0)
        newtrans = vtk.util.numpy_support.vtk_to_numpy(displacement_field_vtk.GetPointData().GetArray(0)).ravel()
        #print newtrans.shape
        print(f"UPDATE NONRIGID GRID: {self.nonrigid_grid_resolution} {len(self.transform_as_array)} ==> {len(newtrans)}", end=' ')
        self.transform_as_array = newtrans
        
    def set_subject(self, polydata, subject_id):
        self.subject_polydata = polydata
        if self.mode == "Nonrigid":
            # This sets up identity transform to initialize.
            res = self.nonrigid_grid_resolution
            trans = np.zeros(res*res*res*3)
            vtktrans = wma.register_two_subjects_nonrigid_bsplines.convert_transform_to_vtk(trans)
            self.transform = vtktrans
            self.transform_as_array = trans
            print(f"ADD PD: {trans}")
        else:
            trans = vtk.vtkTransform()
            self.transform = trans
            self.transform_as_array =  np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0]).astype(float)
        self.subject_id = subject_id

    def set_atlas(self, polydata, atlas_id):
        self.atlas_polydata = polydata
        self.atlas_id = atlas_id

    def iterate(self):
        self.total_iterations += 1

        # make a directory for the current iteration
        dirname = f"iteration_{self.total_iterations:05d}_sigma_{int(self.sigma):05d}"
        outdir = os.path.join(self.output_directory, dirname)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        print("filtering and downsampling atlas")
        mean_brain = wma.filter.preprocess(self.atlas_polydata, self.fiber_length, max_length_mm=self.fiber_length_max, return_indices=False, preserve_point_data=False, preserve_cell_data=False, verbose=False)
        mean_brain, tmp = wma.filter.downsample(mean_brain, self.mean_brain_size, return_indices=True, verbose=False, random_seed=self.random_seed+self.total_iterations) 

        mean_fibers = wma.fibers.FiberArray()
        mean_fibers.convert_from_polydata(mean_brain, self.points_per_fiber)
        fixed = np.array([mean_fibers.fiber_array_r,mean_fibers.fiber_array_a,mean_fibers.fiber_array_s])
        
        print("filtering and downsampling subject")
        subject_brain = wma.filter.preprocess(self.subject_polydata, self.fiber_length, max_length_mm=self.fiber_length_max, return_indices=False, preserve_point_data=False, preserve_cell_data=False, verbose=False)
        subject_brain, tmp = wma.filter.downsample(subject_brain, self.subject_brain_size, return_indices=True, verbose=False, random_seed=self.random_seed+self.total_iterations)

        fibers = wma.fibers.FiberArray()
        fibers.convert_from_polydata(subject_brain, self.points_per_fiber)
        moving = np.array([fibers.fiber_array_r,fibers.fiber_array_a,fibers.fiber_array_s])

        subject_idx = 1
        iteration_count = self.total_iterations
        output_directory = self.output_directory
        step_size = np.array([self.initial_step, self.final_step])
        render = False
        
        (self.transform_as_array, objectives, diff) = wma.congeal_multisubject.congeal_multisubject_inner_loop(fixed, moving, self.transform_as_array, self.mode, self.sigma, subject_idx, iteration_count, self.output_directory, step_size, self.maxfun, render, self.nonrigid_grid_resolution)

        if self.mode == "Nonrigid":
            vtktrans = wma.register_two_subjects_nonrigid_bsplines.convert_transform_to_vtk(self.transform_as_array)
            print(vtktrans)
        else:
            vtktrans = wma.register_two_subjects.convert_transform_to_vtk(self.transform_as_array)
            print(vtktrans.GetMatrix())
        self.transform = vtktrans

        # get the subject's objectives as computed so far
        # also get the total objective right now (future. now just printed to screen)
        # if requested, render all (future)

        # save the current transforms to disk
        tx_list = list()
        tx_list.append(self.transform)
        id_list = list()
        id_list.append(self.subject_id)
        wma.io.write_transforms_to_itk_format(tx_list, outdir, id_list)

        
    def save_transformed_polydata(self, intermediate_save=False):

        print("SAVING POLYDATA")

        if intermediate_save:
            # make a directory for the current iteration
            dirname = f"iteration_{self.total_iterations:05d}_sigma_{int(self.sigma):05d}"
            outdir = os.path.join(self.output_directory, dirname)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            
            # make a directory for output
            outdir_pds = os.path.join(outdir, 'transformed_output')
            if not os.path.exists(outdir_pds):
                os.makedirs(outdir_pds)
            self.save_transformed_polydata_to_disk(outdir_pds)
            
        else:
            # make a directory for the final output
            outdir = os.path.join(self.output_directory, 'output_tractography')
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            print("______________ FINAL__________________")
            if self.mode == "Nonrigid":
                print(self.transform)
            else:
                print(self.transform.GetMatrix())

            self.save_transformed_polydata_to_disk(outdir)

            # save the transform to text files
            tx_list = list()
            tx_list.append(self.transform)
            id_list = list()
            id_list.append(self.subject_id)
            wma.io.write_transforms_to_itk_format(tx_list, outdir, id_list)
 
    def save_transformed_polydata_to_disk(self, outdir):
        out_fname = os.path.join(outdir, f'{self.subject_id}_reg.vtk')
        print(f"{self.subject_id} Transforming {self.input_polydata_filename} -> {out_fname}...")

        pd = wma.io.read_polydata(self.input_polydata_filename)

        transformer = vtk.vtkTransformPolyDataFilter()
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            transformer.SetInputData(pd)
        else:
            transformer.SetInput(pd)
        transformer.SetTransform(self.transform)
        transformer.Update()

        pd2 = transformer.GetOutput()
        wma.io.write_polydata(pd2, out_fname)


""" congeal_multisubject.py

re-implementation of fiber tractography registration (group)

class MultiSubjectRegistration


"""

import os
import numpy
import vtk
import whitematteranalysis as wma

class SubjectToAtlasRegistration:

    def __init__(self):
        # parameters that can be set by user
        self.output_directory = "."
        self.random_seed = None
        self.input_polydata_filename = None
        
        # performance options set by user
        self.verbose = False
        self.points_per_fiber = 15
        self.nonlinear = False
        
        # optimizer parameters set by user
        self.maxfun = 300
        self.sigma = 20
        self.initial_step = 10
        self.final_step = 5
        self.mean_brain_size = 4000
        # A different sample of 1000 is taken every iteration, so we really
        # see more fibers than this over the course of the registration
        self.subject_brain_size = 1000
        self.smooth_mean_brain = False
        #self.smooth_mean_brain = True
        #print "TEST SMOOTHING MEAN BRAIN"
        
        # internal stuff
        self.subject_polydata = None
        self.atlas_polydata = None
        self.transform = None
        self.transform_as_array = None
        self.objectives = list()
        self.total_iterations = 0
        self.subject_id = None

        self.target_landmarks = list()
        landmarks = list()
        for r in [-120, -60, 0, 60, 120]:
            for a in [-120, -60, 0, 60, 120]:
                for s in [-120, -60, 0, 60, 120]:
                    #self.target_landmarks.extend([r, a, s])input
                    landmarks.append([r, a, s])
        # now shuffle the order of these points so that the optimizer doesn't
        # start working in one corner only every time
        # the order was computed as
        # order = numpy.random.permutation(range(0,125))
        # hold it constant for now so we don't have to pass it to the subprocesses.
        order = [ 75,  18,  54,  61,  64,  73,  95,  13, 111, 118,  43,   7,  46, 56,   4, 124,  77,  98,  72,  60,  38,  80,  36,  27, 120, 119, 51,  81,   0,  93,  11,  41,  69,  83, 107,  12, 106,  30,  53, 105,  33,  91,  28,  17,  58,  90,  45,  94,  14,  26,  84,   1, 92,  21,  47,  59, 100,   2,   3,  87,  65, 102,  68,  20,  85, 79,  82,  15,  32,  88, 115,  74,   6,  19,  35,  99, 104, 109, 70, 101,  96,  66,  52,  48,  49,  31,  97, 122,  78, 113,  55, 112,  76,  44,  23, 103,  16,  10, 123,  86,  39,   8,  62, 110, 42, 114,  40, 117,  63,   9,  25,  67,  71,  37,  24, 116,  57, 89, 121,  34,   5,  29, 108,  50,  22]
        for idx in order:
            self.target_landmarks.extend(landmarks[idx])
            #print landmarks[idx]

        #print "LANDMARKS FINAL"
        #print self.target_landmarks
        
        self.target_points = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(self.target_landmarks)
                    
    def set_subject(self, polydata, subject_id):
        self.subject_polydata = polydata
        if self.nonlinear:
            # set source and target points equal for initial identity transform
            trans = wma.register_two_subjects_nonlinear.compute_thin_plate_spline_transform(self.target_points,self.target_points)
            self.transform = trans
            self.transform_as_array = self.target_landmarks
            print "ADD PD:", trans
        else:
            trans = vtk.vtkTransform()
            self.transform = trans
            self.transform_as_array =  numpy.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0]).astype(float)
        self.subject_id = subject_id

    def set_atlas(self, polydata, atlas_id):
        self.atlas_polydata = polydata
        self.atlas_id = atlas_id

    def iterate(self):
        self.total_iterations += 1

        # make a directory for the current iteration
        dirname = "iteration_%05d_sigma_%05d" % (self.total_iterations, self.sigma)
        outdir = os.path.join(self.output_directory, dirname)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        mean_brain = wma.filter.downsample(self.atlas_polydata, self.mean_brain_size, verbose=False, random_seed=self.random_seed)
        mean_fibers = wma.fibers.FiberArray()
        mean_fibers.convert_from_polydata(mean_brain, self.points_per_fiber)
        fixed = numpy.array([mean_fibers.fiber_array_r,mean_fibers.fiber_array_a,mean_fibers.fiber_array_s])
        
        subject_brain = wma.filter.downsample(self.subject_polydata, self.subject_brain_size, verbose=False, random_seed=self.random_seed)
        fibers = wma.fibers.FiberArray()
        fibers.convert_from_polydata(subject_brain, self.points_per_fiber)
        moving = numpy.array([fibers.fiber_array_r,fibers.fiber_array_a,fibers.fiber_array_s])

        if self.nonlinear:
            mode = 'Nonlinear'
        else:
            mode = 'Linear'
        subject_idx = 1
        iteration_count = self.total_iterations
        output_directory = self.output_directory
        step_size = numpy.array([self.initial_step, self.final_step])
        render = False
        
        self.transform_as_array = wma.congeal_multisubject.congeal_multisubject_inner_loop(fixed, moving, self.transform_as_array, mode, self.sigma, subject_idx, iteration_count, self.output_directory, step_size, self.maxfun, render)

        if self.nonlinear:
            vtktrans = wma.register_two_subjects_nonlinear.convert_transform_to_vtk(self.transform_as_array, self.target_points)
            print vtktrans
        else:
            vtktrans = wma.register_two_subjects.convert_transform_to_vtk(self.transform_as_array)
            print vtktrans.GetMatrix()
        self.transform = vtktrans

        # get the subject's objectives as computed so far
        # also get the total objective right now (future. now just printed to screen)
        # if requested, render all (future)

        # save the current transforms to disk
        tx_list = list()
        tx_list.append(self.transform)
        id_list = list()
        id_list.append(self.subject_id)
        wma.registration_functions.write_transforms_to_itk_format(tx_list, outdir, id_list)

        
    def save_transformed_polydata(self, intermediate_save=False):

        print "SAVING POLYDATA"

        if intermediate_save:
            # make a directory for the current iteration
            dirname = "iteration_%05d_sigma_%05d" % (self.total_iterations, self.sigma)
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

            print "______________ FINAL__________________"
            if self.nonlinear:
                print self.transform
            else:
                print self.transform.GetMatrix()

            self.save_transformed_polydata_to_disk(outdir)

            # save the transform to text files
            tx_list = list()
            tx_list.append(self.transform)
            id_list = list()
            id_list.append(self.subject_id)
            wma.registration_functions.write_transforms_to_itk_format(tx_list, outdir, id_list)
 
    def save_transformed_polydata_to_disk(self, outdir):
        out_fname = os.path.join(outdir, self.subject_id + '_reg.vtk')
        print self.subject_id, " Transforming ", self.input_polydata_filename, "->", out_fname, "..."

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


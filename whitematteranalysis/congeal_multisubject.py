""" congeal_multisubject.py

re-implementation of fiber tractography registration (group)

class MultiSubjectRegistration


"""

import os
import numpy
import vtk
from joblib import Parallel, delayed
import whitematteranalysis as wma

class MultiSubjectRegistration:

    def __init__(self):
        # parameters that can be set by user
        self.output_directory = "."
        self.input_directory = None
        self.random_seed = None
        
        # performance options set by user
        self.verbose = False
        self.parallel_jobs = 4
        self.parallel_verbose = 0
        #self.points_per_fiber = 5
        #self.points_per_fiber = 9
        self.points_per_fiber = 15
        self.render = True
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
        
        # internal stuff
        self.polydatas = list()
        self.transforms = list()
        self.transforms_as_array = list()
        self.objectives = list()
        self.total_iterations = 0
        self.subject_ids = list()

        self.target_landmarks = list()
        self.nonlinear_grid_resolution = 5
        self.nonlinear_grid_5 = [-120, -60, 0, 60, 120]
        self.nonlinear_grid_6 = [-120, -60, -20, 20, 60, 120]
        # random order so that optimizer does not start in one corner every time
        # the order was computed as
        # numpy.random.permutation(range(0,125))
        # numpy.random.permutation(range(0,216))     
        self.grid_order_5 = [ 75,  18,  54,  61,  64,  73,  95,  13, 111, 118,  43,   7,  46, 56,   4, 124,  77,  98,  72,  60,  38,  80,  36,  27, 120, 119, 51,  81,   0,  93,  11,  41,  69,  83, 107,  12, 106,  30,  53, 105,  33,  91,  28,  17,  58,  90,  45,  94,  14,  26,  84,   1, 92,  21,  47,  59, 100,   2,   3,  87,  65, 102,  68,  20,  85, 79,  82,  15,  32,  88, 115,  74,   6,  19,  35,  99, 104, 109, 70, 101,  96,  66,  52,  48,  49,  31,  97, 122,  78, 113,  55, 112,  76,  44,  23, 103,  16,  10, 123,  86,  39,   8,  62, 110, 42, 114,  40, 117,  63,   9,  25,  67,  71,  37,  24, 116,  57, 89, 121,  34,   5,  29, 108,  50,  22]
        self.grid_order_6 = [165,  63, 129, 170, 148, 131,   1, 205, 181,  69,  35, 100,   6, 13,  82, 193, 159, 130,  54, 164, 174,  62, 108, 101, 169,  93, 112,  42, 110, 213, 107,  47,  45, 140, 138, 199, 125, 117,   8, 41, 215,  74, 143, 155, 144, 187,  26,  60, 139,  58,  97,  10, 113,  34, 116,  33, 202, 210,  27, 135, 152, 206, 128,  16, 177, 25,  67, 192, 147, 132, 160, 158, 161,  90, 102,  49,  21, 191, 32,  18,  81, 157, 114, 175,  94, 172, 207, 186, 167, 163, 196, 118,  28,  43, 133, 171, 211,  77,  56, 195, 173,  57,  96,  29, 64, 180,  89, 190, 115,  20,  52,  50,   4, 141,  98, 134, 109, 149, 176, 212,  11,   0, 146,  65,  91,  23,  53,  44, 123,  87, 24, 178, 184,  68, 124,  46,  76, 151, 127, 204, 154, 150, 106, 70,  37,  84,  17,  12, 189,   2,  92,  36,  71,  39,  30,  75, 179, 168,  73, 121,  86, 214, 188,  59, 209,  22,  19, 153, 162, 99, 182,  14,  48, 119, 203,  66,  61, 103, 208, 145,  79,  85, 142,  72, 126, 194, 104, 122, 198, 120, 200, 183, 201,   3,  78, 40,  83, 137,  31, 111,  15,  51,   9, 185,  55,  38, 156, 136, 7,  95,  80, 105, 166,  88, 197,   5]
        self.initialize_nonlinear_grid()

    def initialize_nonlinear_grid(self):
        self.target_landmarks = list()
        if self.nonlinear_grid_resolution == 5:
            grid = self.nonlinear_grid_5
            grid_order = self.grid_order_5
        elif self.nonlinear_grid_resolution == 6:
            grid = self.nonlinear_grid_6
            grid_order = self.grid_order_6
        else:
            print "<congeal_multisubject.py> Error: Unknown nonlinear grid mode:", self.nonlinear_grid_resolution
        tmp = list()
        for r in grid:
            for a in grid:
                for s in grid:
                    tmp.append([r, a, s])
        # now shuffle the order of these points to avoid biases
        for idx in grid_order:
            self.target_landmarks.extend(tmp[idx])
        self.target_points = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(self.target_landmarks)

    def add_polydata(self, polydata, subject_id):
        self.polydatas.append(polydata)
        if self.nonlinear:
            # set source and target points equal for initial identity transform
            trans = wma.register_two_subjects_nonlinear.compute_thin_plate_spline_transform(self.target_points,self.target_points)
            self.transforms.append(trans)
            self.transforms_as_array.append(self.target_landmarks)
            print "ADD PD:", trans
        else:
            trans = vtk.vtkTransform()
            self.transforms.append(trans)
            self.transforms_as_array.append(numpy.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0]).astype(float))
        self.subject_ids.append(subject_id)
        
    def remove_mean_from_transforms(self):
        """ Remove mean rotations and mean scaling and mean
         translations from transforms.  the mean rotation and
         translation should not affect the objective function anyway.
         A mean scaling will affect the objective, i.e. if the brain
         shrinks all distances become smaller and the similarity is
         higher. Definitely not the desired effect."""

        if self.nonlinear:
            # remove any average displacement of each source point.
            # this means the mean of the source points must equal the target point
            transforms_array = numpy.array(self.transforms_as_array)
            meansource = numpy.mean(transforms_array, 0)
            landmarks = numpy.array(self.target_landmarks)
            meandisp = meansource - landmarks
            print "MEAN DISPLACEMENT:", meandisp

            if self.verbose:
                print "<congeal.py> TRANSFORMS"
                print numpy.round(transforms_array * 100) / 100
                print "<congeal.py> Removing current (accumulated) mean transform before computing objective:"
                print meandisp

            for transform in self.transforms_as_array:
                transform[:] = transform - meandisp

            transforms_array = numpy.array(self.transforms_as_array)
            meansource = numpy.mean(transforms_array, 0)
            landmarks = numpy.array(self.target_landmarks)
            meandisp = meansource - landmarks
            print "MEAN DISPLACEMENT 2:", meandisp

            print "Zero-mean based on affine part"
            matrix_average = numpy.zeros((3,4))
            target_landmarks = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(self.target_landmarks)
            for trans in self.transforms_as_array:
                affine_part = vtk.vtkLandmarkTransform()
                source_landmarks = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(trans)
                affine_part.SetSourceLandmarks(source_landmarks)
                affine_part.SetTargetLandmarks(target_landmarks)
                affine_part.SetModeToAffine()
                affine_part.Update()
                #print affine_part.GetMatrix()
                for row in [0,1,2]:
                    # the fourth row must be 0, 0, 0, 1 for homogeneous transform so don't need to average it
                    for col in [0,1,2,3]:
                        matrix_average[row,col] += affine_part.GetMatrix().GetElement(row,col)
            matrix_average = numpy.divide(matrix_average, len(self.transforms_as_array))

            # remove this average from the transforms
            # want: txform*meantxform^-1 * sourcept' = targetpt
            # new source points:
            # sourcept' = meantxform * sourcept 
            
            matrix = vtk.vtkMatrix4x4()
            for row in [0,1,2]:
                    # the fourth row must be 0, 0, 0, 1 for homogeneous transform so don't need to average it
                    for col in [0,1,2,3]:
                        matrix.SetElement(row,col, matrix_average[row,col])
            print "MEAN TRANS", matrix
            mean_trans = vtk.vtkMatrixToLinearTransform()
            mean_trans.SetInput(matrix)
            
            # modify all the source points using the tps inverse and the mean transform we will remove
            new_source_pts = list()
            for trans in self.transforms_as_array:
                #source_landmarks = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(trans)
                #tps =  wma.register_two_subjects_nonlinear.compute_thin_plate_spline_transform(source_landmarks,target_landmarks)
                #tps.Inverse()
                trans2 = list()
                for pt in zip(trans[::3], trans[1::3], trans[2::3]):
                        #pt2 = tps.TransformPoint(pt[0], pt[1], pt[2])
                        #pt3 = mean_trans.TransformPoint(pt2[0], pt2[1], pt2[2])
                        pt2 = mean_trans.TransformPoint(pt[0], pt[1], pt[2])
                        trans2.extend(pt2)
                new_source_pts.append(numpy.array(trans2))

            for (trans, newtrans) in zip(self.transforms_as_array, new_source_pts):
                #print "REMOVE MEAN FROM TXFORMS:", trans, "NEW:", newtrans
                print "Mean removed. Source points changed by:", numpy.min(trans-newtrans), numpy.max(trans-newtrans)

            self.transforms_as_array = new_source_pts

            # TEST ONLY
            target_landmarks = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(self.target_landmarks)
            matrix_average = numpy.zeros((3,4))
            for trans in self.transforms_as_array:
                affine_part = vtk.vtkLandmarkTransform()
                source_landmarks = wma.register_two_subjects_nonlinear.convert_numpy_array_to_vtk_points(trans)
                affine_part.SetSourceLandmarks(source_landmarks)
                affine_part.SetTargetLandmarks(target_landmarks)
                affine_part.SetModeToAffine()
                affine_part.Update()
                #print affine_part.GetMatrix()
                for row in [0,1,2]:
                    # the fourth row must be 0, 0, 0, 1 for homogeneous transform so don't need to average it
                    for col in [0,1,2,3]:
                        matrix_average[row,col] += affine_part.GetMatrix().GetElement(row,col)
            matrix_average = numpy.divide(matrix_average, len(self.transforms_as_array))

            # remove this average from the transforms
            # want: txform'*meantxform^-1 * sourcept = targetpt
            # new source points:
            # txform'*sourcept' = targetpt
            # sourcept' = meantxform^-1 * sourcept 
            
            matrix = vtk.vtkMatrix4x4()
            for row in [0,1,2]:
                    # the fourth row must be 0, 0, 0, 1 for homogeneous transform so don't need to average it
                    for col in [0,1,2,3]:
                        matrix.SetElement(row,col, matrix_average[row,col])
            print "MEAN TRANS AFTER CORRECTION:", matrix
            
        else:
            transforms_array = numpy.array(self.transforms_as_array)
            meantrans = numpy.mean(transforms_array, 0)
            if self.verbose:
                print "<congeal.py> TRANSFORMS"
                print numpy.round(transforms_array * 100) / 100
                print "<congeal.py> Removing current (accumulated) mean transform before computing objective:"
                print numpy.round(meantrans * 1000) / 1000        

            for transform in self.transforms_as_array:
                transform[0:6] = transform[0:6] - meantrans[0:6]
                transform[6:9] = transform[6:9] / meantrans[6:9]
                transform[9:15] = transform[9:15] - meantrans[9:15]

    def iterate(self):
        self.total_iterations += 1

        # make a directory for the current iteration
        dirname = "iteration_%05d_sigma_%05d" % (self.total_iterations, self.sigma)
        outdir = os.path.join(self.output_directory, dirname)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        # make a directory for any intermediate rendering by subprocesses
        outdir_render = os.path.join(outdir, 'progress_rendered')
        if not os.path.exists(outdir_render):
            os.makedirs(outdir_render)

        fibers_per_subject = self.mean_brain_size / (len(self.polydatas) - 1)
        print "Fibers per subject for computing mean brain:", fibers_per_subject, "=", self.mean_brain_size, "/",  len(self.polydatas) -1

        mean_list = list()
        subject_list = list()
        mode_list = list()
        sigma_list = list()
        subj_idx_list = list()
        iteration_list = list()
        outdir_list = list()
        stepsize_list = list()
        maxfun_list = list()
        render_list = list()
        grid_resolution_list = list()
        
        # register each subject to the current mean

        # Sample fibers from each subject for use in the "mean brain"
        subject_sampled_fibers = list()
        for (input_pd, trans) in zip(self.polydatas, self.transforms):
            # apply the current transform to this polydata for computation of mean brain
            transformer = vtk.vtkTransformPolyDataFilter()
            if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                transformer.SetInputData(input_pd)
            else:
                transformer.SetInput(input_pd)
            transformer.SetTransform(trans)
            transformer.Update()
            pd = wma.filter.downsample(transformer.GetOutput(), fibers_per_subject, verbose=False, random_seed=self.random_seed)
            subject_sampled_fibers.append(pd)
            del transformer
            
        # Loop over all subjects and prepare lists of inputs for subprocesses
        subj_idx = 0
        for input_pd in self.polydatas:
            # compute mean in a leave-one out fashion. Otherwise, the
            # optimal transform may be identity.
            appender = vtk.vtkAppendPolyData()
            for subj_idx2 in range(len(subject_sampled_fibers)):
                if subj_idx2 != subj_idx:
                    pd = subject_sampled_fibers[subj_idx2]
                    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                        appender.AddInputData(pd)
                    else:
                        appender.AddInput(pd)

            appender.Update()
            mean_brain = appender.GetOutput()
            mean_fibers = wma.fibers.FiberArray()
            mean_fibers.convert_from_polydata(mean_brain, self.points_per_fiber)
            mean_fibers = numpy.array([mean_fibers.fiber_array_r,mean_fibers.fiber_array_a,mean_fibers.fiber_array_s])
            #  R,A,S is the first index
            # then fiber number
            # then points along fiber
            mean_list.append(mean_fibers)

            pd = wma.filter.downsample(input_pd, self.subject_brain_size, verbose=False, random_seed=self.random_seed)
            fibers = wma.fibers.FiberArray()
            fibers.convert_from_polydata(pd, self.points_per_fiber)
            fibers_array = numpy.array([fibers.fiber_array_r,fibers.fiber_array_a,fibers.fiber_array_s])
            subject_list.append(fibers_array)
            sigma_list.append(self.sigma)
            if self.nonlinear:
                mode_list.append('Nonlinear')
            else:
                mode_list.append('Linear')                    
            subj_idx_list.append(subj_idx)
            subj_idx += 1
            iteration_list.append(self.total_iterations)
            outdir_list.append(outdir_render)
            stepsize_list.append(numpy.array([self.initial_step, self.final_step]))
            maxfun_list.append(self.maxfun)
            render_list.append(self.render)
            grid_resolution_list.append(self.nonlinear_grid_resolution)
            
        # multiprocess over subjects
        print "STARTING MULTIPROCESSING. NUMBER OF JOBS:", self.parallel_jobs
        # note we can't pass vtk objects to subprocesses since they can't be pickled.
        ret = Parallel(
            n_jobs=self.parallel_jobs, verbose=self.parallel_verbose)(
                delayed(congeal_multisubject_inner_loop)(fixed, moving, initial_transform, mode, sigma, subject_idx, iteration_count, output_directory, step_size, maxfun, render, grid_resolution)
                for (fixed, moving, initial_transform, mode, sigma, subject_idx, iteration_count, output_directory, step_size, maxfun, render, grid_resolution) in zip(mean_list, subject_list, self.transforms_as_array, mode_list, sigma_list, subj_idx_list, iteration_list, outdir_list, stepsize_list, maxfun_list, render_list, grid_resolution_list))

            
        #print "RETURNED VALUES", ret
        
        # get the current transform for each subject
        self.transforms_as_array = list()
        for trans in ret:
            self.transforms_as_array.append(trans)
            #print "APPEND TRANS OUTPUT:", trans
            
        # remove_mean_from_transforms
        self.remove_mean_from_transforms()

        # update our transforms list for the next iteration
        self.transforms = list()
        for trans in self.transforms_as_array:
            if self.nonlinear:
                vtktrans = wma.register_two_subjects_nonlinear.convert_transform_to_vtk(trans, self.target_points)
                print vtktrans
            else:
                vtktrans = wma.register_two_subjects.convert_transform_to_vtk(trans)
                print vtktrans.GetMatrix()
            self.transforms.append(vtktrans)

        # get the subject's objectives as computed so far
        # also get the total objective right now (future. now just printed to screen)
        # if requested, render all (future)

        # save the current transforms to disk
        wma.registration_functions.write_transforms_to_itk_format(self.transforms, outdir, self.subject_ids)

        
    def save_transformed_polydatas(self, intermediate_save=False, midsag_symmetric=False):

        print "SAVING POLYDATAS"
        
        transform_list = self.transforms
        subject_id_list = self.subject_ids
        if midsag_symmetric:
            transform_list = transform_list[::2]
            subject_id_list = subject_id_list[::2]

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
                
            output_pds = wma.registration_functions.transform_polydatas_from_disk(self.input_directory, transform_list, outdir_pds)

        else:
            # make a directory for the final output
            outdir = os.path.join(self.output_directory, 'output_tractography')
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            print "______________ FINAL__________________"
            for trans in  self.transforms:
                if self.nonlinear:
                    print trans
                else:
                    print trans.GetMatrix()

            output_pds = wma.registration_functions.transform_polydatas_from_disk(self.input_directory, transform_list, outdir)
        
            # save the current atlas representation to disk
            # right now this is all the input fibers from all subjects
            # we could downsample to mean brain size here.
            #wma.registration_functions.save_atlas(output_pds, outdir_current)
            output_pds = list()
            for (pd, trans) in zip (self.polydatas, self.transforms):
                transformer = vtk.vtkTransformPolyDataFilter()
                if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                    transformer.SetInputData(pd)
                else:
                    transformer.SetInput(pd)
                transformer.SetTransform(trans)
                transformer.Update()
                output_pds.append(transformer.GetOutput())

            appender = vtk.vtkAppendPolyData()
            for pd in output_pds:
                if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                    appender.AddInputData(pd)
                else:
                    appender.AddInput(pd)
            appender.Update()
            wma.io.write_polydata(appender.GetOutput(), os.path.join(self.output_directory, 'registration_atlas.vtk'))
            del appender

            # save the transforms to text files
            wma.registration_functions.write_transforms_to_itk_format(transform_list, outdir, subject_id_list)
            
def congeal_multisubject_inner_loop(mean, subject, initial_transform, mode, sigma, subject_idx, iteration_count, output_directory, step_size, maxfun, render, grid_resolution):

    print "\n BEGIN ITERATION", iteration_count, "subject", subject_idx, "sigma:", sigma, "mean brain:", mean.shape, "subject:", subject.shape, "initial transform length:", len(initial_transform), "steps:", step_size[0], step_size[1], "maxfun:", maxfun, type(initial_transform), "Grid:", grid_resolution, "Mode:", mode, "initial transform:", initial_transform,
    
    
    if mode == 'Linear':
        register = wma.register_two_subjects.RegisterTractography()
    elif mode == "Nonlinear":
        register = wma.register_two_subjects_nonlinear.RegisterTractographyNonlinear()
    else:
        print "ERROR: Unknown registration mode"
        
    register.maxfun = maxfun
    register.sigma = sigma        
    register.parallel_jobs = 1
    register.threshold = 0
    #register.distance_method = "Hausdorff"
    register.fixed = mean
    register.moving = subject
    register.initial_transform = initial_transform
    register.verbose = 0
    #register.verbose = 1
    register.output_directory = output_directory
    register.process_id_string = "_subject_%05d_iteration_%05d_sigma_%05d" % (subject_idx, iteration_count, sigma) 
    register.initial_step = step_size[0]
    register.final_step = step_size[1]
    register.render = render
    register.nonlinear_grid_resolution = grid_resolution
    register.initialize_nonlinear_grid()
    
    register.compute()

    # only update the transform if the objective function improved.
    # sometimes it falls into a different minimum that is worse
    obj_diff = register.objective_function_values[-1] - register.objective_function_values[0]
    print "\n END ITERATION", iteration_count, "subject", subject_idx, "OBJECTIVE CHANGE:", obj_diff
    if obj_diff < 0:
        print "UPDATING MATRIX"
        return register.final_transform
    else:
        return initial_transform


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
        
        # internal stuff
        self.polydatas = list()
        self.transforms = list()
        self.transforms_as_array = list()
        self.objectives = list()
        self.total_iterations = 0
        self.subject_ids = list()

    def add_polydata(self, polydata, subject_id):
        self.polydatas.append(polydata)
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
        outdir = os.path.join(self.output_directory, 'iteration_'+str(self.total_iterations))
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
        reg_list = list()
        sigma_list = list()
        subj_idx_list = list()
        iteration_list = list()
        outdir_list = list()
        stepsize_list = list()
        maxfun_list = list()
        render_list = list()
                
        # register each subject to the current mean
        subj_idx = 0
        for input_pd in self.polydatas:
            # compute mean in a leave-one out fashion. Otherwise, the
            # optimal transform may be identity.
            appender = vtk.vtkAppendPolyData()
            subj_idx2 = 0
            for (input_pd2, trans) in zip(self.polydatas, self.transforms):
                if input_pd2 != input_pd:
                    # apply the current transform to this polydata for computation of mean brain
                    transformer = vtk.vtkTransformPolyDataFilter()
                    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                        transformer.SetInputData(input_pd2)
                    else:
                        transformer.SetInput(input_pd2)
                    transformer.SetTransform(trans)
                    transformer.Update()
                    pd = wma.filter.downsample(transformer.GetOutput(), fibers_per_subject, verbose=False, random_seed=self.random_seed)
                    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                        appender.AddInputData(pd)
                    else:
                        appender.AddInput(pd)
                        
            appender.Update()
            mean_brain = appender.GetOutput()
            # this smoothing does seem to help. need to test more.
            if self.smooth_mean_brain:
                (mean_brain, weights) = wma.filter.smooth(mean_brain, fiber_distance_sigma = 20, points_per_fiber=30, n_jobs=self.parallel_jobs, upper_thresh=20)
                # here we could remove fibers with low weights
            
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
            subj_idx_list.append(subj_idx)
            subj_idx += 1
            iteration_list.append(self.total_iterations)
            outdir_list.append(outdir_render)
            stepsize_list.append(numpy.array([self.initial_step, self.final_step]))
            maxfun_list.append(self.maxfun)
            render_list.append(self.render)
            
        # multiprocess over subjects
        print "STARTING MULTIPROCESSING. NUMBER OF JOBS:", self.parallel_jobs
        
        ret = Parallel(
            n_jobs=self.parallel_jobs, verbose=self.parallel_verbose)(
                delayed(congeal_multisubject_inner_loop)(fixed, moving, initial_transform, sigma, subject_idx, iteration_count, output_directory, step_size, maxfun, render)
                for (fixed, moving, initial_transform, sigma, subject_idx, iteration_count, output_directory, step_size, maxfun, render) in zip(mean_list, subject_list, self.transforms_as_array, sigma_list, subj_idx_list, iteration_list, outdir_list, stepsize_list, maxfun_list, render_list))

        # get the current transform for each subject
        self.transforms_as_array = list()
        for trans in ret:
            self.transforms_as_array.append(trans)

        # remove_mean_from_transforms
        self.remove_mean_from_transforms()

        # update our transforms list for the next iteration
        self.transforms = list()
        for trans in self.transforms_as_array:
            vtktrans = wma.register_two_subjects.convert_transform_to_vtk(trans)
            print vtktrans.GetMatrix()
            self.transforms.append(vtktrans)

        # get the subject's objectives as computed so far
        # also get the total objective right now (future. now just printed to screen)
        # if requested, render all (future)

        # save the current transforms to disk
        wma.registration_functions.write_transforms_to_itk_format(self.transforms, outdir, self.subject_ids)

        
    def save_transformed_polydatas(self, intermediate_save=False, midsag_symmetric=False):

        transform_list = self.transforms
        if midsag_symmetric:
            transform_list = transform_list[::2]

        if intermediate_save:
            # make a directory for the current iteration
            outdir = os.path.join(self.output_directory, 'iteration_'+str(self.total_iterations))
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

            wma.registration_functions.write_transforms_to_itk_format(self.transforms, outdir, self.subject_ids)
            
            print "______________ FINAL__________________"
            for trans in  self.transforms:
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
    

def congeal_multisubject_inner_loop(mean, subject, initial_transform, sigma, subject_idx, iteration_count, output_directory, step_size, maxfun, render):

    print "ITERATION", iteration_count, "subject", subject_idx, "sigma:", sigma, "mean brain:", mean.shape, "subject:", subject.shape, "initial transform:", initial_transform, "steps:", step_size[0], step_size[1], "maxfun:", maxfun
    
    scale_mode = "Medium"
    register = wma.register_two_subjects.RegisterTractography()
    register.maxfun = maxfun
    register.sigma = sigma        
    register.parallel_jobs = 1
    register.threshold = 0
    #register.distance_method = "Hausdorff"
    register.fixed = mean
    register.moving = subject
    register.initial_transform = initial_transform
    register.verbose = 0
    register.output_directory = output_directory
    register.process_id_string = "_subject_%05d_iteration_%05d_sigma_%05d" % (subject_idx, iteration_count, sigma) 
    register.rhobeg = step_size[0]
    register.rhoend = step_size[1]
    register.render = render
    
    register.compute()

    # only update the transform if the objective function improved.
    # sometimes it falls into a different minimum that is worse
    obj_diff = register.objective_function_values[-1] - register.objective_function_values[0]
    print "ITERATION", iteration_count, "subject", subject_idx, "OBJECTIVE CHANGE:", obj_diff
    if obj_diff < 0:
        print "UPDATING MATRIX"
        return register.final_transform
    else:
        return initial_transform
    




## def transform_fiber_array(in_array, transform):
##     """Transform in_array (of class FiberArray) by transform (9
##     components, rotation about R,A,S, translation in R, A, S, and
##     scale along R, A, S. Fibers are assumed to be in RAS.
##     Transformed fibers are returned. """

##     out_array = whitematteranalysis.fibers.FiberArray()
##     out_array.number_of_fibers = in_array.number_of_fibers
##     out_array.points_per_fiber = in_array.points_per_fiber
##     # allocate array number of lines by line length
##     out_array.fiber_array_r = numpy.zeros((in_array.number_of_fibers,
##                                             in_array.points_per_fiber))
##     out_array.fiber_array_a = numpy.zeros((in_array.number_of_fibers,
##                                             in_array.points_per_fiber))
##     out_array.fiber_array_s = numpy.zeros((in_array.number_of_fibers,
##                                             in_array.points_per_fiber))

##     vtktrans = convert_transform_to_vtk(transform)

##     # Transform moving fiber array by applying transform to original fibers
##     for lidx in range(0, in_array.number_of_fibers):
##         for pidx in range(0, in_array.points_per_fiber):
##             pt = vtktrans.TransformPoint(in_array.fiber_array_r[lidx, pidx],
##                                          in_array.fiber_array_a[lidx, pidx], 
##                                          in_array.fiber_array_s[lidx, pidx])
##             out_array.fiber_array_r[lidx, pidx] = pt[0]
##             out_array.fiber_array_a[lidx, pidx] = pt[1]
##             out_array.fiber_array_s[lidx, pidx] = pt[2]

##     del vtktrans
##     return out_array

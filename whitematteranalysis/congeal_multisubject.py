""" congeal_multisubject.py

re-implementation of fiber tractography registration (group)

class MultiSubjectRegistration


"""

import os
import time
import numpy
import vtk
from joblib import Parallel, delayed

HAVE_PLT = 1
try:
    import matplotlib
    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except:
    print "<wm_congeal_multisubject.py> Error importing matplotlib.pyplot package, can't plot objectives.\n"
    HAVE_PLT = 0

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
        self.nonrigid = False
        
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
        self.objectives_before = list()
        self.objectives_after = list()
        self.total_iterations = 0
        self.subject_ids = list()

        #self.nonrigid_grid_resolution = 3
        #self.nonrigid_grid_resolution = 5
        self.nonrigid_grid_resolution = 6

    def update_nonrigid_grid(self):
        """This updates the nonrigid grids. Subjects must be added first using

        add_polydata.
        """
        
        # Apply the existing transform to the new grid
        new_transforms = list()
        for trans in self.transforms_as_array:
            # create the bspline transform from this displacement field
            vtktrans = wma.register_two_subjects_nonrigid_bsplines.convert_transform_to_vtk(trans)
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
            newtrans = vtk.util.numpy_support.vtk_to_numpy(displacement_field_vtk.GetPointData().GetArray(0))
            #print newtrans.shape
            new_transforms.append(newtrans.ravel())

        # Update all the relevant variables (the spline transform does not change but all source and target points do)
        print "UPDATE NONRIGID GRID: ", self.nonrigid_grid_resolution, len(trans), "==>", len(new_transforms[-1]),
        self.transforms_as_array = new_transforms

    def add_polydata(self, polydata, subject_id):
        """Add a subject's data to the groupwise registration. self.nonrigid

        must be set before calling this, if nonrigid registration is desired.
        """
        
        self.polydatas.append(polydata)
        if self.nonrigid:
            # This sets up identity transform to initialize. This will
            # be re-calculated with current grid resolution in
            # update_nonrigid_grid.
            res = self.nonrigid_grid_resolution
            trans = numpy.zeros(res*res*res*3)
            vtktrans = wma.register_two_subjects_nonrigid_bsplines.convert_transform_to_vtk(trans)
            self.transforms.append(vtktrans)
            self.transforms_as_array.append(trans)
        else:
            trans = vtk.vtkTransform()
            self.transforms.append(trans)
            self.transforms_as_array.append(numpy.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0]).astype(float))
        self.subject_ids.append(subject_id)
        
    def remove_mean_from_transforms(self):
        """ Remove mean rotations and mean scaling and mean
         translations from transforms.  The mean rotation and
         translation should not affect the objective function anyway.
         A mean scaling will affect the objective, i.e. if the brain
         shrinks all distances become smaller and the similarity is
         higher. This is not the desired effect."""

        if self.nonrigid:
            # remove any average displacement of each source point.
            # this means the mean of the source points must equal the target point
            transforms_array = numpy.array(self.transforms_as_array)
            meandisp = numpy.mean(transforms_array, 0)
            #if self.verbose:
            print "MEAN DISPLACEMENT:", meandisp

            for transform in self.transforms_as_array:
                transform[:] = transform - meandisp

            transforms_array = numpy.array(self.transforms_as_array)

            meandisp = numpy.mean(transforms_array, 0)
            #if self.verbose:
            print "MEAN DISPLACEMENT 2:", meandisp

        else:
            # Here we are in the affine case
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
        """ Run a single iteration of optimization, multiprocessing over input subjects."""
        
        self.total_iterations += 1
        start_time = time.time()

        # set up progress information saving if this is the first iteration
        if self.total_iterations == 1:
            self.progress_filename = os.path.join(self.output_directory, 'registration_performance.txt')
            progress_file = open(self.progress_filename, 'w')
            print >> progress_file, 'iteration','\t', 'sigma', '\t', 'nonrigid', '\t', 'subject_brain_fibers', '\t', 'fibers_per_subject_in_mean_brain','\t', 'mean_brain_fibers','\t', 'maxfun','\t', 'grid_resolution_if_nonrigid','\t', 'initial_step','\t', 'final_step','\t', 'objective_before','\t', 'objective_after', '\t', 'objective_change', '\t', 'objective_percent_change', '\t', 'mean_function_calls_per_subject','\t', 'min_function_calls_per_subject','\t', 'max_function_calls_per_subject','\t', 'subjects_hitting_maxfun','\t', 'total_subjects','\t', 'subjects_decreased','\t', 'mean_subject_change', '\t', 'mean_subject_decrease_if_decreased', '\t', 'time'
            progress_file.close()
            
        # make a directory for the current iteration
        if self.nonrigid:
            dirname = "iteration_%05d_sigma_%03d_grid_%03d" % (self.total_iterations, self.sigma, self.nonrigid_grid_resolution)
        else:
            dirname = "iteration_%05d_sigma_%03d" % (self.total_iterations, self.sigma)

        outdir = os.path.join(self.output_directory, dirname)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        # make a directory for any intermediate rendering by subprocesses
        outdir_render = os.path.join(outdir, 'progress_rendered')
        if not os.path.exists(outdir_render):
            os.makedirs(outdir_render)

        # Calculate how many fibers are needed to sample from each subject to compute the mean brain at the requested size
        fibers_per_subject = self.mean_brain_size / (len(self.polydatas) - 1)
        if self.verbose:
            print "Fibers per subject for computing mean brain:", fibers_per_subject, "=", self.mean_brain_size, "/",  len(self.polydatas) -1

        # Set up lists of data to pass to the per-subject processes
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
        
        # Each subject will be registered to the current model or "mean brain"
        # Sample fibers from each subject for use in the "mean brain"
        subject_sampled_fibers = list()
        for (input_pd, trans) in zip(self.polydatas, self.transforms):
            pd = wma.filter.downsample(input_pd, fibers_per_subject, verbose=False, random_seed=self.random_seed)
            # apply the current transform to this polydata for computation of mean brain
            transformer = vtk.vtkTransformPolyDataFilter()
            if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                transformer.SetInputData(pd)
            else:
                transformer.SetInput(pd)
            transformer.SetTransform(trans)
            transformer.Update()
            subject_sampled_fibers.append(transformer.GetOutput())
            del transformer
            
        # Loop over all subjects and prepare lists of inputs for subprocesses
        subj_idx = 0
        for input_pd in self.polydatas:
            # Compute current atlas model "mean brain" in a leave-one out fashion.
            # Otherwise, the optimal transform may be identity.
            appender = vtk.vtkAppendPolyData()
            for subj_idx2 in range(len(subject_sampled_fibers)):
                if subj_idx2 != subj_idx:
                    pd = subject_sampled_fibers[subj_idx2]
                    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                        appender.AddInputData(pd)
                    else:
                        appender.AddInput(pd)

            # Convert "mean brain" from vtk to numpy format
            appender.Update()
            mean_brain = appender.GetOutput()
            mean_fibers = wma.fibers.FiberArray()
            mean_fibers.convert_from_polydata(mean_brain, self.points_per_fiber)
            mean_fibers = numpy.array([mean_fibers.fiber_array_r,mean_fibers.fiber_array_a,mean_fibers.fiber_array_s])
            #  R,A,S is the first index
            # then fiber number
            # then points along fiber
            mean_list.append(mean_fibers)

            # Now get the current sample of fibers from the subject for registration to the "mean brain"
            pd = wma.filter.downsample(input_pd, self.subject_brain_size, verbose=False, random_seed=self.random_seed)
            fibers = wma.fibers.FiberArray()
            fibers.convert_from_polydata(pd, self.points_per_fiber)
            fibers_array = numpy.array([fibers.fiber_array_r,fibers.fiber_array_a,fibers.fiber_array_s])
            subject_list.append(fibers_array)

            # Append parameter information to lists of parameters for subprocesses
            sigma_list.append(self.sigma)
            if self.nonrigid:
                mode_list.append('Nonrigid')
            else:
                mode_list.append('Linear')                    
            subj_idx_list.append(subj_idx)
            subj_idx += 1
            iteration_list.append(self.total_iterations)
            outdir_list.append(outdir_render)
            stepsize_list.append(numpy.array([self.initial_step, self.final_step]))
            maxfun_list.append(self.maxfun)
            render_list.append(self.render)
            grid_resolution_list.append(self.nonrigid_grid_resolution)
            
        # Multiprocess over subjects
        print "\nITERATION", self.total_iterations, "STARTING MULTIPROCESSING. NUMBER OF JOBS:", self.parallel_jobs, "\n"

        # note we can't pass vtk objects to subprocesses since they can't be pickled.
        ret = Parallel(
            n_jobs=self.parallel_jobs, verbose=self.parallel_verbose)(
                delayed(congeal_multisubject_inner_loop)(fixed, moving, initial_transform, mode, sigma, subject_idx, iteration_count, output_directory, step_size, maxfun, render, grid_resolution)
                for (fixed, moving, initial_transform, mode, sigma, subject_idx, iteration_count, output_directory, step_size, maxfun, render, grid_resolution) in zip(mean_list, subject_list, self.transforms_as_array, mode_list, sigma_list, subj_idx_list, iteration_list, outdir_list, stepsize_list, maxfun_list, render_list, grid_resolution_list))

            
        #print "RETURNED VALUES", ret
        
        # Progress reporting: loop over all registration outputs.
        # Get the current transform for each subject and report the
        # objective values to the user by printing, saving, and plotting.
        self.transforms_as_array = list()
        objective_total_before = 0.0
        objective_total_after = 0.0
        sidx = 0
        functions_per_subject = list()
        objective_changes_per_subject = list()
        decreases = list()
        if HAVE_PLT:
            plt.figure(0)
            plt.title('Iteration '+str(self.total_iterations)+' Objective Values for All Subjects')
            plt.xlabel('objective function computations')
            plt.ylabel('objective value')
            
        for (trans, objectives, diff) in ret:
            self.transforms_as_array.append(trans)
            print "Iteration:", self.total_iterations, "Subject:", sidx, "Objective function computations:", len(objectives), "change", diff
            functions_per_subject.append(len(objectives))
            # Normalize by the number of fibers so this is comparable across iterations if sigma does not change
            objectives = numpy.divide(objectives, self.subject_brain_size)
            # Compute total objective for progress reporting.
            objective_total_before += objectives[0]
            if diff < 0:
                objective_total_after += objectives[-1]
                decreases.append(diff)
            else:
                objective_total_after += objectives[0]
            objective_changes_per_subject.append(diff)
            sidx += 1
            if HAVE_PLT:
                plt.figure(0)
                plt.plot(objectives, 'o-', label=sidx)

        number_of_subjects = sidx
        functions_per_subject = numpy.array(functions_per_subject)
        objective_changes_per_subject = numpy.array(objective_changes_per_subject)
        decreases = numpy.array(decreases)

        self.objectives_before.append(objective_total_before)
        self.objectives_after.append(objective_total_after)
        total_change =  self.objectives_after[-1] - self.objectives_before[-1]
        percent_change = total_change / self.objectives_before[-1]
        print "Iteration:", self.total_iterations, "TOTAL objective change:",  total_change
        print "Iteration:", self.total_iterations, "PERCENT objective change:",  percent_change

        if HAVE_PLT:
            plt.figure(0)
            if self.nonrigid:
                fname_fig_base = "iteration_%05d_sigma_%03d_grid_%03d" % (self.total_iterations, self.sigma, self.nonrigid_grid_resolution)
            else:
                fname_fig_base = "iteration_%05d_sigma_%03d_" % (self.total_iterations, self.sigma)
            # Place the legend below the plot so it does not overlap it when there are many subjects
            #lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=False, shadow=False, ncol=1)
            fname_fig = 'objectives_per_subject_' + fname_fig_base + '.pdf'
            # save everything even if the legend is long and goes off the plot
            #plt.savefig(os.path.join(outdir, fname_fig), bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.savefig(os.path.join(outdir, fname_fig))
            plt.close(0)

        progress_file = open(self.progress_filename, 'a')
        elapsed_time = time.time() - start_time
        if len(decreases) == 0:
            mean_decreases = 0.0
        else:
            mean_decreases = numpy.mean(decreases)
        print >> progress_file, self.total_iterations,'\t', self.sigma, '\t', self.nonrigid, '\t', self.subject_brain_size, '\t', fibers_per_subject,'\t', self.mean_brain_size,'\t', self.maxfun,'\t', self.nonrigid_grid_resolution,'\t', self.initial_step,'\t', self.final_step,'\t', self.objectives_before[-1],'\t', self.objectives_after[-1],'\t', total_change,'\t',  percent_change,'\t', numpy.mean(functions_per_subject), '\t', numpy.min(functions_per_subject), '\t', numpy.max(functions_per_subject), '\t', numpy.sum(functions_per_subject >= self.maxfun), '\t', number_of_subjects,'\t', len(decreases),'\t', numpy.mean(objective_changes_per_subject), '\t', mean_decreases, '\t', elapsed_time
        progress_file.close()

        # remove_mean_from_transforms
        self.remove_mean_from_transforms()

        # update our transforms list for the next iteration
        self.transforms = list()
        for trans in self.transforms_as_array:
            if self.nonrigid:
                vtktrans = wma.register_two_subjects_nonrigid_bsplines.convert_transform_to_vtk(trans)
                #print vtktrans
            else:
                vtktrans = wma.register_two_subjects.convert_transform_to_vtk(trans)
                #print vtktrans.GetMatrix()
            self.transforms.append(vtktrans)

        # save the current transforms to disk
        # save the transforms to text files
        start_time = time.time()
        wma.io.write_transforms_to_itk_format(self.transforms, outdir, self.subject_ids)
        elapsed_time = time.time() - start_time
        print "WRITE TXFORMS:", elapsed_time
        
    def save_transformed_polydatas(self, intermediate_save=False, midsag_symmetric=False):
        """ Output polydatas for final or intermediate iterations. """
        
        # this can be slow so notify the user what is happening
        print "\nSAVING TRANSFORMED TRACTOGRAPHY FROM ITERATION", self.total_iterations, "\n"
        
        transform_list = self.transforms
        subject_id_list = self.subject_ids
        if midsag_symmetric:
            transform_list = transform_list[::2]
            subject_id_list = subject_id_list[::2]

        if intermediate_save:
            # Make a directory for the current iteration.
            # This directory name must match the one created above in the iteration.
            if self.nonrigid:
                dirname = "iteration_%05d_sigma_%03d_grid_%03d" % (self.total_iterations, self.sigma, self.nonrigid_grid_resolution)
            else:
                dirname = "iteration_%05d_sigma_%03d" % (self.total_iterations, self.sigma)
            outdir = os.path.join(self.output_directory, dirname)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            
            # make a directory for output polydatas
            outdir_pds = os.path.join(outdir, 'transformed_output')
            if not os.path.exists(outdir_pds):
                os.makedirs(outdir_pds)

            wma.io.transform_polydatas_from_disk(self.input_directory, transform_list, outdir_pds, parallel_jobs=self.parallel_jobs)

        else:
            # make a directory for the final output
            outdir = os.path.join(self.output_directory, 'output_tractography')
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            if self.verbose:
                for trans in  self.transforms:
                    if self.nonrigid:
                        print trans
                    else:
                        print trans.GetMatrix()

            wma.io.transform_polydatas_from_disk(self.input_directory, transform_list, outdir, parallel_jobs=self.parallel_jobs)

            # Save the current atlas representation to disk.
            # Right now this is all the input fibers from all subjects.
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
            start_time = time.time()
            wma.io.write_transforms_to_itk_format(transform_list, outdir, subject_id_list)
            elapsed_time = time.time() - start_time
            print "WRITE TXFORMS:", elapsed_time


def congeal_multisubject_inner_loop(mean, subject, initial_transform, mode, sigma, subject_idx, iteration_count, output_directory, step_size, maxfun, render, grid_resolution):

    """This is the code executed by each subprocess that launches the

    registration of one subject to the current atlas model or mean brain.
    """
    
    #print "\n BEGIN ITERATION", iteration_count, "subject", subject_idx, "sigma:", sigma, "mean brain:", mean.shape, "subject:", subject.shape, "initial transform length:", len(initial_transform), "steps:", step_size[0], step_size[1], "maxfun:", maxfun, type(initial_transform), "Grid:", grid_resolution, "Mode:", mode, "initial transform:", initial_transform,
    
    # Set up registration objects and parameters that are specific to affine vs nonrigid
    if mode == 'Linear':
        register = wma.register_two_subjects.RegisterTractography()
        register.process_id_string = "_subject_%05d_iteration_%05d_sigma_%03d" % (subject_idx, iteration_count, sigma) 

    elif mode == "Nonrigid":
        register = wma.register_two_subjects_nonrigid_bsplines.RegisterTractographyNonrigid()
        register.nonrigid_grid_resolution = grid_resolution
        register.initialize_nonrigid_grid()
        register.process_id_string = "_subject_%05d_iteration_%05d_sigma_%03d_grid_%03d" % (subject_idx, iteration_count, sigma, grid_resolution) 

    else:
        print "ERROR: Unknown registration mode"

    # Set up parameters that are used for both affine and nonrigid
    register.maxfun = maxfun
    register.sigma = sigma        
    register.parallel_jobs = 1
    register.threshold = 0
    register.fixed = mean
    register.moving = subject
    register.initial_transform = initial_transform
    register.verbose = False
    register.output_directory = output_directory
    register.initial_step = step_size[0]
    register.final_step = step_size[1]
    register.render = render

    # Run the current iteration of optimization.    
    register.compute()

    # Only update the transform if the objective function improved.
    # With affine registration, some subjects may have converged already to the current model.
    if mode == 'Linear':
        obj_diff = register.objective_function_values[-1] - register.objective_function_values[0]
    elif mode == "Nonrigid":
        obj_diff = numpy.min(register.objective_function_values) - register.objective_function_values[0]
    #print "\n END ITERATION", iteration_count, "subject", subject_idx, "OBJECTIVE CHANGE:", obj_diff, register.objective_function_values[-1] - register.objective_function_values[0], register.final_transform
    if obj_diff < 0:
        #print "UPDATING MATRIX"
        return register.final_transform, register.objective_function_values, obj_diff
    else:
        return initial_transform, register.objective_function_values, obj_diff


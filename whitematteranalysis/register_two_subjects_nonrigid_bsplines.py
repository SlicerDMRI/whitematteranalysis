# -*- coding: utf-8 -*-

""" register.py

implementation of fiber tractography registration (group)

class RegisterTractographyNonrigid


"""

import os
import sys
import time

import numpy as np
import scipy.optimize
import vtk
import vtk.util.numpy_support

import whitematteranalysis as wma

# debug only
#import resource


class RegisterTractographyNonrigid(wma.register_two_subjects.RegisterTractography):

    def constraint(self, x_current):
        # Make sure the optimizer is searching in a reasonable region.
        # TEST: Don't let the translations grow too large
        #penalty = 10.0 - np.mean(np.abs(x_current))
        penalty = 30.0 - np.mean(np.abs(x_current * 0.01))

        # progress report sometimes
        
        iters = self.objective_computations
        #len(self.objective_function_values)
        print_every = int(self.maxfun  / 10)
        if iters % print_every == 0:
            elapsed_time = time.time() - self.last_time
            self.last_time = time.time()
            self.total_time = time.time() - self.start_time
            #print iters, "/", self.maxfun, "Total time:", self.total_time, "Last iters:", elapsed_time, "Per iter:", elapsed_time / print_every, "Memory:", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            progress_file = open(self.progress_filename, 'a')
            print(f"{iters} / {self.maxfun} Total time: {self.total_time} Last iters: {elapsed_time} Per iter: {elapsed_time / print_every} Memory: {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}", file=progress_file)
            progress_file.close()    
        return penalty
    
    def __init__(self):
        # parameters that should be set by user
        self.sigma = 5
        self.process_id_string = ""
        self.output_directory = None
        
        # performance options that should be set by user
        self.verbose = False
        self.render = False
        
        # optimizer parameters that should be set by user
        self.maxfun = 300

        # output of registration
        self.objective_function_values = list()
        #self.objective_function_values = [1.0, 0.0]
        self.objective_computations = 0
        self.final_transform = None
        self.progress_filename = None
        self.start_time = None
        self.total_time = 0.0
        self.last_time = None

        # subject data that must be input
        #self.fixed = None
        #self.moving = None
        #self.initial_step = 5
        #self.final_step = 2

        # set up default grid
        #self.nonrigid_grid_resolution = 3
        #self.nonrigid_grid_resolution = 5
        self.nonrigid_grid_resolution = 6
        self.initialize_nonrigid_grid()

        # transform we optimize over
        self.initial_transform = self.displacement_field_numpy

        # internal recordkeeping
        self.iterations = 0

        # works for Cobyla. May not be needed.
        #self.scaling = 100

        # works well for BFGS
        self.scaling = 0.1
        # not ok for bfgs
        #self.scaling = 0.05
        
        # keep track of the best objective we have seen so far to return that when computation stops.
        self.minimum_objective = np.inf

        # choice of optimization method
        #self.optimizer = "Powell"
        #self.optimizer = "Cobyla"
        self.optimizer = "BFGS"

    def initialize_nonrigid_grid(self):
        res = self.nonrigid_grid_resolution
        self.displacement_field_numpy = np.zeros(res*res*res*3)
        #self.displacement_field_vtk = numpy_support.numpy_to_vtk(num_array=NumPy_data.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
        
    def objective_function(self, current_x):
        """ The actual objective used in registration.  Function of
        the current x in search space, as well as parameters of the
        class: threshold, sigma. Compares sampled fibers from moving
        input, to all fibers of fixed input."""
        self.objective_computations += 1
    
        #scaled_current_x = current_x * 0.01
        scaled_current_x = current_x / self.scaling
        
        # get and apply transforms from current_x
        moving = self.transform_fiber_array_numpy(self.moving, scaled_current_x)

        # compute objective
        obj = wma.register_two_subjects.inner_loop_objective(self.fixed, moving, self.sigma * self.sigma)

        # keep track of minimum objective so far and its matching transform
        if obj < self.minimum_objective:
            #print "OBJECTIVE:", obj, "PREV MIN",  self.minimum_objective
            self.minimum_objective = obj
            # must copy current_x into allocated memory space to keep the value
            self.final_transform[:] = scaled_current_x

        # save objective function value for analysis of performance
        self.objective_function_values.append(obj)

        # temporary for testing
        if self.optimizer == "BFGS":
            self.constraint(current_x)

        if self.verbose:
            print(f"O: {obj} X: {scaled_current_x}")
        #print "X:", self._x_opt

        # Stop the optimizer if needed
        if self.objective_computations > self.maxfun:
            raise RuntimeError("The optimizer needs to stop.")

        return obj

    def transform_fiber_array_numpy(self, in_array, transform):
        """Transform in_array of R,A,S by transform (a list of source points).  Transformed fibers are returned.
        """
        (dims, number_of_fibers, points_per_fiber) = in_array.shape
        out_array = np.zeros(in_array.shape)

        vtktrans = convert_transform_to_vtk(transform)
        #print "2:", vtktrans
        #vtktrans = vtk.vtkTransform()
        
        # Transform moving fiber array by applying transform to original fibers
        for lidx in range(0, number_of_fibers):
            for pidx in range(0, points_per_fiber):
                pt = vtktrans.TransformPoint(in_array[0, lidx, pidx],
                                            in_array[1, lidx, pidx], 
                                            in_array[2, lidx, pidx])
                out_array[0, lidx, pidx] = pt[0]
                out_array[1, lidx, pidx] = pt[1]
                out_array[2, lidx, pidx] = pt[2]

        #print in_array[0, lidx, pidx], in_array[1, lidx, pidx], in_array[2, lidx, pidx], "===>>>", out_array[0, lidx, pidx], out_array[1, lidx, pidx], out_array[2, lidx, pidx]
        del vtktrans

        ## uncomment for testing only
        ## # convert it back to a fiber object and render it
        ## global __render_count
        ## if (np.mod(__render_count, 500) == 0) & False:
        ##     fiber_array = wma.fibers.FiberArray()
        ##     fiber_array.fiber_array_r = out_array[0,:,:]
        ##     fiber_array.fiber_array_a = out_array[1,:,:]
        ##     fiber_array.fiber_array_s = out_array[2,:,:]
        ##     fiber_array.points_per_fiber = points_per_fiber
        ##     fiber_array.number_of_fibers = number_of_fibers
        ##     pd = fiber_array.convert_to_polydata()
        ##     ren = wma.render.render(pd, number_of_fibers, verbose=False)
        ##     ren.save_views('.', 'moving_{0:05d}_'.format(__render_count)+str(time.clock())[-5:-1])
        ##     del ren
        ## __render_count += 1

        return out_array

    def compute(self):

        """ Run the registration.  Add subjects first (before calling
        compute). Then call compute several times, using different
        parameters for the class, for example first just for
        translation."""
        print(f"OPTIMIZER: {self.optimizer}")
        self.start_time = time.time()
        self.total_time = 0.0
        self.progress_filename = os.path.join(self.output_directory, f"log{self.process_id_string}log")
        print(self.progress_filename)
        progress_file = open(self.progress_filename, 'w')
        print(f"Starting computation, time: {self.start_time}", file=progress_file)
        progress_file.close()
        # subject data must be input first. No check here for speed
        #self.fixed = None
        #self.moving = None
        #self.initial_transform = None

        # This is left if needed in future for debugging.
        # convert it back to a fiber object and render it
        ## (dims, number_of_fibers_moving, points_per_fiber) = self.moving.shape
        ## fiber_array = wma.fibers.FiberArray()
        ## fiber_array.fiber_array_r = self.moving[0,:,:]
        ## fiber_array.fiber_array_a = self.moving[1,:,:]
        ## fiber_array.fiber_array_s = self.moving[2,:,:]
        ## fiber_array.points_per_fiber = points_per_fiber
        ## fiber_array.number_of_fibers = number_of_fibers_moving
        ## pd = fiber_array.convert_to_polydata()
        ## ren = wma.render.render(pd, number_of_fibers_moving, verbose=False)
        ## ren.save_views('.', 'moving_brain_{0:05d}_'.format(self.iterations)+str(time.clock())[-5:-1])
        ## #ren.save_views('.', 'moving_brain_{0:05d}'.format(self.iterations))
        ## del ren

        # For debugging/monitoring of progress
        if self.render:
            (dims, number_of_fibers_fixed, points_per_fiber) = self.fixed.shape
            fiber_array = wma.fibers.FiberArray()
            fiber_array.fiber_array_r = self.fixed[0,:,:]
            fiber_array.fiber_array_a = self.fixed[1,:,:]
            fiber_array.fiber_array_s = self.fixed[2,:,:]
            fiber_array.points_per_fiber = points_per_fiber
            fiber_array.number_of_fibers = number_of_fibers_fixed
            pd2 = fiber_array.convert_to_polydata()
            ren = wma.render.render(pd2, number_of_fibers_fixed, verbose=False)
            # save low-res images for speed
            ren.magnification = 3
            ren.save_views(self.output_directory, f'fixed_brain_{self.process_id_string}')
            del ren
                
        self.iterations += 1
        self.final_transform = np.zeros(self.initial_transform.shape)

        if self.verbose:
            print(f"<{os.path.basename(__file__)}> Initial value for X: {self.initial_transform}")

        progress_file = open(self.progress_filename, 'a')
        self.total_time = time.time() - self.start_time
        print(f"INITIAL transform shape {self.initial_transform.shape} Init time: {self.total_time}", file=progress_file)
        progress_file.close()

        # initialize time for objective function computations
        self.last_time = time.time()

        if self.optimizer == "Cobyla":

            print(f"INITIAL transform shape {self.initial_transform.shape}")
            # Optimize using cobyla. Allows definition of initial and
            # final step size scales (rhos), as well as constraints.  Here
            # we use the constraints to encourage that the transform stays a transform.
            # note disp 0 turns off all display
            not_used = scipy.optimize.fmin_cobyla(self.objective_function,
                                                  self.initial_transform * self.scaling, self.constraint,
                                                  maxfun=self.maxfun, rhobeg=self.initial_step * self.scaling,
                                                  rhoend=self.final_step * self.scaling, disp=1)
        elif self.optimizer == "BFGS":
            # Test optimization with BFGS
            # (Broyden-Fletcher-Goldfarb-Shanno algorithm) refines at each
            # step an approximation of the Hessian.  L-BFGS:
            # Limited-memory BFGS Sits between BFGS and conjugate
            # gradient: in very high dimensions (> 250) the Hessian matrix
            # is too costly to compute and invert. L-BFGS keeps a low-rank
            # version. In addition, the scipy version,
            # scipy.optimize.fmin_l_bfgs_b(), includes box bounds.
            # Note If you do not specify the gradient to the L-BFGS
            # solver, you need to add approx_grad=1
            # list of (min,max) pairs for the values being optimized. Assume we never should move by >30mm
            #bounds = list()
            #for lm in self.target_landmarks:
            #    bounds.append((lm-30,lm+30))
            ## (self.final_transform, f, dict) = scipy.optimize.fmin_l_bfgs_b(self.objective_function,
            ##                                                                self.initial_transform,
            ##                                                                approx_grad = True,
            ##                                                                maxfun=self.maxfun,
            ##                                                                maxiter=self.maxfun,
            ##                                                                factr=1e12,
            ##                                                                epsilon=self.final_step,
            ##                                                                iprint=0,
            ##                                                                bounds=bounds)
            #maxiter=1,
            try:
                (not_used, f, dict) = scipy.optimize.fmin_l_bfgs_b(self.objective_function,
                                                                           self.initial_transform * self.scaling,
                                                                           approx_grad = True,
                                                                           maxfun=self.maxfun,
                                                                           maxiter=2,
                                                                           factr=1e12,
                                                                           epsilon=self.final_step * self.scaling,
                                                                           iprint=1)
            except:
                print(f"EXCEPTION WAS CAUGHT, total objectives: {self.objective_computations}")
            #print f, dict

        elif self.optimizer == "Powell":
            # Test optimization with Powell's method
            # Powell's method is a conjugate direction method.
            #(self.final_transform, fopt, direc, iters, funcalls, warnflag, allvecs)
            (self.final_transform, fopt, direc, iters, funcalls, warnflag) = scipy.optimize.fmin_powell(self.objective_function,
                                                                            self.initial_transform,
                                                                            xtol=self.initial_step,
                                                                            ftol=self.final_step,
                                                                            maxfun=self.maxfun,
                                                                            maxiter=self.maxfun,
                                                                            disp=1, full_output=True)

            print(f"TRANS: {self.final_transform} FLAG: {warnflag}")

        else:
            raise NotImplementedError(
                f"Workflow not implemented for optimizer: {self.optimizer}.")

        progress_file = open(self.progress_filename, 'a')
        self.total_time = time.time() - self.start_time
        print(f"Done optimizing. TOTAL TIME: {self.total_time}", file=progress_file)
        progress_file.close()

        if self.verbose:
            print(f"O: {self.objective_function_values}")

        # Return output transforms from this iteration
        return self.final_transform

def convert_numpy_array_to_vtk_points(inarray):
    """ Convert numpy array or flat list of points to vtkPoints."""
    
    number_of_points = len(inarray)/3
    vtk_points = vtk.vtkPoints()
    vtk_points.SetNumberOfPoints(number_of_points)
    idx = 0
    for pt in zip(inarray[::3], inarray[1::3], inarray[2::3]):
        #print pt
        vtk_points.SetPoint(idx, pt[0], pt[1], pt[2])
        idx += 1
    return vtk_points

def convert_transform_to_vtk(transform):
    """Produce an output vtkBSplineTransform corresponding to the

    registration results. Input is a numpy array corresponding to the displacement field.
    """
    displacement_field_vtk = vtk.util.numpy_support.numpy_to_vtk(num_array=transform, deep=True, array_type=vtk.VTK_FLOAT)
    displacement_field_vtk.SetNumberOfComponents(3)
    displacement_field_vtk.SetName('DisplacementField')
    grid_image = vtk.vtkImageData()
    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        grid_image.AllocateScalars(vtk.VTK_FLOAT, 3)
        grid_image.GetPointData().SetScalars(displacement_field_vtk)
    else:
        grid_image.SetScalarTypeToFloat()
        grid_image.SetNumberOfScalarComponents(3)
        grid_image.GetPointData().SetScalars(displacement_field_vtk)
        grid_image.Update()
    #print "CONVERT TXFORM 1:", grid_image.GetExtent(), displacement_field_vtk.GetSize()

    # this is a hard-coded assumption about where the polydata is located in space.
    # other code should check that it is centered.
    # This code uses a grid of 240mm x 240mm x 240mm
    #spacing origin extent
    num_vectors = len(transform) / 3
    dims = round(np.power(num_vectors, 1.0/3.0))
    # This MUST correspond to the size used in congeal_multisubject update_nonrigid_grid
    #size_mm = 240.0
    size_mm = 200.0
    origin = -size_mm / 2.0
    # assume 240mm x 240mm x 240mm grid
    spacing = size_mm / (dims - 1)
    grid_image.SetOrigin(origin, origin, origin)
    grid_image.SetSpacing(spacing, spacing, spacing)
    #grid_image.SetExtent(0, dims-1.0, 0, dims-1.0, 0, dims-1.0)
    grid_image.SetDimensions(int(dims), int(dims), int(dims))
    #print "CONVERT TXFORM:", num_vectors, dims, int(dims), dims-1.0, grid_image.GetExtent(), 
    
    #print "GRID:", grid_image
    coeff = vtk.vtkImageBSplineCoefficients()
    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        coeff.SetInputData(grid_image)
    else:
        coeff.SetInput(grid_image)

    coeff.Update()
    # this was in the test code.
    coeff.UpdateWholeExtent()
    #print "TX:", transform.shape, transform, displacement_field_vtk, grid_image.GetExtent(), coeff.GetOutput().GetExtent()

    vtktrans = vtk.vtkBSplineTransform()
    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        vtktrans.SetCoefficientData(coeff.GetOutput())
    else:
        vtktrans.SetCoefficients(coeff.GetOutput())
    vtktrans.SetBorderModeToZero()

    ## print "~~~~~~~~~~~~~~~~~~~~~~~~"
    ## print "COEFF:",  coeff.GetOutput()
    ## print "*********"
    ## print "COEFF2:", vtktrans.GetCoefficients()
    ## print "======="
    
    return vtktrans


 


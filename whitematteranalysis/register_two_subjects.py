# -*- coding: utf-8 -*-

""" register.py

implementation of fiber tractography registration (group)

class RegisterTractography


"""

import os
import sys
import time

import numpy as np
import scipy.optimize
import vtk

import whitematteranalysis as wma


class RegisterTractography:

    def constraint(self, x_current):
        penalty = 0
        
        # The values should be a registration matrix.
        # Make sure the optimizer is searching in a reasonable region.
        transform = np.divide(x_current,self.transform_scaling)
        
        # want translation reasonable, say -100 to 100 range
        max_trans = 100
        if np.abs(transform[0]) > max_trans:
            penalty -= 1
        if np.abs(transform[1]) > max_trans:
            penalty -= 1
        if np.abs(transform[2]) > max_trans:
            penalty -= 1
        # want rotation reasonable, say +/- 90 degrees is the max we could reasonably handle
        max_rot = 90
        if np.abs(transform[3]) > max_rot:
            penalty -= 1
        if np.abs(transform[4]) > max_rot:
            penalty -= 1
        if np.abs(transform[5]) > max_rot:
            penalty -= 1
        # want scaling between 0.9 and 1.2. In practice, only shrinking is an issue.
        # Don't let it shrink by more than 0.75
        #min_scale = 0.75
        min_scale = 0.95
        if transform[6] < min_scale:
            penalty -= 1
        if transform[7] < min_scale:
            penalty -= 1
        if transform[8] < min_scale:
            penalty -= 1
        # want shear angles reasonable like rotation
        max_shear_angle = 90
        if np.abs(transform[9]) > max_shear_angle:
            penalty -= 1
        if np.abs(transform[10]) > max_shear_angle:
            penalty -= 1
        if np.abs(transform[11]) > max_shear_angle:
            penalty -= 1
        if np.abs(transform[12]) > max_shear_angle:
            penalty -= 1
        if np.abs(transform[13]) > max_shear_angle:
            penalty -= 1
        if np.abs(transform[14]) > max_shear_angle:
            penalty -= 1
            
        return penalty
    
    def __init__(self):
        # parameters that can be set by user
        self.sigma = 10
        self.process_id_string = ""
        self.output_directory = None
        
        # performance options set by user
        self.verbose = 0
        self.render = False
        
        # optimizer parameters set by user
        self.maxfun = 300

        # output of registration
        self.objective_function_values = list()
        self.final_transform = None
        
        # subject data to be input
        self.fixed = None
        self.moving = None
        self.initial_step = 5
        self.final_step = 2

        # Order the components as translate then rotate so optimizer will translate first.
        # trans, rot, scale, shear:
        self.initial_transform = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0])

        # the scaling components should have much smaller changes than the others, so scale them accordingly in optimizer search space
        #self.transform_scaling = np.array([1, 1, 1, .5, .5, .5, 300, 300, 300,  1, 1, 1, 1, 1, 1])
        self.transform_scaling = np.array([1, 1, 1, .5, .5, .5, 200, 200, 200,  1, 1, 1, 1, 1, 1])

        # the registration mode includes trans, rot, scale, shear.
        # so for rigid, use mode =  [1, 1, 0, 0]
        # so for translation only, use mode =  [1, 1, 0, 0]
        self.mode=[1,1,1,1]

        # internal recordkeeping
        self._x_opt = None
        self.iterations = 0
        ## self.moving_points = None
        ## self.number_of_fibers_moving = None
        ## self.points_per_fiber = None

        # choice of optimization method
        #self.optimizer = "Powell"
        self.optimizer = "Cobyla"
        
    def objective_function(self, current_x):
        """ The actual objective used in registration.  Function of
        the current x in search space, as well as parameters of the
        class: sigma. Compares sampled fibers from moving
        input, to all fibers of fixed input."""

        # keep track of current x value
        self._x_opt = current_x

        # get and apply transforms from current_x
        ## t1 = time.time()
        moving = transform_fiber_array_numpy(self.moving, current_x, self.mode)
        ## t2 = time.time()
        ##moving = transform_fiber_array_numpy(self.moving_points, self.number_of_fibers_moving, self.points_per_fiber, current_x)
        ## t3 = time.time()
        ## diff = np.abs(movingOLD - moving)
        ## print "DIFFERENCE IN TXFORMS:", np.max(diff), "TIME:", t2-t1, t3-t2

        # compute objective
        obj = inner_loop_objective(self.fixed, moving, self.sigma * self.sigma)

        # save objective function value for analysis of performance
        self.objective_function_values.append(obj)

        if self.verbose:
            print(f"O: {obj} X: {self._x_opt}")

        return obj

    def compute(self):

        """ Run the registration.  Add subjects first (before calling
        compute). Then call compute several times, using different
        sigma, to perform multiscale registration."""

        print("OPTIMIZER:", self.optimizer)

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

        # In the end this was much slower so stick with all data in numpy arrays.
        # Create a vtkPoints object with the input points
        ## self.moving_points = vtk.vtkPoints()
        ## (dims, self.number_of_fibers_moving, self.points_per_fiber) = self.moving.shape
        ## count = 0
        ## for lidx in range(0, self.number_of_fibers_moving):
        ##     for pidx in range(0, self.points_per_fiber):
        ##         self.moving_points.InsertNextPoint(self.moving[0, lidx, pidx], self.moving[1, lidx, pidx], self.moving[2, lidx, pidx])
        ##         count += 1

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

        if self.verbose:
            print(f"<{os.path.basename(__file__)}> Initial value for X: {self.initial_transform}")


        if self.optimizer == "Cobyla":

            # Optimize using cobyla. Allows definition of initial and
            # final step size scales (rhos), as well as constraints.  Here
            # we use the constraints to encourage that the transform stays a transform.
            # note disp 0 turns off all display
            self.final_transform = scipy.optimize.fmin_cobyla(self.objective_function,
                                                      np.multiply(self.initial_transform,self.transform_scaling), self.constraint,
                                                      maxfun=self.maxfun, rhobeg=self.initial_step,
                                                      rhoend=self.final_step, disp=0)
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
            (self.final_transform, f, dict) = scipy.optimize.fmin_l_bfgs_b(self.objective_function,
                                                                           np.multiply(self.initial_transform,self.transform_scaling),
                                                                           approx_grad = True,
                                                                           maxfun=self.maxfun,
                                                                           maxiter=self.maxfun,
                                                                           factr=1e12,
                                                                           epsilon=self.final_step,
                                                                           iprint=0)
            print(f, dict)

        elif self.optimizer == "Powell":
            # Test optimization with Powell's method
            # Powell's method is a conjugate direction method.
            #(self.final_transform, fopt, direc, iters, funcalls, warnflag, allvecs)
            (self.final_transform, fopt, direc, iters, funcalls, warnflag) = scipy.optimize.fmin_powell(self.objective_function,
                                                                            np.multiply(self.initial_transform,self.transform_scaling),
                                                                            xtol=self.initial_step,
                                                                            ftol=self.final_step,
                                                                            maxfun=self.maxfun,
                                                                            maxiter=self.maxfun,
                                                                            disp=1, full_output=True)

            print(f"FLAG: {warnflag}")

        else:
            raise NotImplementedError(
                f"Workflow not implemented for optimizer: {self.optimizer}.")

        self.final_transform = np.divide(self.final_transform, self.transform_scaling)

        # modify the output according to the mode. Note: ideally for
        # rigid or translation only, should make the optimizer search
        # space smaller instead of ignoring some dimensions. This is a
        # test for now.
        if not self.mode[0] == 1:
            self.final_transform[0] = 0.0
            self.final_transform[1] = 0.0
            self.final_transform[2] = 0.0
        if not self.mode[1] == 1:
            self.final_transform[3] = 0.0
            self.final_transform[4] = 0.0
            self.final_transform[5] = 0.0
        if not self.mode[2] == 1:
            self.final_transform[6] = 1.0
            self.final_transform[7] = 1.0
            self.final_transform[8] = 1.0
        if not self.mode[3] == 1:
            self.final_transform[9] = 0.0
            self.final_transform[10] = 0.0
            self.final_transform[11] = 0.0
            self.final_transform[12] = 0.0
            self.final_transform[13] = 0.0
            self.final_transform[14] = 0.0

        tx = self.final_transform
        print(f"TRANS: {tx[0]} {tx[1]} {tx[2]} ROT: {tx[3]} {tx[4]} {tx[5]} SCALE: {tx[6]} {tx[7]} {tx[8]} SHEAR: {tx[9]} {tx[10]} {tx[11]} {tx[12]} {tx[13]} tx[14] MODE: {self.mode} MODE0: {self.mode[0]}")
                                
        # Return output transform from this iteration
        return self.final_transform


def inner_loop_objective(fixed, moving, sigmasq):
    """The code called within the objective_function to find the negative log

    probability of one brain given all other brains.
    """

    (dims, number_of_fibers_moving, points_per_fiber) = moving.shape
    # number of compared fibers (normalization factor)
    (dims, number_of_fibers_fixed, points_per_fiber) = fixed.shape

    probability = np.zeros(number_of_fibers_moving) + 1e-20
    
    # Loop over fibers in moving. Find total probability of
    # each fiber using all fibers from fixed.
    for idx in range(number_of_fibers_moving):
        probability[idx] += total_probability_numpy(moving[:,idx,:], fixed,
                sigmasq)

    # Divide total probability by number of fibers in the atlas ("mean
    # brain").  This neglects Z, the normalization constant for the
    # pdf, which would not affect the optimization.
    probability /= number_of_fibers_fixed
    #print np.min(probability), np.max(probability)
    # add negative log probabilities of all fibers in this brain.
    entropy = np.sum(- np.log(probability))
    return entropy

def total_probability_numpy(moving_fiber, fixed_fibers, sigmasq):
    """Compute total probability for moving fiber when compared to all fixed

    fibers.
    """
    distance = fiber_distance_numpy(moving_fiber, fixed_fibers)
    probability = np.exp(np.divide(-distance, sigmasq))
    return np.sum(probability)

def fiber_distance_numpy(moving_fiber, fixed_fibers):
    """
    Find pairwise fiber distance from fixed fiber to all moving fibers.
    """
    # compute pairwise fiber distances along fibers
    distance_1 = _fiber_distance_internal_use_numpy(moving_fiber, fixed_fibers)
    distance_2 = _fiber_distance_internal_use_numpy(moving_fiber, fixed_fibers, reverse_fiber_order=True)
    
    # choose the lowest distance, corresponding to the optimal fiber
    # representation (either forward or reverse order)
    return np.minimum(distance_1, distance_2)

def _fiber_distance_internal_use_numpy(moving_fiber, fixed_fibers, reverse_fiber_order=False):
    """Compute the total fiber distance from one fiber to an array of many

    fibers.  This function does not handle equivalent fiber
    representations. For that use fiber_distance, above.
    """
    
    # compute the distance from this fiber to the array of other fibers
    if reverse_fiber_order:
        ddx = fixed_fibers[0,:,:] - moving_fiber[0,::-1]
        ddy = fixed_fibers[1,:,:] - moving_fiber[1,::-1]
        ddz = fixed_fibers[2,:,:] - moving_fiber[2,::-1]
    else:
        ddx = fixed_fibers[0,:,:] - moving_fiber[0,:]
        ddy = fixed_fibers[1,:,:] - moving_fiber[1,:]
        ddz = fixed_fibers[2,:,:] - moving_fiber[2,:]

    #print "MAX abs ddx:", np.max(np.abs(ddx)), "MAX ddy:", np.max(np.abs(ddy)), "MAX ddz:", np.max(np.abs(ddz))
    #print "MIN abs ddx:", np.min(np.abs(ddx)), "MIN ddy:", np.min(np.abs(ddy)), "MIN ddz:", np.min(np.abs(ddz))
    
    distance = np.square(ddx)
    distance += np.square(ddy)
    distance += np.square(ddz)

    # Use the mean distance as it works better than Hausdorff-like distance
    return np.mean(distance, 1)

    # This is how to test Hausdorff-like distance. Left here for documentation.
    # Hausdorff
    # take max along fiber
    #return np.max(distance, 1)

    
def transform_fiber_array_numpy(in_array, transform, mode=[1,1,1,1]):
    """Transform in_array of R,A,S by transform (15 components, rotation about
    R,A,S, translation in R, A, S,  scale along R, A, S, and
    shear. Fibers are assumed to be in RAS (or LPS as long as all inputs
    are consistent).  Transformed fibers are returned.
    """
    (dims, number_of_fibers, points_per_fiber) = in_array.shape
    out_array = np.zeros(in_array.shape)

    vtktrans = convert_transform_to_vtk(transform, scaled=True, mode=mode)

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

def transform_fiber_array_numpyNOTUSED(moving_points, number_of_fibers, points_per_fiber, transform):
    """Transform in_array of R,A,S by transform (15 components, rotation about
    R,A,S, translation in R, A, S,  scale along R, A, S, and
    shear. Fibers are assumed to be in RAS (or LPS as long as all inputs
    are consistent).  Transformed fibers are returned.

    This overall slows down computations. So keep all data in numpy double arrays for now.
    """
    out_array = np.zeros((3, number_of_fibers, points_per_fiber))

    vtktrans = convert_transform_to_vtk(transform, scaled=True)

    # Transform moving fiber array by applying transform to original fibers
    pts = vtk.vtkPoints()
    vtktrans.TransformPoints(moving_points, pts)

    # Copy this into numpy world
    count = 0
    for lidx in range(0, number_of_fibers):
        for pidx in range(0, points_per_fiber):
            (out_array[0, lidx, pidx], out_array[1, lidx, pidx], out_array[2, lidx, pidx]) = pts.GetPoint(count)
            ## pt = pts.GetPoint(count)
            ## out_array[0, lidx, pidx] = pt[0]
            ## out_array[1, lidx, pidx] = pt[1]
            ## out_array[2, lidx, pidx] = pt[2]
            count += 1

    del vtktrans
    return out_array

def convert_transform_to_vtk(transform, scaled=False, mode=[1,1,1,1]):
    """ Produce an output vtkTransform corresponding to the
    registration results. Input is a 15-component
    transform vector."""

    if scaled:
        # transform_scaling (must be same as defined above in class)
        #transform = np.divide(transform, np.array([1, 1, 1, .5, .5, .5, 300, 300, 300, 1, 1, 1, 1, 1, 1]))
        transform = np.divide(transform, np.array([1, 1, 1, .5, .5, .5, 200, 200, 200, 1, 1, 1, 1, 1, 1]))
                                 
    vtktrans = vtk.vtkTransform()

    if mode[0]:
        # translate first so optimizer starts there
        vtktrans.Translate(transform[0],
                        transform[1], transform[2])

    if mode[1]:
        # degrees
        vtktrans.RotateX(transform[3])
        vtktrans.RotateY(transform[4])
        vtktrans.RotateZ(transform[5])

    if mode[2]:
        vtktrans.Scale(transform[6],
                    transform[7], transform[8])

    if mode[3]:
        #// Update affine transformation: Add shearing
        #vtkMatrix4x4 *skewx= vtkMatrix4x4::New();  skewx->Identity();
        #skewx->SetElement(2, 1,tan(_szy*(pi/180.0)));   skewx->SetElement(1, 2,tan(_syz*(pi/180.0))); 
        #vtkMatrix4x4 *skewy= vtkMatrix4x4::New();   skewy->Identity();
        #skewy->SetElement(2, 0,tan(_szx*(pi/180.0)));   skewy->SetElement(0, 2,tan(_sxz*(pi/180.0)));
        #vtkMatrix4x4 *skewz= vtkMatrix4x4::New();   skewz->Identity();
        #skewz->SetElement(1, 0,tan(_sxy*(pi/180.0)));    skewz->SetElement(0, 1,tan(_syx*(pi/180.0))); 
        #tr->Concatenate(skewx);   tr->Concatenate(skewy);   tr->Concatenate(skewz);
        sxy = transform[9] * np.pi/180.0
        sxz = transform[10] * np.pi/180.0
        syx = transform[11] * np.pi/180.0
        syz = transform[12] * np.pi/180.0
        szx = transform[13] * np.pi/180.0
        szy = transform[14] * np.pi/180.0
        skewx = vtk.vtkMatrix4x4()
        skewy = vtk.vtkMatrix4x4()
        skewz = vtk.vtkMatrix4x4()
        skewx.SetElement(2, 1, np.tan(szy))
        skewx.SetElement(1, 2, np.tan(syz))
        skewy.SetElement(2, 0, np.tan(szx))
        skewy.SetElement(0, 2, np.tan(sxz))
        skewz.SetElement(1, 0, np.tan(sxy))
        skewz.SetElement(0, 1, np.tan(syx))
        vtktrans.Concatenate(skewx)
        vtktrans.Concatenate(skewy)
        vtktrans.Concatenate(skewz)

    return vtktrans


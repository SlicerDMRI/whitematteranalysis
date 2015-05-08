""" register.py

implementation of fiber tractography registration (group)

class RegisterTractography


"""

try:
    import scipy.optimize
    USE_SCIPY = 1
except ImportError:
    USE_SCIPY = 0
    print "<congeal.py> Failed to import  scipy.optimize, cannot align or register."
    print "<congeal.py> Please install  scipy.optimize for this functionality."

import numpy
import sys
import time
import vtk
try:
    from joblib import Parallel, delayed
    USE_PARALLEL = 1
except ImportError:
    USE_PARALLEL = 0
    print "<congeal.py> Failed to import joblib, cannot multiprocess."
    print "<congeal.py> Please install joblib for this functionality."

import whitematteranalysis as wma

class RegisterTractography:

    def constraint(self, x_current):
        penalty = 0
        
        # The values should be a registration matrix.
        # Make sure the optimizer is searching in a reasonable region.
        transform = numpy.divide(x_current,self.transform_scaling)
        
        # want translation reasonable, say -100 to 100 range
        max_trans = 100
        if numpy.abs(transform[0]) > max_trans:
            penalty -= 1
        if numpy.abs(transform[1]) > max_trans:
            penalty -= 1
        if numpy.abs(transform[2]) > max_trans:
            penalty -= 1
        # want rotation reasonable, say +/- 90 degrees is the max we could reasonably handle
        max_rot = 90
        if numpy.abs(transform[3]) > max_rot:
            penalty -= 1
        if numpy.abs(transform[4]) > max_rot:
            penalty -= 1
        if numpy.abs(transform[5]) > max_rot:
            penalty -= 1
        # want scaling between 0.9 and 1.2. In practice, only shrinking is an issue.
        # Don't let it shrink by more than 0.75
        min_scale = 0.75
        if transform[6] < min_scale:
            penalty -= 1
        if transform[7] < min_scale:
            penalty -= 1
        if transform[8] < min_scale:
            penalty -= 1
        # want shear angles reasonable like rotation
        #max_shear_angle = 45
        max_shear_angle = 90
        if numpy.abs(transform[9]) > max_shear_angle:
            penalty -= 1
        if numpy.abs(transform[10]) > max_shear_angle:
            penalty -= 1
        if numpy.abs(transform[11]) > max_shear_angle:
            penalty -= 1
        if numpy.abs(transform[12]) > max_shear_angle:
            penalty -= 1
        if numpy.abs(transform[13]) > max_shear_angle:
            penalty -= 1
        if numpy.abs(transform[14]) > max_shear_angle:
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

        # translate then rotate so optimizer will translate first
        # trans, rot, scale, shear:
        self.initial_transform = numpy.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0])

        # the scaling components should have much smaller changes than the others, so scale them accordingly in optimizer search space
        #self.transform_scaling = numpy.array([1, 1, 1, .5, .5, .5, 500, 500, 500, 2, 2, 2, 2, 2, 2])
        # try more flexible shear
        #self.transform_scaling = numpy.array([1, 1, 1, .5, .5, .5, 500, 500, 500,1, 1, 1, 1, 1, 1])
        # try more flexible scaling
        self.transform_scaling = numpy.array([1, 1, 1, .5, .5, .5, 300, 300, 300,  1, 1, 1, 1, 1, 1])

        # internal recordkeeping
        self._x_opt = None
        self.iterations = 0
        
    def objective_function(self, current_x):
        """ The actual objective used in registration.  Function of
        the current x in search space, as well as parameters of the
        class: threshold, sigma. Compares sampled fibers from moving
        input, to all fibers of fixed input."""

        # keep track of current x value
        self._x_opt = current_x

        # get and apply transforms from current_x
        moving = transform_fiber_array_numpy(self.moving, current_x)

        # compute objective
        obj = inner_loop_objective(self.fixed, moving, self.sigma * self.sigma)

        # save objective function value for analysis of performance
        self.objective_function_values.append(obj)

        if self.verbose:
            print "O:",  obj, "X:", self._x_opt
       #print "X:", self._x_opt
        return obj

    def compute(self):

        """ Run the registration.  Add subjects first (before calling
        compute). Then call compute several times, using different
        parameters for the class, for example first just for
        translation."""

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
            ren.save_views(self.output_directory, 'fixed_brain_' + self.process_id_string)
            del ren
                
        self.iterations += 1

        if self.verbose:
            print "<congeal.py> Initial value for X:", self.initial_transform


        # These optimizers were tested, less successfully than cobyla.
        # This code is left here as documentation.
        #scipy.optimize.fmin(self.objective_function,
        #                    self._x_opt, params, maxiter=self.maxiter,
        #                    ftol=self.ftol, maxfun=self.maxfun)
        #self.final_transform = scipy.optimize.fmin_powell(self.objective_function,
        #                                            numpy.multiply(self.initial_transform,self.transform_scaling),
        #                                            xtol=self.rhobeg,
        #                                            ftol = 0.05,
        #                                            maxfun=self.maxfun)
        #maxiter=5)

        # Optimize using cobyla. Allows definition of initial and
        # final step size scales (rhos), as well as constraints.  Here
        # we use the constraints to encourage that the transform stays a transform.
        # note disp 0 turns off all display
        self.final_transform = scipy.optimize.fmin_cobyla(self.objective_function,
                                                  numpy.multiply(self.initial_transform,self.transform_scaling), self.constraint,
                                                  maxfun=self.maxfun, rhobeg=self.initial_step,
                                                  rhoend=self.final_step, disp=1)
        self.final_transform = numpy.divide(self.final_transform,self.transform_scaling)

        print "O:", self.objective_function_values

        # Return output transforms from this iteration
        return self.final_transform


def inner_loop_objective(fixed, moving, sigmasq):
    """ The code called within the objective_function to find the
    negative log probability of one brain given all other brains. Used
    with Entropy objective function."""

    (dims, number_of_fibers_moving, points_per_fiber) = moving.shape
    # number of compared fibers (normalization factor)
    (dims, number_of_fibers_fixed, points_per_fiber) = fixed.shape

    probability = numpy.zeros(number_of_fibers_moving) + 1e-20
    
    # Loop over fibers in moving, find total probability of
    # fiber using all fibers from fixed.
    for idx in range(number_of_fibers_moving):
        probability[idx] += total_probability_numpy(moving[:,idx,:], fixed,
                sigmasq)

    # divide total probability by number of fiber comparisons
    # the output must be between 0 and 1
    number_comparisons = numpy.multiply(number_of_fibers_moving, number_of_fibers_fixed)
    #print "PROBABILITY:", numpy.min(probability), numpy.max(probability),
    probability /= number_comparisons
    #print numpy.min(probability), numpy.max(probability)
    # add negative log probabilities of all fibers in this brain.
    entropy = numpy.sum(- numpy.log(probability))
    #print "NUMBER COMPARISONS:", number_comparisons, number_of_fibers_moving, number_of_fibers_fixed
    return entropy

def total_probability_numpy(moving_fiber, fixed_fibers, sigmasq):
    distance = fiber_distance_numpy(moving_fiber, fixed_fibers)
    return numpy.sum(numpy.exp(-distance / (sigmasq)))

def fiber_distance_numpy(moving_fiber, fixed_fibers):
    """
    Find pairwise fiber distance from fiber to all fibers in fiber_array.
    The Mean and MeanSquared distances are the average distance per
    fiber point, to remove scaling effects (dependence on number of
    points chosen for fiber parameterization). The Hausdorff distance
    is the maximum distance between corresponding points.
    input fiber should be class Fiber. fibers should be class FiberArray
    """
    # compute pairwise fiber distances along fibers
    distance_1 = _fiber_distance_internal_use_numpy(moving_fiber, fixed_fibers)
    distance_2 = _fiber_distance_internal_use_numpy(moving_fiber, fixed_fibers,reverse_fiber_order=True)
        
    # choose the lowest distance, corresponding to the optimal fiber
    # representation (either forward or reverse order)
    return numpy.minimum(distance_1, distance_2)

def _fiber_distance_internal_use_numpy(moving_fiber, fixed_fibers, reverse_fiber_order=False):
    """ Compute the total fiber distance from one fiber to an array of
    many fibers.
    This function does not handle equivalent fiber representations,
    for that use fiber_distance, above.
    """
    #print "SHAPE:", fixed_fibers[0,:,:].shape, moving_fiber[0,::-1].shape
    
    # compute the distance from this fiber to the array of other fibers
    if reverse_fiber_order:
        ddx = fixed_fibers[0,:,:] - moving_fiber[0,::-1]
        ddy = fixed_fibers[1,:,:] - moving_fiber[1,::-1]
        ddz = fixed_fibers[2,:,:] - moving_fiber[2,::-1]
    else:
        ddx = fixed_fibers[0,:,:] - moving_fiber[0,:]
        ddy = fixed_fibers[1,:,:] - moving_fiber[1,:]
        ddz = fixed_fibers[2,:,:] - moving_fiber[2,:]

    #print "MAX abs ddx:", numpy.max(numpy.abs(ddx)), "MAX ddy:", numpy.max(numpy.abs(ddy)), "MAX ddz:", numpy.max(numpy.abs(ddz))
    #print "MIN abs ddx:", numpy.min(numpy.abs(ddx)), "MIN ddy:", numpy.min(numpy.abs(ddy)), "MIN ddz:", numpy.min(numpy.abs(ddz))
    
    dx = numpy.square(ddx)
    dy = numpy.square(ddy)
    dz = numpy.square(ddz)

    distance = dx + dy + dz

    #elif distance_method == 'Hausdorff':
    # take max along fiber
    #return numpy.max(distance, 1)
    # to test mean distance (had a less steep objective)
    return numpy.mean(distance, 1)

        

##__render_count = 0

def transform_fiber_array_numpy(in_array, transform):
    """Transform in_array of R,A,S by transform (15 components, rotation about
    R,A,S, translation in R, A, S,  scale along R, A, S, and
    shear. Fibers are assumed to be in RAS (or LPS as long as all inputs
    are consistent).  Transformed fibers are returned.
    """
    (dims, number_of_fibers, points_per_fiber) = in_array.shape
    out_array = numpy.zeros(in_array.shape)

    vtktrans = convert_transform_to_vtk(transform, scaled=True)
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
    ## if (numpy.mod(__render_count, 500) == 0) & False:
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


def convert_transform_to_vtk(transform, scaled=False, mode=[1,1,1,1]):
    """ Produce an output vtkTransform corresponding to the
    registration results. Input is a 15-component
    transform vector."""

    if scaled:
        # transform_scaling (must be same as defined above in class)
        transform = numpy.divide(transform, numpy.array([1, 1, 1, .5, .5, .5, 300, 300, 300, 1, 1, 1, 1, 1, 1]))
                                 
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
        sxy = transform[9] * numpy.pi/180.0
        sxz = transform[10] * numpy.pi/180.0
        syx = transform[11] * numpy.pi/180.0
        syz = transform[12] * numpy.pi/180.0
        szx = transform[13] * numpy.pi/180.0
        szy = transform[14] * numpy.pi/180.0
        skewx = vtk.vtkMatrix4x4()
        skewy = vtk.vtkMatrix4x4()
        skewz = vtk.vtkMatrix4x4()
        skewx.SetElement(2, 1, numpy.tan(szy))
        skewx.SetElement(1, 2, numpy.tan(syz))
        skewy.SetElement(2, 0, numpy.tan(szx))
        skewy.SetElement(0, 2, numpy.tan(sxz))
        skewz.SetElement(1, 0, numpy.tan(sxy))
        skewz.SetElement(0, 1, numpy.tan(syx))
        vtktrans.Concatenate(skewx)
        vtktrans.Concatenate(skewy)
        vtktrans.Concatenate(skewz)

    return vtktrans


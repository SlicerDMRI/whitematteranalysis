""" register.py

implementation of fiber tractography registration (group)

class RegisterTractographyNonlinear


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

class RegisterTractographyNonlinear(wma.register_two_subjects.RegisterTractography):

    def constraint(self, x_current):
        # Make sure the optimizer is searching in a reasonable region.
        # aim to preserve volume using the determinant of the overall affine transform
        affine_part = vtk.vtkLandmarkTransform()
        source_landmarks = convert_numpy_array_to_vtk_points(x_current)
        target_landmarks = convert_numpy_array_to_vtk_points(self.target_landmarks)
        affine_part.SetSourceLandmarks(source_landmarks)
        affine_part.SetTargetLandmarks(target_landmarks)
        affine_part.SetModeToAffine()
        affine_part.Update()
        det = affine_part.GetMatrix().Determinant()
        penalty = -100 * (1 - det)
        return penalty
    
    def __init__(self):
        # parameters that should be set by user
        self.sigma = 5
        self.sigma2 = 1.5
        self.process_id_string = ""
        self.output_directory = None
        
        # performance options that should be set by user
        self.verbose = False
        self.render = False
        
        # optimizer parameters that should be set by user
        self.maxfun = 300

        # output of registration
        self.objective_function_values = list()
        self.final_transform = None
        
        # subject data that must be input
        #self.fixed = None
        #self.moving = None
        #self.initial_step = 5
        #self.final_step = 2

        # set up default landmarks
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

        # transform we optimize over is the source landmarks (initialize to identity, equal to target landmarks)
        self.initial_transform = self.target_landmarks

        # internal recordkeeping
        self._x_opt = None
        self.iterations = 0
        
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
        self.target_points = convert_numpy_array_to_vtk_points(self.target_landmarks)

    def objective_function(self, current_x):
        """ The actual objective used in registration.  Function of
        the current x in search space, as well as parameters of the
        class: threshold, sigma. Compares sampled fibers from moving
        input, to all fibers of fixed input."""

        # keep track of current x value
        self._x_opt = current_x

        # get and apply transforms from current_x
        moving = self.transform_fiber_array_numpy(self.moving, current_x)

        # compute objective
        obj = inner_loop_objective(self.fixed, self.length_fixed, moving, [self.sigma * self.sigma, self.sigma2 * self.sigma2])

        # save objective function value for analysis of performance
        self.objective_function_values.append(obj)

        if self.verbose:
            print "O:",  obj, "X:", self._x_opt
        #print "X:", self._x_opt
        return obj

    def transform_fiber_array_numpy(self, in_array, source_landmarks):
        """Transform in_array of R,A,S by transform (a list of source points).  Transformed fibers are returned.
        """
        (dims, number_of_fibers, points_per_fiber) = in_array.shape
        out_array = numpy.zeros(in_array.shape)

        vtktrans = convert_transform_to_vtk(source_landmarks, self.target_points)
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

        # Find step length of all fibers in the fixed brain
        (dims, number_of_fibers_fixed, points_per_fiber) = self.fixed.shape
        ddx = self.fixed[0,:,0:-2] - self.fixed[0,:,1:-1]
        ddy = self.fixed[1,:,0:-2] - self.fixed[1,:,1:-1]
        ddz = self.fixed[2,:,0:-2] - self.fixed[2,:,1:-1]
        dx = numpy.square(ddx)
        dy = numpy.square(ddy)
        dz = numpy.square(ddz)
        distance = numpy.sqrt(dx + dy + dz)
        #self.length_fixed = numpy.divide(numpy.sum(distance, 1), points_per_fiber)
        self.length_fixed = numpy.sum(distance, 1)
        
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
                                                  self.initial_transform, self.constraint,
                                                  maxfun=self.maxfun, rhobeg=self.initial_step,
                                                  rhoend=self.final_step, disp=0)

        #self.final_transform = numpy.divide(self.final_transform,self.transform_scaling)

        print "O:", self.objective_function_values

        # Return output transforms from this iteration
        return self.final_transform


def inner_loop_objective(fixed, length_fixed, moving, sigmasq):
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
                sigmasq, length_fixed)

    # divide total probability by number of fibers in the atlas ("mean
    # brain").  this neglects Z, the normalization constant for the
    # pdf, which would not affect the optimization.
    probability /= number_of_fibers_fixed
    print "PROBABILITY:", numpy.min(probability), numpy.max(probability),
    # add negative log probabilities of all fibers in this brain.
    entropy = numpy.sum(- numpy.log(probability))
    return entropy

def total_probability_numpyOLD(moving_fiber, fixed_fibers, sigmasq):
    distance = fiber_distance_numpy(moving_fiber, fixed_fibers)
    # also allow alignment to a parallel but shorter fiber
    # try to decrease sensitivity to how far was scanned in brainstem and tracking differences
    #distance_flex = fiber_distance_numpy(moving_fiber[:,1:-1], fixed_fibers[:,:,1:-1])
    #distance = numpy.minimum(distance, distance_flex)
    # twice changes objective for the worse. once has little change, but a bit steeper translation.
    #distance_flex = fiber_distance_numpy(moving_fiber[:,2:-2], fixed_fibers[:,:,2:-2])
    #distance = numpy.minimum(distance, distance_flex)
    return numpy.sum(numpy.exp(-distance / (sigmasq)))

def total_probability_numpy(moving_fiber, fixed_fibers, sigmasq, length_fixed):
    distance = fiber_distance_numpy(moving_fiber, fixed_fibers)

    probability1 = numpy.exp(numpy.divide(-distance, (sigmasq[0])))
    
    # test to reduce shrinkage: compare step size between points of fixed vs moving brain
    #(dims, number_of_fibers_fixed, points_per_fiber) = fixed_fibers.shape
    (dims, points_per_fiber) = moving_fiber.shape
    ddx = moving_fiber[0,0:-2] - moving_fiber[0,1:-1]
    ddy = moving_fiber[1,0:-2] - moving_fiber[1,1:-1]
    ddz = moving_fiber[2,0:-2] - moving_fiber[2,1:-1]
    dx = numpy.square(ddx)
    dy = numpy.square(ddy)
    dz = numpy.square(ddz)
    dists = numpy.sqrt(dx + dy + dz)
    #step_length_moving = numpy.divide(numpy.sum(dists), points_per_fiber)
    length_moving = numpy.sum(dists)

    diffs = length_moving - length_fixed
    probability2 = numpy.exp(-numpy.divide(numpy.square(length_moving - length_fixed), sigmasq[1]))

    probability = numpy.multiply(probability1, probability2)
    #print "SIGMAS:", numpy.sqrt(sigmasq), "LENGTHS:", numpy.min(length_fixed), numpy.max(length_fixed), length_moving, "DIFFS:", numpy.min(diffs), numpy.max(diffs), "PROBABILITY 2:", numpy.min(probability2), numpy.max(probability2), probability2.shape, "DISTANCES:", numpy.min(numpy.sqrt(distance)), numpy.max(numpy.sqrt(distance)), distance.shape, "PROBABILITY 1:", numpy.min(probability1), numpy.max(probability1), probability1.shape, "PROBABILITY:", numpy.min(probability), numpy.max(probability), probability.shape, "FIBERS:", moving_fiber.shape, fixed_fibers.shape, length_fixed.shape

    return numpy.sum(probability)

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

    # to test mean distance (had a less steep objective but looks similar)
    return numpy.mean(distance, 1)

    #elif distance_method == 'Hausdorff':
    # take max along fiber
    #return numpy.max(distance, 1)

def convert_numpy_array_to_vtk_points(inarray):
    number_of_points = len(inarray)/3
    vtk_points = vtk.vtkPoints()
    vtk_points.SetNumberOfPoints(number_of_points)
    idx = 0
    for pt in zip(inarray[::3], inarray[1::3], inarray[2::3]):
        #print pt
        vtk_points.SetPoint(idx, pt[0], pt[1], pt[2])
        idx += 1
    return vtk_points

def convert_transform_to_vtk(source_landmarks, target_points):
    """Produce an output vtkThinPlateSplineTransform corresponding to the

    registration results. Input is a numpy array of of source (moving)
    landmarks and a vtkPoints object of target (fixed) landmarks.
    """

    source_points = convert_numpy_array_to_vtk_points(source_landmarks)
    #number_of_points = len(source_landmarks)/3
    #print "CONVERT:", len(source_landmarks), number_of_points, target_points.GetNumberOfPoints(), source_points.GetNumberOfPoints()

    return compute_thin_plate_spline_transform(source_points, target_points)

def compute_thin_plate_spline_transform(source_points, target_points):
    """Produce an output vtkThinPlateSplineTransform.
    Input is a vtkPoints object of source (moving) landmarks and a vtkPoints
    object of target (fixed) landmarks.
    """
    vtktrans = vtk.vtkThinPlateSplineTransform()
    vtktrans.SetSourceLandmarks(source_points)
    vtktrans.SetTargetLandmarks(target_points)
    vtktrans.SetBasisToR()
    
    return vtktrans


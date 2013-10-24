""" register.py

implementation of fiber tractography registration (pairwise)

class RegisterTractography


"""

try:
    import scipy.optimize
    USE_SCIPY = 1
except ImportError:
    USE_SCIPY = 0
    print "<register.py> Failed to import  scipy.optimize, cannot align or register."
    print "<register.py> Please install  scipy.optimize for this functionality."

import numpy

import vtk
try:
    from joblib import Parallel, delayed
    USE_PARALLEL = 1
except ImportError:
    USE_PARALLEL = 0
    print "<register.py> Failed to import joblib, cannot multiprocess."
    print "<register.py> Please install joblib for this functionality."


import whitematteranalysis.fibers
import whitematteranalysis.similarity


class RegistrationInformation:
    def __init__(self):
        self._original_fibers = whitematteranalysis.fibers.FiberArray()
        self._moving_fibers = whitematteranalysis.fibers.FiberArray()
        # indices of moving fibers to compute the objective function
        self._moving_fiber_sample = []
        self.points_per_fiber = 5
        self.fiber_sample_size = 200

        # transformation matrices for internal use
        # (vtkTransform is returned by compute) 
        # rot x,y,z trans x,y,z scale x,y,z
        #self.transform = numpy.array([0, 0, 0, 0, 0, 0, 1, 1, 1]).astype(float)
        self.transform = numpy.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0]).astype(float)
        self.modified = False
        self.x = []

    def initialize(self, polydata):
        # internal representation for fast similarity computation
        self._original_fibers.convert_from_polydata(polydata,
                                                self.points_per_fiber)
        #self._moving_fibers.convert_from_polydata(movingpd,
        #                                         self.points_per_fiber)
        self.modified = True
        self.apply_transform()

    def initialize_fiber_sample(self):
        # indices of moving fibers to compute the objective function
        self._moving_fiber_sample = numpy.random.random_integers(
            0, self._original_fibers.number_of_fibers - 1,
            self.fiber_sample_size)
        self.modified = True

    def apply_transform(self):
        # apply transform to moving fiber data IF the transform is modified
        if self.modified:
            self._moving_fibers = \
                self.transform_fiber_array(self._original_fibers,
                                           self.transform)

    def transform_fiber_array(self, in_array, transform):
        """Transform in_array (of class FiberArray) by transform (9
        components, rotation about R,A,S, translation in R, A, S, and
        scale along R, A, S. Fibers are assumed to be in RAS.  Transformed
        fibers are returned."""
        
        out_array = whitematteranalysis.fibers.FiberArray()
        pd_in = in_array.convert_to_polydata()
        transformer = vtk.vtkTransformPolyDataFilter()
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                transformer.SetInputData(pd_in)
            else:
                transformer.SetInput(pd_in) 
        vtktrans = self.convert_transform_to_vtk(transform)
        transformer.SetTransform(vtktrans)
        transformer.Update()
        pd_out = transformer.GetOutput()
        out_array.convert_from_polydata(pd_out)
        return out_array

    def set_transform(self, input_transform):
        input_transform = numpy.array(input_transform)
        # decide whether transform was modified
        if numpy.max(numpy.abs(self.transform - input_transform)) == 0.0:
            self.modified = False
        else:
            # directly set it. assume it has 9 components 
            self.transform = numpy.copy(input_transform)
            self.modified = True
        #print "Set_transform, modified:", self.modified, "T0:", self.transform, "T1:", input_transform

    def convert_transform_to_vtk(self, transform=None):
        """ Produce an output vtkTransform corresponding to the
        registration results. Optionally can input a 9-component
        transform vector."""
        
        if transform is None:
            transform = self.transform
        
        vtktrans = vtk.vtkTransform()

        vtktrans.RotateX(transform[0] * (180 / numpy.pi))
        vtktrans.RotateY(transform[1] * (180 / numpy.pi))
        vtktrans.RotateZ(transform[2] * (180 / numpy.pi))

        vtktrans.Translate(transform[3],
                           transform[4], transform[5])

        vtktrans.Scale(transform[6],
                       transform[7], transform[8])

        #// Update affine transformation: Add shearing
        #vtkMatrix4x4 *skewx= vtkMatrix4x4::New();  skewx->Identity();
        #skewx->SetElement(2, 1,tan(_szy*(pi/180.0)));   skewx->SetElement(1, 2,tan(_syz*(pi/180.0))); 
        #vtkMatrix4x4 *skewy= vtkMatrix4x4::New();   skewy->Identity();
        #skewy->SetElement(2, 0,tan(_szx*(pi/180.0)));   skewy->SetElement(0, 2,tan(_sxz*(pi/180.0)));
        #vtkMatrix4x4 *skewz= vtkMatrix4x4::New();   skewz->Identity();
        #skewz->SetElement(1, 0,tan(_sxy*(pi/180.0)));    skewz->SetElement(0, 1,tan(_syx*(pi/180.0))); 
        #tr->Concatenate(skewx);   tr->Concatenate(skewy);   tr->Concatenate(skewz);
        sxy = transform[9]
        sxz = transform[10]
        syx = transform[11]
        syz = transform[12]
        szx = transform[13]
        szy = transform[14]
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
        #del skewx
        #del skewy
        #del skewz
        return vtktrans

class RegisterTractography:
    """ Input 2 polydata using initialize method. Then run compute
    method to perform registration. Compute may be run several times,
    and the state (transform matrix) is preserved between runs. This
    allows coarse-to-fine registration strategy, using the sigma and
    threshold parameters for fiber distance and conversion to
    similarity, as well as the fiber_sample_size parameter that can be
    increased from a low value (i.e. 200) to a reasonable value for
    final tuning (>=1000)."""
    
    def translate_only(self):
        """ Set registration mode to only translation. """
        self._registration_mode = "TranslateOnly"

    def rotate_only(self):
        """ Set registration mode to only rotation. """
        self._registration_mode = "RotateOnly"

    def translate_and_rotate(self):
        """ Set registration mode to both translation and rotation. """
        self._registration_mode = "TranslateAndRotate"

    def translate_and_rotate_and_scale(self): 
        """ Set registration mode to 9 dof, translate, rotate, and scale. """
        self._registration_mode = "TranslateRotateScale"

    def __init__(self):
        # parameters that can be set by user
        self.sigma = 10
        self.threshold = 5
        self.points_per_fiber = 5
        self.fiber_sample_size = 200

        # performance options set by user
        #self.verbose = 0
        self.parallel_verbose = 0
        self.parallel_jobs = 4

        # optimizer parameters set by user
        #self.maxiter = 30
        #self.ftol = 0.0001
        self.maxfun = 300

        # output of registration
        self.objective_function_values = []
        # transformation matrices for internal use
        # (vtkTransform is returned by compute) 
        # rot x,y,z trans x,y,z scale x,y,z
        self._transform = [0, 0, 0, 0, 0, 0, 1, 1, 1]

        # internal parameter, please set with functions above
        self._registration_mode = "TranslateOnly"
        
        # internal data storage
        self._fixed_fibers = whitematteranalysis.fibers.FiberArray()
        self._moving_fibers = whitematteranalysis.fibers.FiberArray()
        self._tx_fibers = whitematteranalysis.fibers.FiberArray()
        # indices of moving fibers to compute the objective function
        self._moving_fiber_sample = []

        # square sigma for use in objective function
        self._sigmasq = self.sigma * self.sigma

        # internal variables for optimization
        # We want dimensions of radians and mm search space
        # comparable, so that increments taken by optimizer are
        # similar in each direction.  Set scaling of x according to
        # maximum expected transform values: 30 degrees around any
        # axis, 20 mm translation, and scale factor delta of .1
        #maxrot = (60.0 / 180.0) * numpy.pi
        #maxtrans = 20.0
        ##self.maxscale = .05
        #maxscale = 1.1
        #self._x_scaling = [maxrot, maxrot, maxrot,
        #                 maxtrans, maxtrans, maxtrans,
        #                 maxscale, maxscale, maxscale]

        # initial increments: 5mm, 10 degrees, .05 scale
        # these should each be a step of 1 in search space       
        #inc_rot = (20.0 / 180.0) * numpy.pi
        #inc_trans = 10.0
        #inc_scale = 1
        #self._x_scaling = [inc_rot, inc_rot, inc_rot,
        #                 inc_trans, inc_trans, inc_trans,
        #                 inc_scale, inc_scale, inc_scale]

        # make search space "larger"
        #self._x_scaling = numpy.divide(self._x_scaling, 100)

        # cobyla
        self._x_scaling = [1, 1, 1, 1, 1, 1, 1, 1, 1]

        inc_rot = (30.0 / 180.0) * numpy.pi
        inc_trans = 20.0
        inc_scale = .2

        self.rhobeg = [inc_rot, inc_rot, inc_rot,
                       inc_trans, inc_trans, inc_trans,
                       inc_scale, inc_scale, inc_scale]

        inc_rot = (0.1 / 180.0) * numpy.pi
        inc_trans = 0.01
        inc_scale = 0.001

        self.rhoend = [inc_rot, inc_rot, inc_rot,
                       inc_trans, inc_trans, inc_trans,
                       inc_scale, inc_scale, inc_scale]

        # this takes more reasonable sized steps when scaled, so it
        # converges faster (for fmin only)
        #self._x_scaling = numpy.multiply(self._x_scaling, 1000)

        # x is the quantity we are optimizing.
        self._x_opt = self.convert_transform_to_x(self._transform)

    def convert_x_to_transform(self, input_x):
        """ Convert the input_x (in optimizer search space
        coordinates) to a transform. Pays attention to the current
        registration mode and only modifies the corresponding
        components of the transform, allowing preservation of previous
        registration information."""
        
        transform = self._transform

        # only modify the part of the transform that x represents
        # since we may have a previous/current/initial transform set

        if self._registration_mode == "TranslateOnly":
            transform[3:6] = numpy.multiply(input_x[0:3], self._x_scaling[3:6])
        elif self._registration_mode == "RotateOnly":
            transform[0:3] = numpy.multiply(input_x[0:3], self._x_scaling[0:3])
        elif self._registration_mode == "TranslateAndRotate":
            transform[0:6] = numpy.multiply(input_x, self._x_scaling[0:6])
        elif self._registration_mode == "TranslateRotateScale":
            transform = numpy.multiply(input_x, self._x_scaling)

        return transform
    
    def convert_transform_to_x(self, transform=None):
        """ Convert a transform (9 components) to the search space (x)
        coordinates."""
        
        if transform is None:
            transform = self._transform

        output_x = []
        rhobeg = []
        rhoend = []

        if self._registration_mode == "TranslateOnly":
            output_x = numpy.divide(transform[3:6], self._x_scaling[3:6])
            rhobeg = self.rhobeg[3:6]
            rhoend = self.rhoend[3:6]
        elif self._registration_mode == "RotateOnly":
            output_x = numpy.divide(transform[0:3], self._x_scaling[0:3])
            rhobeg = self.rhobeg[0:3]
            rhoend = self.rhoend[0:3]
        elif self._registration_mode == "TranslateAndRotate":
            output_x = numpy.divide(transform[0:6], self._x_scaling[0:6])
            rhobeg = self.rhobeg[0:6]
            rhoend = self.rhoend[0:6]
        elif self._registration_mode == "TranslateRotateScale":
            output_x = numpy.divide(transform, self._x_scaling)
            rhobeg = self.rhobeg[0:9]
            rhoend = self.rhoend[0:9]
 
        return output_x, rhobeg, rhoend

    def convert_transform_to_vtk(self, transform=None):
        """ Produce an output vtkTransform corresponding to the
        registration results. Optionally can input a 9-component
        transform vector."""
        
        if transform is None:
            transform = self._transform
        
        output_vtk_transform = vtk.vtkTransform()

        output_vtk_transform.RotateX(transform[0] * (180 / numpy.pi))
        output_vtk_transform.RotateY(transform[1] * (180 / numpy.pi))
        output_vtk_transform.RotateZ(transform[2] * (180 / numpy.pi))

        output_vtk_transform.Translate(transform[3],
                                     transform[4], transform[5])

        output_vtk_transform.Scale(transform[6],
                                 transform[7], transform[8])

        return output_vtk_transform
    
    def objective_function(self, current_x):
        """ The actual objective used in registration.  Function of
        the current x in search space, as well as parameters of the
        class: threshold, sigma. Compares sampled fibers from moving
        input, to all fibers of fixed input."""
        
        # keep track of current x value
        self._x_opt = current_x

        # get transform from current_x
        # x has units that better scale search space
        transform = self.convert_x_to_transform(self._x_opt)
        self._transform = transform
        
        # apply transform to moving fiber data
        # rotate and translate fiber Array
        # start with original array and output transformed version in tx_fibers
        whitematteranalysis.fibers.transform_fiber_array(self._moving_fibers,
                                                         self._tx_fibers,
                                                         transform)

        # compute total similarity from each moving fiber to all fixed fibers
        similarity1 = Parallel(
            n_jobs=self.parallel_jobs, verbose=self.parallel_verbose)(
            delayed(whitematteranalysis.similarity.total_similarity)(
                self._tx_fibers.get_fiber(lidx),
                self._fixed_fibers,
                self.threshold,
                self._sigmasq)
            for lidx in self._moving_fiber_sample)

        # other things that were tested:
        # compute total similarity from each fixed fiber to all moving fibers
        # combine the two directions of similarity
        #similarity1s =  numpy.sum(similarity1) / self.numFixedLines
        #similarity2s =  numpy.sum(similarity2) / self.numMovingLines
        #obj = 100000.0 / (similarity1s + similarity2s)
        # combine the two directions of similarity
        ###epsilon = 0.01
        ###similarity1s =  numpy.product(  numpy.add( similarity1 , epsilon))
        ###similarity2s =  numpy.product(  numpy.add( similarity2 , epsilon))
        ###obj = 1 / (similarity1s * similarity2s * 1e100)
        #obj = obj * numpy.abs(similarity1s-similarity2s)

        # divide total similarity by number of fiber comparisons
        # in order to have consistent values across runs (unless sigma changes)
        total_similarity = 1000 * numpy.sum(similarity1) / (
            self.fiber_sample_size * self._fixed_fibers.number_of_fibers)

        # return total similarity as objective function
        # since we want to minimize, return 1 over total similarity
        # this total is in the thousands so 1000/total is clearer
        #obj = 100000.0 / numpy.sum(similarity1)
        obj = 1 / total_similarity

        # save it for analysis of performance
        self.objective_function_values.append(obj)

        print "S:",  total_similarity
        print "O:",  obj
        print "            R:", self._transform[0:3]
        print "            T:", self._transform[3:6]
        print "            S:", self._transform[6:9]
        print "            X:", self._x_opt
        return obj

    def initialize(self, fixedpd, movingpd):
        """ Set input to the class. Do this once before starting
        multiple compute cycles."""
        
        # internal representation for fast similarity computation
        self._fixed_fibers.convert_from_polydata(fixedpd,
                                                self.points_per_fiber)
        self._moving_fibers.convert_from_polydata(movingpd,
                                                 self.points_per_fiber)

        # tell user we are doing something
        print "Input number of lines (fixed/moving):", \
              self._fixed_fibers.number_of_fibers, \
              self._moving_fibers.number_of_fibers, \
              "line length:", self.points_per_fiber
        print "SIGMA:", self.sigma
        print "Starting registration..."

        # initialize optimization using initial transform
        self._x_opt = self.convert_transform_to_x(self._transform)

        # initialize output of convergence info
        self.objective_function_values = []

        # indices of moving fibers to compute the objective function
        self._moving_fiber_sample = numpy.random.random_integers(
            0, self._moving_fibers.number_of_fibers - 1,
            self.fiber_sample_size)

    def compute_current_objective(self):
        """ For testing. Returns the value of the objective function,
        given the current state of the class."""
        
        # initialize optimization  using current transform
        self._x_opt = self.convert_transform_to_x(self._transform)

        # compute objective function once
        obj = self.objective_function(self._x_opt)

        return(obj)

    # temporary version. really should check bounds on transformation
    def constraint(self, current_x):
        """ Constraint function for cobyla optimizer. Currently does nothing."""
        # convert current_x to transform and check is within reasonable bounds
        # then return results
        return 1

    def compute(self):
        """ Run the registration.  Call initialize first. Then call
        compute several times, using different parameters for the
        class, for example first just for translation."""
        
        # initialize optimization  using current transform
        self._x_opt, rhobeg, rhoend = self.convert_transform_to_x(self._transform)

        # square current sigma parameter for later Gaussian
        self._sigmasq = self.sigma * self.sigma

        print "Initial values: T: ", self._transform, "X:", self._x_opt

        #scipy.optimize.fmin(self.objective_function,
        #                    self._x_opt, params, maxiter=self.maxiter,
        #                    ftol=self.ftol, maxfun=self.maxfun)

        scipy.optimize.fmin_cobyla(self.objective_function,
                                   self._x_opt, self.constraint,
                                   maxfun=self.maxfun, rhobeg=rhobeg,
                                   rhoend=rhoend)

        return self.convert_transform_to_vtk()

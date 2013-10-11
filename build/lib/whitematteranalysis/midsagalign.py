""" midsagalign.py

Automatic midsagittal alignment of tractography.

Works by optimizing an objective function over 3 parameters: rotation
about A-P axis, rotation about S-I axis, and translation along L-R
axis. Other transformations do not affect symmetry about the
midsagittal plane.

The objective function relates to the total fiber similarity of the
whole-brain tractography with its reflection.


"""

import numpy
try:
    import scipy.optimize
    USE_SCIPY = 1
except ImportError:
    USE_SCIPY = 0
    print "<midsagalign.py> Failed to import  scipy.optimize, cannot align or register."
    print "<midsagalign.py> Please install  scipy.optimize for this functionality."


import vtk
try:
    from joblib import Parallel, delayed
    USE_PARALLEL = 1
except ImportError:
    USE_PARALLEL = 0
    print "<midsagalign.py> Failed to import joblib, cannot multiprocess."
    print "<midsagalign.py> Please install joblib for this functionality."


import whitematteranalysis.similarity
import whitematteranalysis.register

class MidsagittalAlignment:

    """Optimize midsagittal plane by aligning tractography with its
    mirror image"""

    def __init__(self):
        # parameters that can be set by user
        self.sigma = 3
        self.threshold = 0
        self.points_per_fiber = 5
        self.fiber_sample_size = 300

        # performance options
        self.verbose = 1
        self.parallel_verbose = 0
        self.parallel_jobs = 1

        # data storage
        self.subject = whitematteranalysis.register.RegistrationInformation()
        self.objective_function_values = []

        # internal variables for optimization x has units of fractions
        # of max expected values, in order to better scale search
        # space
        self.x_opt = [0, 0, 0]

        # internal variables for optimization
        # cobyla
        #self._x_scaling = [1, 1, 1]
        #self.maxfun = 100
        
        #inc_rot = (10.0 / 180.0) * numpy.pi
        #inc_trans = 10.0
        #self.rhobeg = [inc_rot, inc_rot, inc_trans]

        #inc_rot = (0.5 / 180.0) * numpy.pi
        #inc_trans = 0.5
        #self.rhoend = [inc_rot, inc_rot, inc_trans]        
        # for optimization want dimensions of radians and mm search
        # space comparable max expected transform values are 10
        # degrees around y (A-P) axis, 20 degrees around z (S-I) axis,
        # and 10 mm shift L-R
        #self.xScaling = [(10/180.0)*numpy.pi, (20/180.0)*numpy.pi, 10]
        #self._x_scaling = [(5 / 180.0) * numpy.pi, (20 / 180.0) * numpy.pi, 10]
        self._x_scaling = [(10 / 180.0) * numpy.pi, (20 / 180.0) * numpy.pi, 10]
        
        # this takes more reasonable sized steps when scaled, so it
        # converges faster
        self._x_scaling = numpy.multiply(self._x_scaling, 1000)

    def objective_function(self, x_opt):

        """ 1 / total fiber similarity

        Similarity is computed from each fiber's reflection to all
        other fibers.  Reflection is across current midsagittal (Y-Z,
        i.e. AP-IS) plane.  A subset of fibers is used for reflection.

        """
        midsag_trans = numpy.multiply(x_opt, self._x_scaling)
        # rotation, rotation, translation
        self.subject.transform[1] = midsag_trans[0]
        self.subject.transform[2] = midsag_trans[1]
        self.subject.transform[3] = midsag_trans[2]

        # apply transform to data
        self.subject.apply_transform()

        # allocate outputs
        fiber_similarity = numpy.zeros(self.subject._original_fibers.number_of_fibers)

        # compute similarity from each fiber's reflection to all other fibers
        #fiber_similarity = Parallel(
        #    n_jobs=self.parallel_jobs, verbose=self.parallel_verbose)(
        #    delayed(whitematteranalysis.similarity.total_similarity)(
        #            self.subject._moving_fibers.get_fiber(lidx).get_reflected_fiber(), 
        #            self.subject._moving_fibers, 
        #            self.threshold, 
        #            self.sigmasq) 
        #     for lidx in self.subject._moving_fiber_sample)

        for lidx in self.subject._moving_fiber_sample:
                         fiber_similarity[lidx] = \
                             whitematteranalysis.similarity.total_similarity(
                                 self.subject._moving_fibers.get_fiber(lidx).get_reflected_fiber(), 
                                 self.subject._moving_fibers, 
                                 self.threshold, 
                                 self.sigmasq)
        
        # return total similarity as objective function. since we want
        # to minimize, return 1 over total similarity. this total is in
        # the thousands so 1000/total is a nicer number to look at
        obj = 100000.0 / numpy.sum(fiber_similarity)

        # save it for analysis of performance
        self.objective_function_values.append(obj)

        if self.verbose:
            print "OBJECTIVE:",  obj, \
                " T:", \
                midsag_trans[0] * 180 / numpy.pi, \
                midsag_trans[1] * 180 / numpy.pi, \
                midsag_trans[2]

        return obj

    def convert_transform_to_vtk(self):
        """ Convert self.transform to vtkTransform object, and return
        that object. """

        return self.subject.convert_transform_to_vtk()
    
    # temporary version. really should check bounds on transformation
    def constraint(self, current_x):
        """ Constraint function for cobyla optimizer. Currently does nothing."""

        return 1
    
    def compute(self, input_vtk_polydata):
        """ Perform the optimization, i.e. run the computation. """

        # internal representation for fast similarity computation
        # this also detects which hemisphere fibers are in
        self.subject.points_per_fiber = self.points_per_fiber
        self.subject.initialize(input_vtk_polydata)
        self.subject.fiber_sample_size =  self.fiber_sample_size
        self.subject.initialize_fiber_sample()

        # square sigma for later Gaussian
        self.sigmasq = self.sigma * self.sigma

        # tell user we are doing something
        if self.verbose:
            print "<midsagalign.py> Fibers in each hemisphere.", \
                "L:", self.subject._original_fibers.number_left_hem, \
                "R:", self.subject._original_fibers.number_right_hem, \
                "/ Total:", self.subject._original_fibers.number_of_fibers
            print "<midsagalign.py> Starting to optimize alignment."

        # run the computation, either in parallel or not
        scipy.optimize.fmin(self.objective_function, self.x_opt)

        #self.x_final = scipy.optimize.fmin_cobyla(self.objective_function,
        #                           self.x_opt, self.constraint,
        #                           maxfun=self.maxfun, rhobeg=self.rhobeg,
        #                           rhoend=self.rhoend)

        if self.verbose:
            print self.subject.transform

        self.vtktransform = self.convert_transform_to_vtk()

        transformer = vtk.vtkTransformPolyDataFilter()
        transformer.SetInput(input_vtk_polydata)
        transformer.SetTransform(self.vtktransform)
        transformer.Update()
        outpd = transformer.GetOutput()

        return outpd, self.vtktransform

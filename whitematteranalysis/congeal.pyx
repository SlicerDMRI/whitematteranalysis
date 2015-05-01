""" register.py

implementation of fiber tractography registration (group)

class CongealTractography


"""


try:
    import scipy.optimize
    USE_SCIPY = 1
except ImportError:
    USE_SCIPY = 0
    print "<congeal.py> Failed to import  scipy.optimize, cannot align or register."
    print "<congeal.py> Please install  scipy.optimize for this functionality."

import numpy

import vtk
try:
    from joblib import Parallel, delayed
    USE_PARALLEL = 1
except ImportError:
    USE_PARALLEL = 0
    print "<congeal.py> Failed to import joblib, cannot multiprocess."
    print "<congeal.py> Please install joblib for this functionality."


import whitematteranalysis.fibers
import whitematteranalysis.similarity
import whitematteranalysis.register
import whitematteranalysis.cluster


class CongealTractography:


    def translate_only(self):
        """ Set registration mode to only translation. """
        self._registration_mode = "TranslateOnly"

    def rotate_only(self):
        """ Set registration mode to only rotation. """
        self._registration_mode = "RotateOnly"

    def scale_only(self):
        """ Set registration mode to only rotation. """
        self._registration_mode = "ScaleOnly"

    def shear_only(self):
        """ Set registration mode to only rotation. """
        self._registration_mode = "ShearOnly"
        
    def translate_and_rotate(self):
        """ Set registration mode to both translation and rotation. """
        self._registration_mode = "TranslateAndRotate"

    def translate_and_rotate_and_scale(self):
        """ Set registration mode to 9 dof, translate, rotate, and scale. """
        self._registration_mode = "TranslateRotateScale"

    def use_total_fiber_similarity_objective(self):
        """ Set objective function mode to use fiber similarity total. This was an initial experiment and is not used."""
        self._objective_function_mode = "TotalFiberSimilarity"

    def use_entropy_objective(self):
        """ Set objective function mode to entropy. This was tested and published and is default."""
        self._objective_function_mode = "Entropy"

    def set_rhobeg(self, rot, trans, scale, shear):
        self.rhobeg = [rot, rot, rot, trans, trans, trans, scale, scale, scale, shear, shear, shear, shear, shear, shear]

    def set_rhoend(self, rot, trans, scale, shear):
        self.rhoend = [rot, rot, rot, trans, trans, trans, scale, scale, scale, shear, shear, shear, shear, shear, shear]

    def __init__(self):
        # parameters that can be set by user
        self.sigma = 10
        self.threshold = 0
        self.points_per_fiber = 5
        self.fiber_sample_size = 200
        #self.distance_method = 'MeanSquared'
        self.distance_method = 'Hausdorff'
        
        # performance options set by user
        self.verbose = 0
        self.parallel_verbose = 0
        self.parallel_jobs = 4

        # optimizer parameters set by user
        self.maxfun = 300

        # output of registration
        self.objective_function_values = []
        self.vtk_transform_list = []

        # internal parameters, please set with functions above
        self._registration_mode = "TranslateOnly"
        self._objective_function_mode = "Entropy"

        # squared sigma for use in objective function
        self._sigmasq = self.sigma * self.sigma

        # internal variables for optimization
        # cobyla
        self._x_scaling = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        #inc_rot = (10.0 / 180.0) * numpy.pi #larger than needed
        #inc_trans = 10.0 # this was too large looking at objective function plot
        #inc_scale = .10 # this was too large
        inc_rot = (5.0 / 180.0) * numpy.pi
        inc_trans = 5.0
        inc_scale = 0.01
        inc_shear = (2.0 / 180.0) * numpy.pi
        self.rhobeg = [inc_rot, inc_rot, inc_rot,
                       inc_trans, inc_trans, inc_trans,
                       inc_scale, inc_scale, inc_scale,
                       inc_shear, inc_shear, inc_shear,
                       inc_shear, inc_shear, inc_shear]

        inc_rot = (3 / 180.0) * numpy.pi
        inc_trans = 2.0
        inc_scale = 0.005
        # this was taking longer to converge than other optimizations
        # so it was too fine of a change.
        #inc_shear = (0.5 / 180.0) * numpy.pi
        inc_shear = (1.0 / 180.0) * numpy.pi
        self.rhoend = [inc_rot, inc_rot, inc_rot,
                       inc_trans, inc_trans, inc_trans,
                       inc_scale, inc_scale, inc_scale,
                       inc_shear, inc_shear, inc_shear,
                       inc_shear, inc_shear, inc_shear]
        
        # x is the quantity we are optimizing.
        # It's initialized by self.convert_transform_to_x.
        self._x_opt = []

        # subject data to be input
        self._subjects = list()
        self.number_of_subjects = 0
                
    def add_subject(self, polydata):
        subj = whitematteranalysis.register.RegistrationInformation()
        self._subjects.append(subj)
        subj.points_per_fiber = self.points_per_fiber
        subj.fiber_sample_size = self.fiber_sample_size
        subj.initialize(polydata)
        self.number_of_subjects = len(self._subjects)
	
    def initialize_fiber_samples(self):
        for subj in self._subjects:
            subj.fiber_sample_size = self.fiber_sample_size
            subj.initialize_fiber_sample()

    def convert_xs_to_transforms(self):
        """ For each part of x, corresponding to a subject, convert to
        transform. This code is called each time that the objective
        function is calculated. Mean rotation/translation/scale are removed."""

        if  numpy.isnan(self._x_opt).any():
            print "<congeal.py> ISNAN (conv vs to trans)", numpy.isnan(self._x_opt)

        incr = len(self._x_opt)/self.number_of_subjects

        idx = 0

        transform_list = list()
        
        for subj in self._subjects:
            subj.x = self._x_opt[idx : idx + incr]
            subj.set_transform(self.convert_x_to_transform(subj.x, subj.transform))
            idx = idx + incr
            transform_list.append(subj.transform)
            
        transforms_array = numpy.array(transform_list)

        if self.verbose:
            print "<congeal.py> TRANSFORMS"
            print numpy.round(transforms_array * 100) / 100
            meantrans = numpy.mean(transforms_array, 0)
            print "<congeal.py> Mean transform:"
            print numpy.round(meantrans * 1000) / 1000        

    def remove_mean_from_transforms(self):
        """ Remove mean rotations and mean scaling and mean
         translations from transforms.  the mean rotation and
         translation should not affect the objective function anyway.
         A mean scaling will affect the objective, i.e. if the brain
         shrinks all distances become smaller and the similarity is
         higher. Definitely not the desired effect."""
        
        meantrans = numpy.zeros(15)
        #array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        transform_list = list()
        for subj in self._subjects:
            transform_list.append(subj.transform)

        transforms_array = numpy.array(transform_list)
        meantrans = numpy.mean(transforms_array, 0)

        if self.verbose:
            print "<congeal.py> TRANSFORMS"
            print numpy.round(transforms_array * 100) / 100
            print "<congeal.py> Removing current (accumulated) mean transform before computing objective:"
            print numpy.round(meantrans * 1000) / 1000        

        for subj in self._subjects:
            subj.transform[0:6] = subj.transform[0:6] - meantrans[0:6]
            subj.transform[6:9] = subj.transform[6:9] / meantrans[6:9]
            subj.transform[9:15] = subj.transform[9:15] - meantrans[9:15]
            
    def convert_transforms_to_xs(self):
        """ For each part of x, corresponding to a subject, convert to
        transform. This code is called at the beginning of each
        compute() cycle to initialize x from the existing transforms."""

        x_opt_list = list()

        for subj in self._subjects:
            subj.x = self.convert_transform_to_x(subj.transform)
            x_opt_list.append(subj.x)

        self._x_opt = numpy.array(x_opt_list)

        # make it into a vector rather than a 2D array
        self._x_opt = self._x_opt.reshape(1,self._x_opt.size)[0]

        if  numpy.isnan(self._x_opt).any():
            print "<congeal.py> ISNAN (conv trans to x)", numpy.isnan(self._x_opt)

    def convert_transforms_to_vtk(self):
        """ Convert all subjects' transforms to vtk objects. """
        
        self.vtk_transform_list = list()

        for subj in self._subjects:
            self.vtk_transform_list.append(subj.convert_transform_to_vtk())

        return self.vtk_transform_list

    def convert_x_to_transform(self, input_x, input_transform):
        """ Convert the input_x (in optimizer search space
        coordinates) to a transform. Pays attention to the current
        registration mode and only modifies the corresponding
        components of the previous transform, allowing preservation of
        previous registration information."""

        # copy the transform, don't modify it directly
        transform = numpy.copy(input_transform)

        # only modify the part of the transform that x represents
        # since we may have a previous/current/initial transform set

        if self._registration_mode == "TranslateOnly":
            transform[3:6] = numpy.multiply(input_x[0:3], self._x_scaling[3:6])
        elif self._registration_mode == "RotateOnly":
            transform[0:3] = numpy.multiply(input_x[0:3], self._x_scaling[0:3])
        elif self._registration_mode == "ScaleOnly":
            transform[6:9] = numpy.multiply(input_x[0:3], self._x_scaling[6:9])
        elif self._registration_mode == "ShearOnly":
            transform[9:15] = numpy.multiply(input_x[0:6], self._x_scaling[9:15])
        elif self._registration_mode == "TranslateAndRotate":
            transform[0:6] = numpy.multiply(input_x, self._x_scaling[0:6])
        elif self._registration_mode == "TranslateRotateScale":
            transform = numpy.multiply(input_x, self._x_scaling)

        return transform

    def convert_transform_to_x(self, transform):
        """ Convert a transform (9 components) to the search space (x)
        coordinates."""

        output_x = []

        if self._registration_mode == "TranslateOnly":
            output_x = numpy.divide(transform[3:6], self._x_scaling[3:6])
        elif self._registration_mode == "RotateOnly":
            output_x = numpy.divide(transform[0:3], self._x_scaling[0:3])
        elif self._registration_mode == "ScaleOnly":
            output_x = numpy.divide(transform[6:9], self._x_scaling[6:9])
        elif self._registration_mode == "ShearOnly":
            output_x = numpy.divide(transform[9:15], self._x_scaling[9:15])
        elif self._registration_mode == "TranslateAndRotate":
            output_x = numpy.divide(transform[0:6], self._x_scaling[0:6])
        elif self._registration_mode == "TranslateRotateScale":
            output_x = numpy.divide(transform, self._x_scaling)

        return output_x

    def get_rho_parameters_multisubject(self):
        """ Get step sizes in optimizer search space, for all
        subjects, in one vector. """
        
        rhobeg_list = list()
        rhoend_list = list()

        rhobeg, rhoend = self.get_rho_parameters()

        for subj in self._subjects:
            rhobeg_list.append(rhobeg)
            rhoend_list.append(rhoend)

        # convert to array
        rhobeg = numpy.array(rhobeg_list)
        rhoend = numpy.array(rhoend_list)

        # make it into a vector rather than a 2D array
        rhobeg = rhobeg.reshape(1, rhobeg.size)[0]
        rhoend = rhoend.reshape(1, rhoend.size)[0]

        return rhobeg, rhoend

    def get_rho_parameters(self):
        """ Get step sizes in optimizer search space"""

        rhobeg = []
        rhoend = []

        if self._registration_mode == "TranslateOnly":
            rhobeg = self.rhobeg[3:6]
            rhoend = self.rhoend[3:6]
        elif self._registration_mode == "RotateOnly":
            rhobeg = self.rhobeg[0:3]
            rhoend = self.rhoend[0:3]
        elif self._registration_mode == "ScaleOnly":
            rhobeg = self.rhobeg[6:9]
            rhoend = self.rhoend[6:9]
        elif self._registration_mode == "ShearOnly":
            rhobeg = self.rhobeg[9:15]
            rhoend = self.rhoend[9:15]
        elif self._registration_mode == "TranslateAndRotate":
            rhobeg = self.rhobeg[0:6]
            rhoend = self.rhoend[0:6]
        elif self._registration_mode == "TranslateRotateScale":
            rhobeg = self.rhobeg[0:9]
            rhoend = self.rhoend[0:9]

        return rhobeg, rhoend
         
    def objective_function(self, current_x):
        """ The actual objective used in registration.  Function of
        the current x in search space, as well as parameters of the
        class: threshold, sigma. Compares sampled fibers from moving
        input, to all fibers of fixed input."""

        # keep track of current x value
        self._x_opt = current_x

        # get and apply transforms from current_x
        # x is all transforms concatenated for search space
        # this tests if the transform was modified for each subject
        self.convert_xs_to_transforms()

        for subj in self._subjects:
            # apply transform to moving fiber data
            subj.apply_transform()

        subj_pairs = list()
        idx1 = 0
        for subj1 in self._subjects:
            idx2 = 0
            for subj2 in self._subjects:
                subj_pairs.append([idx1, idx2])
                idx2 += 1
            idx1 += 1

        # The Entropy objective was tested and published
        if self._objective_function_mode == "Entropy":
            ret = Parallel(
                n_jobs=self.parallel_jobs, verbose=self.parallel_verbose)(
                    delayed(inner_loop_objective)(self._subjects[idx],
                                              idx, self._subjects,
                                              self.threshold, self._sigmasq,
                                              self.pairwise_subject_similarity, self.distance_method)
                                              for idx in range(self.number_of_subjects))

            for idx in range(self.number_of_subjects):
                self.pairwise_subject_similarity[idx,:] = ret[idx][0]
                self.computed_pairwise_similarity[idx] = ret[idx][1][0]
        else:
            print "<congeal.py> Chosen objective function NOT IMPLEMENTED"
		
        if self.verbose:
            print "<congeal.py> Pairwise similarity:"
            print numpy.round(self.pairwise_subject_similarity * 100) / 100
            print "<congeal.py> Computed similarity:"
            print self.computed_pairwise_similarity.astype(int)

        # The Entropy objective was tested and published
        if self._objective_function_mode == "Entropy":
            # sum negative log probabilities of all fibers, across all brains.
            # sum for total probability of each brain, in leave-one-out model
            #total_similarity = numpy.sum(self.pairwise_subject_similarity, axis=0)
            #obj = numpy.sum(total_similarity)
            obj = numpy.sum(self.pairwise_subject_similarity)
        else:
            print "<congeal.py> Chosen objective function NOT IMPLEMENTED"

        # save objective function value for analysis of performance
        self.objective_function_values.append(obj)

        if self.verbose:
            #print "<congeal.py> S:",  total_similarity
            print "<congeal.py> O:",  obj
            print "<congeal.py>             X:", self._x_opt
            
        return obj

    def c0(self, current_x):
        return self.constraint_component(current_x, 0)

    def c1(self, current_x):
        return self.constraint_component(current_x, 1)    

    def c2(self, current_x):
        return self.constraint_component(current_x, 2)    

    def c3(self, current_x):
        return self.constraint_component(current_x, 3)    

    def c4(self, current_x):
        return self.constraint_component(current_x, 4)    

    def c5(self, current_x):
        return self.constraint_component(current_x, 5)    

    def c6(self, current_x):
        return self.constraint_component(current_x, 6)    

    def c7(self, current_x):
        return self.constraint_component(current_x, 7)    

    def c8(self, current_x):
        return self.constraint_component(current_x, 8)    

    def c9(self, current_x):
        return self.constraint_component(current_x, 9)
    def c10(self, current_x):
        return self.constraint_component(current_x, 10)
    def c11(self, current_x):
        return self.constraint_component(current_x, 11)
    def c12(self, current_x):
        return self.constraint_component(current_x, 12)
    def c13(self, current_x):
        return self.constraint_component(current_x, 13)
    def c14(self, current_x):
        return self.constraint_component(current_x, 14)
    
    def constraint_component(self, current_x, component):
        """ Return appropriate constraint to keep optimizer searching
        in good region where there is no mean transform added
        meaninglessly to all subjects. """
	
        if (component < 6) | (component > 8):
            # 0-mean: translation, rotation, shear
            r0 = 0
            for subj in self._subjects:
                r0 = r0 + subj.transform[component]
            # positive return value is "good", negative is bad
            # scaling is to penalize such that actual means are close
            # enough to 0 for the param (angle, mm, scale)
            if component< 3:
                # angle
                return -10 * numpy.abs(r0)
            else:
                # mm
                return -numpy.abs(r0)
        else:
            # mean of 1.0 (scale)
            r0 = 1
            for subj in self._subjects:
                r0 = r0 * subj.transform[component]
            # positive return value is "good", negative is bad
            return -100 * numpy.abs(1.0 - r0)

    def compute(self):

        """ Run the registration.  Add subjects first (before calling
        compute). Then call compute several times, using different
        parameters for the class, for example first just for
        translation."""

        # remove any mean value (from initialization?) from transforms.
        self.remove_mean_from_transforms()
        
        # initialize optimization  using current transform
        self.convert_transforms_to_xs()

        # initialize similarity matrix
        self.pairwise_subject_similarity = numpy.zeros((self.number_of_subjects, self.number_of_subjects))
        self.computed_pairwise_similarity = numpy.zeros((self.number_of_subjects, self.number_of_subjects))
        # set up randomly sampled fibers according to current settings
        self.initialize_fiber_samples()

        # get parameters corresponding to the type of transform we are optimizing
        rhobeg, rhoend = self.get_rho_parameters_multisubject()

        # square current sigma parameter for later Gaussian
        self._sigmasq = self.sigma * self.sigma

        if self.verbose:
            print "<congeal.py> Initial value for X:", self._x_opt


        print "<congeal.py> Computing. Registration mode:", self._registration_mode

        # These optimizers were tested, less successfully than cobyla.
        # This code is left here as documentation.
        #scipy.optimize.fmin(self.objective_function,
        #                    self._x_opt, params, maxiter=self.maxiter,
        #                    ftol=self.ftol, maxfun=self.maxfun)
        #self.x_final = scipy.optimize.fmin_powell(self.objective_function,
        #                                          self._x_opt,
        #                                          xtol=rhobeg[0],
        #                                          ftol = 0.05,
        #                                          maxfun=self.maxfun)

        # Optimize using cobyla. Allows definition of initial and
        # final step size scales (rhos), as well as constraints.  Here
        # we use the constraints to encourage that meaningless
        # transforms are not tested (all brains rotate or shrink
        # together, etc.)
        self.x_final = scipy.optimize.fmin_cobyla(self.objective_function,
                                                  self._x_opt,
                                                  [self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6, self.c7, self.c8],
                                                  maxfun=self.maxfun, rhobeg=rhobeg,
                                                  rhoend=rhoend
                                                  )

        # remove any mean value from transforms.
        # Prevents all brains translating, rotating, shifting together.
        self.remove_mean_from_transforms()

        # Return output transforms from this iteration
        return self.convert_transforms_to_vtk()


def inner_loop_objective(subj1, idx1, subjects, threshold, sigmasq, global_subject_similarity, distance_method):
    """ The code called within the objective_function to find the
    negative log probability of one brain given all other brains. Used
    with Entropy objective function."""

    pairwise_similarity = numpy.zeros((len(subjects), 1))
    computed_similarity = numpy.zeros((len(subjects), 1))
    
    idx2 = 0

    # data structure to hold probability of every fiber in subject 1
    # initialize to a tiny probability to avoid the case where a 0
    # probability fiber kills all information in the objective
    # function.
    probability = numpy.zeros(len(subj1._moving_fiber_sample)) + 1e-20

    # number of compared fibers (normalization factor)
    number_comparisons = 0
    
    # loop over all subjects to find probability of each fiber in subj1
    for subj2 in subjects:
        # leave self out of model
        if subj1 != subj2:
            ## compute if changed, or first iteration
            ##if (subj1.modified | subj2.modified)  | (global_subject_similarity[idx1, idx2] == 0):
            ## caching not implemented yet always compute
            
            # Loop over fibers in brain 1, find total probability of
            # fiber using all other fibers from other brains.
            for idx in range(len(subj1._moving_fiber_sample)):
                probability[idx] += whitematteranalysis.similarity.total_similarity(
                    subj1._moving_fibers.get_fiber(subj1._moving_fiber_sample[idx]),
                    subj2._moving_fibers,
                    threshold,
                    sigmasq, distance_method)
                    
                   
            # Record that the similarity was computed for this brain
            # pair.
            computed_similarity[idx2] = 1

            # Record how many fiber comparisons were made so far
            number_comparisons += subj2._moving_fibers.number_of_fibers
            
            # end if subjects are not the same

        idx2 = idx2 + 1

    # divide total probability by number of fiber comparisons
    # the output must be between 0 and 1
    probability /= number_comparisons

    # add negative log probabilities of all fibers in this brain.
    pairwise_similarity = numpy.sum(- numpy.log(probability))
    
    return pairwise_similarity, computed_similarity

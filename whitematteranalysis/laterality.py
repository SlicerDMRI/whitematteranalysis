# -*- coding: utf-8 -*-

""" laterality.py

This module calculates white matter fiber laterality indices.

functions

def laterality_index(left, right, idx=None)

Just calculates hemisphere laterality index from array inputs left and right.
LI= (R-L)/(R+L). If R==L==0, output is 0.

class ComputeWhiteMatterLaterality

Class that actually does the computation. Input is polydata,
parameters include sigma, points used to parametrize the fiber, and
lower threshold on fiber distance.  Output is of class
LateralityResults (io.py)

"""

import os
import warnings

import numpy as np
import vtk

from whitematteranalysis.utils.opt_pckg import optional_package

from . import filter, similarity
from .fibers import FiberArray
from .io import LateralityResults

joblib, have_joblib, _ = optional_package("joblib")
Parallel, _, _ = optional_package("joblib.Parallel")
delayed, _, _ = optional_package("joblib.delayed")

if not have_joblib:
    warnings.warn(joblib._msg)
    warnings.warn("Cannot multiprocess.")


def compute_laterality_index(left, right, idx=None):
    ''' Compute laterality index from left and right hemisphere quantities.'''

    laterality_index = np.zeros(len(left))

    if idx == None:
        # if L=R=0, output 0. (avoid divide by 0, skip masked out data)
        idx = np.nonzero(left)[0] & np.nonzero(right)[0]

    # otherwise use input idx to select relevant data, to avoid doing
    # the above nonzero every time in a loop

    laterality_index[idx] = (right[idx] - left[idx]) / (right[idx] + left[idx])

    return laterality_index


class WhiteMatterLaterality:

    """Laterality computation from fiber tracts."""

    def __init__(self):
        # parameters
        # medium sigma
        self.sigma = 20.0
        # stricter and more accurate than 5 points
        self.points_per_fiber = 20
        # aid in matching across hemispheres with threshold
        self.threshold = 5.0
        # same number of fibers measured from each hemisphere
        # avoid biasing due just to number of fibers
        self.use_equal_fibers = True
        # If using equal fibers, can set how many to use same number across subjects
        self.fibers_per_hemisphere = None
        
        # performance options
        self.verbose = True
        # set parallel_jobs to 0 to turn off multiprocessing
        self.parallel_jobs = 2
        self.parallel_verbose = 0

        # internal data storage
        self.fibers = FiberArray()

    def __str__(self):
        output = f" sigma\t\t\t{str(self.sigma)}\n points_per_fiber\t{str(self.points_per_fiber)}\n threshold\t\t{str(self.threshold)}\n verbose\t\t{str(self.verbose)} \n parallel_jobs\t\t{str(self.parallel_jobs)}\n parallel_verbose\t{str(self.parallel_verbose)}\n fibers\n\t\t\t{str(self.fibers)}"

        return output

    def compute(self, input_vtk_polydata):
        """ Actually calculate the laterality index for every input
        fiber.

        Input polydata is required. This polydata is modified by
        adding a cell data array containing laterality indices.

        Output from this class is a struct: <io.py> class
        LateralityResults

        Parameters in the class can also be modified for experiments:
        sigma (in Gaussian on inter-fiber-point distance),
        points_per_fiber (for parameterization), threshold (below
        which inter-fiber-point distance is set to 0).

        Performance options are the number of parallel_jobs, and
        verbose (whether to print progress).

        """

        # internal representation for fast similarity computation
        # this also detects which hemisphere fibers are in
        self.fibers.points_per_fiber = self.points_per_fiber
        # must request hemisphere computation from object
        self.fibers.hemispheres = True
        # Now convert to array with points and hemispheres as above
        self.fibers.convert_from_polydata(input_vtk_polydata)
        
        # get the same number from each hemisphere if requested
        # -------------------------
        if self.use_equal_fibers:
            num_fibers = min(self.fibers.number_left_hem, self.fibers.number_right_hem)        
            if self.fibers_per_hemisphere is not None:
                if self.fibers_per_hemisphere <= num_fibers:
                    num_fibers = self.fibers_per_hemisphere
                else:
                    raise Exception(f"Fibers per hemisphere is set too high for the dataset. Current subject maximum is {str(num_fibers)}")
        
            # grab num_fibers fibers from each hemisphere.
            # use the first n since they were randomly sampled from the whole dataset
            selected_right = self.fibers.index_right_hem[0:num_fibers]
            selected_left = self.fibers.index_left_hem[0:num_fibers]
            mask = np.zeros(input_vtk_polydata.GetNumberOfLines())
            mask[selected_right] = 1
            mask[selected_left] = 1
            # go back to the input data and use just those fibers
            input_vtk_polydata = filter.mask(input_vtk_polydata, mask)
            # Now convert to array with points and hemispheres as above
            self.fibers.convert_from_polydata(input_vtk_polydata)
            if self.verbose:
                print(f"<{os.path.basename(__file__)}> Using {num_fibers} fibers per hemisphere.")
                
        # square sigma for later Gaussian
        sigmasq = self.sigma * self.sigma

        # allocate outputs
        nf = self.fibers.number_of_fibers
        laterality_index = np.zeros(nf)
        right_hem_total = np.zeros(nf)
        left_hem_total = np.zeros(nf)
        #right_hem_distance = np.zeros([nf, nf])
        #left_hem_distance = np.zeros([nf, nf])


        # grab all fibers from each hemisphere
        fiber_array_right = self.fibers.get_fibers(self.fibers.index_right_hem)
        fiber_array_left = self.fibers.get_fibers(self.fibers.index_left_hem)

        # tell user we are doing something
        if self.verbose:
            print(f"<{os.path.basename(__file__)}> Fibers in each hemisphere. L: {self.fibers.number_left_hem} R: {self.fibers.number_right_hem} / Total: {self.fibers.number_of_fibers}")
            print(f"<{os.path.basename(__file__)}> Starting to compute laterality indices")

        # run the computation, either in parallel or not
        if (have_joblib & (self.parallel_jobs > 1)):
            if self.verbose:
                print(f"<{os.path.basename(__file__)}> Starting parallel code. Processes: {self.parallel_jobs}")

            # compare to right hemisphere (reflect fiber first if in left hem)
            ret = Parallel(
                n_jobs=self.parallel_jobs, verbose=self.parallel_verbose)(
                delayed(similarity.total_similarity_for_laterality)(
                    self.fibers.get_fiber(lidx),
                    fiber_array_right,
                    self.fibers.is_left_hem[lidx],
                    self.threshold,
                    sigmasq)
                for lidx in self.fibers.index_hem)

            #ret = zip(*ret)
            right_hem_total[self.fibers.index_hem] = ret
            #right_hem_distance = ret[1]

            # compare to left hemisphere (reflect fiber first if in right hem)
            ret = Parallel(
                n_jobs=self.parallel_jobs, verbose=self.parallel_verbose)(
                delayed(similarity.total_similarity_for_laterality)(
                    self.fibers.get_fiber(lidx),
                    fiber_array_left,
                    self.fibers.is_right_hem[lidx],
                    self.threshold,
                    sigmasq)
                for lidx in self.fibers.index_hem)
            #ret = zip(*ret)
            left_hem_total[self.fibers.index_hem] = ret
            #left_hem_distance = ret[1]

        else:
            right_hem_distance = np.zeros([nf, len(self.fibers.index_right_hem)])
            left_hem_distance = np.zeros([nf, len(self.fibers.index_left_hem)])

            # compare to right hemisphere (reflect fiber first if in left hem)
            for lidx in self.fibers.index_hem:
                ret = similarity.total_similarity_for_laterality(
                    self.fibers.get_fiber(lidx),
                    fiber_array_right,
                    self.fibers.is_left_hem[lidx],
                    self.threshold,
                    sigmasq)
                right_hem_total[lidx] = ret
                #right_hem_total[lidx] = ret[0]
                #right_hem_distance[lidx,:] = ret[1]

            # compare to left hemisphere (reflect fiber first if in right hem)
            for lidx in self.fibers.index_hem:
                ret = similarity.total_similarity_for_laterality(
                    self.fibers.get_fiber(lidx),
                    fiber_array_left,
                    self.fibers.is_right_hem[lidx],
                    self.threshold,
                    sigmasq)
                left_hem_total[lidx] = ret
                #left_hem_distance[lidx,:] = ret[1]

        laterality_index = compute_laterality_index(left_hem_total,
                                                    right_hem_total,
                                                    self.fibers.index_hem)


        # output the LI as cell data in the polydata
        # for visualization and/or further analyses
        cell_data = vtk.vtkFloatArray()
        cell_data.SetName('Laterality')
        for lidx in range(0, self.fibers.number_of_fibers):
            cell_data.InsertNextTuple1(laterality_index[lidx])
            input_vtk_polydata.GetCellData().SetScalars(cell_data)

        # output everything
        results = LateralityResults()
        results.laterality_index = laterality_index
        results.polydata = input_vtk_polydata
        #results.right_hem_distance = right_hem_distance
        #results.left_hem_distance = left_hem_distance
        results.sigma = self.sigma
        results.points_per_fiber = self.points_per_fiber
        results.threshold = self.threshold
        results.left_hem_similarity = left_hem_total
        results.right_hem_similarity = right_hem_total
        results.hemisphere = self.fibers.fiber_hemisphere
        return results

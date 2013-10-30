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

import numpy

import vtk
try:
    from joblib import Parallel, delayed
    USE_PARALLEL = 1
except ImportError:
    USE_PARALLEL = 0
    print "<laterality.py> Failed to import joblib, cannot multiprocess."
    print "<laterality.py> Please install joblib for this functionality."

from fibers import FiberArray
import similarity
from io import LateralityResults


def compute_laterality_index(left, right, idx=None):
    ''' Compute laterality index from left and right hemisphere quantities.'''

    laterality_index = numpy.zeros(len(left))

    if idx == None:
        # if L=R=0, output 0. (avoid divide by 0, skip masked out data)
        idx = numpy.nonzero(left)[0] & numpy.nonzero(right)[0]

    # otherwise use input idx to select relevant data, to avoid doing
    # the above nonzero every time in a loop

    laterality_index[idx] = (right[idx] - left[idx]) / (right[idx] + left[idx])

    return laterality_index


class WhiteMatterLaterality:

    """Laterality computation from fiber tracts."""

    def __init__(self):
        # parameters
        self.sigma = 10
        self.points_per_fiber = 5
        self.threshold = 5

        # performance options
        self.verbose = True
        # set parallel_jobs to 0 to turn off multiprocessing
        self.parallel_jobs = 2
        self.parallel_verbose = 0

        # internal data storage
        self.fibers = FiberArray()

    def __str__(self):
        output = " sigma\t\t\t" + str(self.sigma) \
            + "\n points_per_fiber\t" + str(self.points_per_fiber) \
            + "\n threshold\t\t" + str(self.threshold) \
            + "\n verbose\t\t" + str(self.verbose) \
            + "\n parallel_jobs\t\t" + str(self.parallel_jobs) \
            + "\n parallel_verbose\t" + str(self.parallel_verbose) \
            + "\n fibers\n\t\t\t" \
            + str(self.fibers)

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

        # square sigma for later Gaussian
        sigmasq = self.sigma * self.sigma

        # allocate outputs
        nf = self.fibers.number_of_fibers
        laterality_index = numpy.zeros(nf)
        right_hem_total = numpy.zeros(nf)
        left_hem_total = numpy.zeros(nf)
        #right_hem_distance = numpy.zeros([nf, nf])
        #left_hem_distance = numpy.zeros([nf, nf])

        # grab fibers from each hemisphere
        fiber_array_right = self.fibers.get_fibers(self.fibers.index_right_hem)
        fiber_array_left = self.fibers.get_fibers(self.fibers.index_left_hem)

        # tell user we are doing something
        if self.verbose:
            print "<laterality.py> Fibers in each hemisphere.", \
                "L:", self.fibers.number_left_hem, \
                "R:", self.fibers.number_right_hem, \
                "/ Total:", self.fibers.number_of_fibers
            print "<laterality.py> Starting to compute laterality indices"

        # run the computation, either in parallel or not
        if (USE_PARALLEL & (self.parallel_jobs > 1)):
            if self.verbose:
                print "<laterality.py> Starting parallel code. Processes:", \
                    self.parallel_jobs

            # compare to right hemisphere (reflect fiber first if in left hem)
            ret = Parallel(
                n_jobs=self.parallel_jobs, verbose=self.parallel_verbose)(
                delayed(similarity.total_similarity_and_distances)(
                    self.fibers.get_fiber(lidx),
                    fiber_array_right,
                    self.fibers.is_left_hem[lidx],
                    self.threshold,
                    sigmasq)
                for lidx in self.fibers.index_hem)

            ret = zip(*ret)
            right_hem_total[self.fibers.index_hem] = ret[0]
            right_hem_distance = ret[1]

            # compare to left hemisphere (reflect fiber first if in right hem)
            ret = Parallel(
                n_jobs=self.parallel_jobs, verbose=self.parallel_verbose)(
                delayed(similarity.total_similarity_and_distances)(
                    self.fibers.get_fiber(lidx),
                    fiber_array_left,
                    self.fibers.is_right_hem[lidx],
                    self.threshold,
                    sigmasq)
                for lidx in self.fibers.index_hem)
            ret = zip(*ret)
            left_hem_total[self.fibers.index_hem] = ret[0]
            left_hem_distance = ret[1]

        else:
            right_hem_distance = numpy.zeros([nf, len(self.fibers.index_right_hem)])
            left_hem_distance = numpy.zeros([nf, len(self.fibers.index_left_hem)])

            # compare to right hemisphere (reflect fiber first if in left hem)
            for lidx in self.fibers.index_hem:
                ret = similarity.total_similarity_and_distances(
                    self.fibers.get_fiber(lidx),
                    fiber_array_right,
                    self.fibers.is_left_hem[lidx],
                    self.threshold,
                    sigmasq)
                right_hem_total[lidx] = ret[0]
                right_hem_distance[lidx,:] = ret[1]

            # compare to left hemisphere (reflect fiber first if in right hem)
            for lidx in self.fibers.index_hem:
                ret = similarity.total_similarity_and_distances(
                    self.fibers.get_fiber(lidx),
                    fiber_array_left,
                    self.fibers.is_right_hem[lidx],
                    self.threshold,
                    sigmasq)
                left_hem_total[lidx] = ret[0]
                left_hem_distance[lidx,:] = ret[1]

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

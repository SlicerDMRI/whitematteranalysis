""" relative_distance.py

implementation of relative distance model, to test representation of
fiber tracts represented in coordinates relative to fMRI and
anatomical points

class RelativeDistanceModel


"""


import numpy

import vtk
try:
    from joblib import Parallel, delayed
    USE_PARALLEL = 1
except ImportError:
    USE_PARALLEL = 0
    print("<relative_distance.py> Failed to import joblib, cannot multiprocess.")
    print("<relative_distance.py> Please install joblib for this functionality.")


import whitematteranalysis.fibers
import whitematteranalysis.similarity


class RelativeDistanceModel:

    
    def __init__(self):
        # parameters that can be set by user
        self.points_per_fiber = 5

        # performance options set by user
        self.parallel_verbose = 0
        self.parallel_jobs = 4

        # data stored in model
        self.distances = []
        self.points = []

        # internal data storage, not part of model
        self._fibers = whitematteranalysis.fibers.FiberArray()

    def compute(self, input_polydata, input_points):

        # internal representation for fast similarity computation
        self._fibers.convert_from_polydata(input_polydata,
                                           self.points_per_fiber)

        # loop over all points.
        # 3D array output: distances[ fiber_number, point_number, arc_length]

        self.points = input_points
        self.distances = self._compute_relative_distance(self._fibers, input_points)

    def _compute_relative_distance(self, fibers, point):
        """
        Find distance from fibers to a point.
        
        input fiber should be class Fiber. point should be class ????
        
        """
        
        nfib = fibers.number_of_fibers 
        
        # like repmat. copy point to full matrix size of all lines
        #point_r = numpy.tile(point.r, (nfib, 1))
        #point_a = numpy.tile(point.a, (nfib, 1))
        #point_s = numpy.tile(point.s, (nfib, 1))
        
        # compute the distance from each point to the array of fibers
        dx = fibers.fiber_array_r - point[0]
        dy = fibers.fiber_array_a - point[1]
        dz = fibers.fiber_array_s - point[2]

        dx = numpy.power(dx, 2)
        dy = numpy.power(dy, 2)
        dz = numpy.power(dz, 2)
        
        # sum dx dx dz at each point on the fiber and sqrt for threshold
        distance = numpy.sqrt(dx + dy + dz)
    
        return distance



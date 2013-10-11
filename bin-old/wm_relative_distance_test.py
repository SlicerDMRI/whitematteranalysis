#!/usr/bin/env python

import whitematteranalysis as wma
print 'Read and preprocess'


number_of_fibers = 2000


# read list of points
# why can't these be polydata also?
import numpy
point = numpy.zeros(3)
point[0] = 1
point[1] = 2
point[2] = 3

pd = wma.io.read_polydata('/Users/lauren/Data/TBI/Tracts/01030-dwi-filt-Ed-DTI-tract.vtp')
pd2 = wma.filter.preprocess(pd,80)
pd2 = wma.filter.downsample(pd2,number_of_fibers)


# compute model for each subject
rd = wma.relative_distance.RelativeDistanceModel()
rd.compute(pd2, point)

# compute model across subjects


# predict in leave one out, measure overlap


# detect in leave one out, measure overlap


# compute model for each pt


# predict in each pt


# detect in each pt



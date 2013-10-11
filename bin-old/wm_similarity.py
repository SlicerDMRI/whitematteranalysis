#!/usr/bin/env python

import whitematteranalysis as wma
import numpy

pd = wma.io.read_polydata('/Users/lauren/Data/TBI/Tracts/01035-dwi-filt-Ed-DTI-tract.vtp')
pd2 = wma.filter.preprocess(pd,80)
pd2 = wma.filter.downsample(pd2,1000)

# should I downsample pd rather than subset of fiber array?????

fiber_array = wma.fibers.FiberArray()
fiber_array.convert_from_polydata(pd2)
#fiber = fiber_array.get_fiber(100)
fidx=111
fiber = fiber_array.get_fiber(fidx)
threshold = 5
distance = wma.similarity.fiber_distance(fiber, fiber_array, threshold)
similarity = wma.similarity.distance_to_similarity(distance, 100)

fiber_mask = numpy.ones(len(distance))
fiber_mask = similarity > 0.2

pd3 = wma.filter.mask(pd2, fiber_mask, similarity)

#wma.render.render(pd3,[0,  numpy.max(distance)])

wma.render.render(pd3)

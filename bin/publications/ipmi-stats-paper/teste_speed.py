import whitematteranalysis as wma
import numpy
import vtk

pd = wma.io.read_polydata('/Users/odonnell/Dropbox/Coding/OUTPUTS/Dec2012/ipmi-reg-gender/iteration_6/white_matter_0000.vtk')
print 'Converting fibers to array representation for dist and averaging'
fiber_array = wma.fibers.FiberArray()
fiber_array.verbose = 1
points_per_fiber = 30
fiber_array.convert_from_polydata(pd, points_per_fiber)


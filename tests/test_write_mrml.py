from __future__ import print_function
import numpy
import whitematteranalysis as wma

pd_filenames = list()
pd_filenames.append('test_data_small/brain_0001.vtk')
pd_filenames.append('test_data_small/brain_0002.vtk')
colors = list()
colors.append([255, 0, 100])
colors.append([0, 200, 150])
colors = numpy.array(colors)
wma.mrml.write(pd_filenames, colors , 'try.mrml')

print('Opening try.mrml in slicer should display a red and a green fiber bundle.')
print('try.mrml must be in the test directory for this to work. It references data in test_data_small.')

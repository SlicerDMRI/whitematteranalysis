import os
import glob
import matplotlib.pyplot as plt
import numpy
import scipy.stats

import vtk

import whitematteranalysis as wma

import multiprocessing


indir = '/Users/odonnell/Desktop/OLD-Results-MICCAI-notused/tbi_with_scalars'

input_mask = "{0}/*.vtk".format(indir)
input_poly_datas = glob.glob(input_mask)

print input_poly_datas


input_pds = list()
input_pds_downsampled = list()

number_of_fibers_per_subject = 2000
number_of_subjects = len(input_poly_datas)
points_per_fiber = 10


# read in ones with scalars already
# this is SLOW
for fname in input_poly_datas:
    print fname
    pd = wma.io.read_polydata(fname)
    input_pds.append(pd)

# downsample for analysis
input_mean_fas_per_subject = list()
input_pds_downsampled = list()
downsample_indices = list()
for pd in input_pds:
    pd2, fiber_indices = wma.filter.downsample(pd, number_of_fibers_per_subject,return_indices=True)
    input_pds_downsampled.append(pd2)
    downsample_indices.append(fiber_indices)

input_pds_as_arrays = list()
for pd2 in input_pds_downsampled:
    # convert to array representation
    print 'Converting fibers to array representation for dist and averaging'
    fiber_array = wma.fibers.FiberArray()
    fiber_array.convert_from_polydata(pd2, points_per_fiber)
    input_pds_as_arrays.append(fiber_array)

print 'Done converting fibers to array representation for dist and averaging'




group_indices =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,\
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

group_indices = numpy.array(group_indices)

g0 = numpy.nonzero(group_indices==0)
g1 = numpy.nonzero(group_indices)

# now look at laterality information
diff_lines = list()
li_lines = list()
for pda in input_pds_as_arrays:
    print numpy.divide(numpy.sum(pda.fiber_hemisphere), numpy.sum(numpy.abs(pda.fiber_hemisphere)))
    diff_lines.append(numpy.sum(pda.fiber_hemisphere))
    li_lines.append(numpy.divide(numpy.sum(pda.fiber_hemisphere), numpy.sum(numpy.abs(pda.fiber_hemisphere))))
diff_lines = numpy.array(diff_lines)
li_lines = numpy.array(li_lines)
abs_li_lines = numpy.abs(li_lines)
abs_diff_lines = numpy.abs(diff_lines)

scipy.stats.ttest_ind(diff_lines[g0],diff_lines[g1])
scipy.stats.ttest_ind(li_lines[g0],li_lines[g1])
scipy.stats.ttest_ind(abs_li_lines[g0],abs_li_lines[g1])
scipy.stats.ttest_ind(abs_diff_lines[g0],abs_diff_lines[g1])

plt.boxplot([abs_li_lines[g0],abs_li_lines[g1]])
plt.boxplot([abs_diff_lines[g0],abs_diff_lines[g1]])

num_lines = list()
for pd in input_pds:
    print pd.GetNumberOfLines()
    num_lines.append(pd.GetNumberOfLines())


num_lines = numpy.array(num_lines)
numpy.mean(num_lines[g0])
numpy.mean(num_lines[g1])

plt.boxplot([num_lines[g0],num_lines[g1]])
scipy.stats.ttest_ind(num_lines[g0],num_lines[g1])

age_control = [29, 43, 38, 31, 23, 29, 40, 24, 42, 26, 47, 23, 40]
age_tbi = [44, 37, 43, 27, 29, 42, 27, 24, 25, 29, 24, 39, 44]
age_control= numpy.array(age_control)
age_tbi= numpy.array(age_tbi)
age =  numpy.array(list(age_control) + list(age_tbi))

slope, intercept, r, p, err = scipy.stats.linregress(num_lines[g0], age[g0])

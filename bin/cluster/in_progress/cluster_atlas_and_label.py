#!/usr/bin/env python

import whitematteranalysis as wma
import vtk
import numpy
import matplotlib.pyplot as plt
import multiprocessing
import glob

number_of_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', number_of_jobs

print 'Read and preprocess'

indir = '/Users/lauren/Dropbox/Coding/OUTPUTS/Dec2012/ipmi-reg-gender/iteration_6'
outdir = '/Users/lauren/Dropbox/Coding/OUTPUTS/March2013'

# parameters for clustering/atlas creation
number_of_clusters = 400
number_of_fibers_per_subject = 1000
minimum_length = 40
number_of_sampled_fibers = 500
number_of_eigenvectors = 15
sigma = 50

# parameters for atlas labeling
number_of_fibers_to_label = 4000

# read all polydata in the directory
input_mask = "{0}/*.vtk".format(indir)
input_poly_datas = glob.glob(input_mask)
#input_poly_datas = input_poly_datas[0:1]
input_poly_datas = input_poly_datas[0:6]
print input_poly_datas

# read in data
input_pds = list()
for fname in input_poly_datas:
    print fname
    pd = wma.io.read_polydata(fname)
    pd2 = wma.filter.preprocess(pd, minimum_length)
    pd3 = wma.filter.downsample(pd2, number_of_fibers_per_subject)
    input_pds.append(pd3)
    del pd
    del pd2
    del pd3

# append input data
appender = vtk.vtkAppendPolyData()
for pd in input_pds:
    appender.AddInput(pd)

appender.Update()
input_data = appender.GetOutput()

del input_pds

# CLUSTERING step
nystrom_mask = numpy.random.permutation(input_data.GetNumberOfLines()) < number_of_sampled_fibers

output_polydata_s, cluster_numbers_s, color, embed, distortion, atlas = \
    wma.cluster.spectral(input_data, number_of_clusters=number_of_clusters, \
                             number_of_jobs=number_of_jobs, use_nystrom=True, \
                             nystrom_mask = nystrom_mask, \
                             number_of_eigenvectors = number_of_eigenvectors, \
                             sigma = sigma)

print 'View results'
plt.figure()
plt.hist(cluster_numbers_s, 200)
plt.savefig('cluster_hist.pdf')
plt.close()

# view the whole thing
ren1 = wma.render.render(output_polydata_s, 1000)

# view it cluster by cluster (run line 2 several times)
cidx = 1
#fiber_mask = cluster_numbers_s == cidx; pd4 = wma.filter.mask(output_polydata_s, fiber_mask, color); wma.render.render(pd4); cidx = cidx+1

ren2 = wma.cluster.view_cluster_number(output_polydata_s, cidx, cluster_numbers_s); cidx = cidx+1

#wma.cluster.view_cluster_number(output_polydata_s, cidx)
# cidx = cidx+1; wma.cluster.view_cluster_number(output_polydata_s, cidx)
# evaluate

# to view clusters
if 0:
    output_polydata_s.GetCellData().SetActiveScalars('ClusterNumber')
    ren2 = wma.render.render(output_polydata_s, 1000)

# compute distances to centroids
diff = atlas.centroids[cluster_numbers_s] - embed
centroid_distance = numpy.sum(numpy.multiply(diff,diff), 1)
fiber_mask = centroid_distance > 10.0
pd_dist = wma.filter.mask(output_polydata_s, fiber_mask, centroid_distance)
ren_dmax = wma.render.render(pd_dist, 1000)
plt.hist(centroid_distance,1000); plt.savefig('centroid_distances.pdf')

fiber_mask = centroid_distance < 1.0
pd_dist = wma.filter.mask(output_polydata_s, fiber_mask, centroid_distance)
ren_dmin = wma.render.render(pd_dist, 1000)



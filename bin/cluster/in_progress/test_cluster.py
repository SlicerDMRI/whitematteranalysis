#!/usr/bin/env python

import whitematteranalysis as wma
import vtk
import numpy
import matplotlib.pyplot as plt
print 'Read and preprocess'

number_of_clusters = 100
number_of_fibers = 2000
number_of_fibers = 1000
import multiprocessing
number_of_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', number_of_jobs
#number_of_jobs = 12
#pd = wma.io.read_polydata('/Users/lauren/Dropbox/Coding/OUTPUTS/Dec2012/ipmi-gender5/downsampled/white_matter_0000_downsampled.vtk')
pd = wma.io.read_polydata('/Users/lauren/Dropbox/Coding/OUTPUTS/Miccai2013/test_tbi_1/atlas_info.vtp')
#pd = wma.io.read_polydata('/Users/odonnell/Dropbox/Data/TBI_FE_PNL/controls/01231-dwi-filt-Ed-DTI-tract.vtp')
#pd = wma.io.read_polydata('/Users/odonnell/Dropbox/Data/fMRIDTI/controls/fibers/subject_1_fibers.vtk')

#pd2 = wma.filter.preprocess(pd,80)
#pd3 = wma.filter.downsample(pd2,number_of_fibers)

pd3 = pd
#wma.render.render(pd3)
#wma.render.render(pd2)

if 0:
    # run the clustering
    #output_polydata_h, cluster_idx_h = wma.cluster.hierarchical(pd3)
    output_polydata_h, cluster_idx_h = wma.cluster.hierarchical(pd3,300,2,0.95,number_of_jobs)
    wma.render.render(output_polydata_h)
    import matplotlib.pyplot as plt
    print 'View results'
    plt.figure()
    not_used = plt.hist(cluster_idx_h, max(cluster_idx_h))
    ren1 = wma.render.render(output_polydata_h)
    #wma.cluster.view_cluster_number(output_polydata_h,1)

# make another polydata to put the spectral information into
pd4 = vtk.vtkPolyData()
pd4.DeepCopy(pd3)

#output_polydata_s, cluster_numbers_s, color, embed, distortion = wma.cluster.spectral(pd4,number_of_jobs=number_of_jobs)

number_of_sampled_fibers = 400
#nystrom_mask = numpy.ones(pd4.GetNumberOfLines())
nystrom_mask = numpy.random.permutation(pd4.GetNumberOfLines()) < number_of_sampled_fibers

output_polydata_s, cluster_numbers_s, color, embed, distortion, atlas = \
    wma.cluster.spectral(pd4,number_of_clusters=number_of_clusters, number_of_jobs=number_of_jobs, use_nystrom=True, nystrom_mask = nystrom_mask)

print 'View results'
plt.figure()
plt.hist(cluster_numbers_s, 200)
plt.savefig('cluster_hist.pdf')
plt.close()

# view the whole thing
ren2 = wma.render.render(output_polydata_s, 1000)

# view it cluster by cluster (run line 2 several times)
cidx = 1
#fiber_mask = cluster_numbers_s == cidx; pd4 = wma.filter.mask(output_polydata_s, fiber_mask, color); wma.render.render(pd4); cidx = cidx+1

#wma.cluster.view_cluster_number(output_polydata_s,1,cluster_numbers_s)

#wma.cluster.view_cluster_number(output_polydata_s, cidx)
# cidx = cidx+1; wma.cluster.view_cluster_number(output_polydata_s, cidx)
# evaluate

# to view clusters
if 0:
    output_polydata_s.GetCellData().SetActiveScalars('ClusterNumber')
    ren2 = wma.render.render(output_polydata_s, 1000)

# now try using the above as the atlas and re-labeling to test. Should give same answers.
pd5 = vtk.vtkPolyData()
pd5.DeepCopy(pd3)
output_polydata_2, cluster_numbers_2, color2, embed2 = \
    wma.cluster.spectral_atlas_label(pd5, atlas)
ren3 = wma.render.render(output_polydata_2, 1000)

# test if the cluster numbers match
print "Test cluster indices are the same in both clustering for atlas creation and labeling of the same data."
print "This value should be 0:", numpy.max(cluster_numbers_s - cluster_numbers_2)

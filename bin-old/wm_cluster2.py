#!/usr/bin/env python

import whitematteranalysis as wma

print 'Read and preprocess'

pd = wma.io.read_polydata('/Users/lauren/Data/TBI/Tracts/01035-dwi-filt-Ed-DTI-tract.vtp')
pd2 = wma.filter.preprocess(pd,60)
#pd2 = wma.filter.downsample(pd2,6000)
#pd2 = wma.filter.downsample(pd2,4000)
# note the below is expensive and can be done as part of clustering also
#pd2 = wma.filter.remove_outliers(pd2, 2)
#pd2 = wma.filter.downsample(pd2,1000)
pd3 = wma.filter.downsample(pd2,2000)
#pd3 = wma.filter.downsample(pd2,5000)

output_polydata, cluster_idx = wma.cluster.hierarchical(pd3)
output_polydata, cluster_idx = wma.cluster.hierarchical(pd3,300,2,0.95)
wma.render.render(output_polydata)
import matplotlib.pyplot as plt
print 'View results'
plt.figure()
not_used = plt.hist(cluster_idx, max(cluster_idx))

wma.cluster.view_cluster_number(output_polydata,1)



output_polydata, cluster_numbers, color, embed = wma.cluster.spectral(pd3)


# view it cluster by cluster (run line 2 several times)
cidx = 1
fiber_mask = cluster_numbers == cidx; pd4 = wma.filter.mask(output_polydata, fiber_mask, color); wma.render.render(pd4); cidx = cidx+1

wma.cluster.view_cluster_number(output_polydata,1,cluster_numbers)

wma.cluster.view_cluster_number(output_polydata,1)
# evaluate

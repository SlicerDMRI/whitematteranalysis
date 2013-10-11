#!/usr/bin/env python

import whitematteranalysis as wma
import numpy

import scipy.cluster.vq as cluster
import matplotlib.pyplot as plt

try:
    from joblib import Parallel, delayed
    USE_PARALLEL = 1
except ImportError:
    USE_PARALLEL = 0
    print "<laterality.py> Failed to import joblib, cannot multiprocess."
    print "<io.py> Please install joblib for this functionality."


print 'Read and preprocess'

pd = wma.io.read_polydata('/Users/lauren/Data/TBI/Tracts/01035-dwi-filt-Ed-DTI-tract.vtp')
pd2 = wma.filter.preprocess(pd,40)
pd2 = wma.filter.downsample(pd2,6000)
#pd2 = wma.filter.downsample(pd2,4000)
pd2 = wma.filter.remove_outliers(pd2, 2)
#pd2 = wma.filter.downsample(pd2,1000)
#pd3 = wma.filter.downsample(pd2,3000)
pd3 = wma.filter.downsample(pd2,5000)

print 'Convert to array'
fiber_array = wma.fibers.FiberArray()
fiber_array.convert_from_polydata(pd3)

#cluster just one brain (or many...)

# group register (as neeeded)

# for each brain (subsample as needed, keep track of it)

# pairwise distance matrix        
# compare to right hemisphere (reflect fiber first if in left hem)
all_fibers = range(0, fiber_array.number_of_fibers)
threshold = 5

print 'Compute fiber distances'
distances = Parallel(n_jobs=3, verbose=1)(
    delayed(wma.similarity.fiber_distance)(
        fiber_array.get_fiber(lidx),
        fiber_array,
        threshold)
    for lidx in all_fibers)

distances = numpy.array(distances)

# remove outliers if desired

print 'Convert to similarity, and embed'
# similarity matrix
sigmasq = 30*30
#sigmasq = 10*10
similarity = wma.similarity.distance_to_similarity(distances, sigmasq)

# sanity check that on-diagonal elements are all 1
print "This should be 1.0: ", numpy.min(numpy.diag(similarity))

# See the paper:
# "Spectral Grouping Using the Nystrom Method"
# (D^-1/2 W D^-1/2) V = V Lambda
# embedding_i,j = V_i+i,j./sqrt(D_j,j)
# norm cuts transform of similarity matrix
row_sum = numpy.sum(similarity, axis=0)
dhat = numpy.divide(1, numpy.sqrt(row_sum))
affinity_matrix = numpy.multiply(similarity, numpy.outer(dhat, dhat.T))

# embed
w,v = numpy.linalg.eigh(affinity_matrix)

print "Eigenvalue range should be within -1 to 1:", w[0], w[-1]
print "Max eigenvalue should be equal to 1:", w[-1]

# keep top few eigenvectors
power= numpy.cumsum(w[::-1])/numpy.sum(w)
numeig=10
print "power: ", power[numeig]
embed = v[:, -numeig-1:-1]
embed = numpy.divide(embed.T, v[:,-1]).T
# reverse order of embedding so highest eigenvalue
# information is first
embed = embed[:,::-1]

print 'K-means clustering'
# k-means (or improved)
#maxClusters = 400
maxClusters = 200
#k_iter = 1
#embed_w = cluster.whiten(embed)
#centroids, errors = cluster.kmeans(embed_w, maxClusters, k_iter)
centroids, errors = cluster.kmeans(embed, maxClusters)
c, dist = cluster.vq(embed,centroids)

print 'View results'
plt.figure()
not_used = plt.hist(c,200)

# get colors from embedding (first 3 coordinates corresponding to max
# eigenvalues). Need to shift so minimum color is 0,0,0
# find minimum value of each coordinate 
embed_min = numpy.min(embed[:,0:3],0)
#color = embed[:,0:3] - embed_min
color = embed[:,0:3]
# normalize all colors to length 1
color_len = numpy.sqrt(numpy.sum(numpy.power(color,2),1))
color = numpy.divide(color.T, color_len).T
# now use colors from 0 to 255 for unsigned char (no LUT)
#color_min = numpy.min(color,0)
#color = color - color_min
# color components ranged from -1 to +1, now from 0 to 255
# actually from 10 to 245, in order not to get really dark or bright fibers
# 245-10 = 235, and 235/2 = 117
color = 127.5 + (color * 127.5)

import colorsys
h = numpy.zeros(color.shape[0])
s = numpy.zeros(color.shape[0])
v = numpy.zeros(color.shape[0])
for c_idx in range(0, color.shape[0]):
    h[c_idx], s[c_idx], v[c_idx] = colorsys.rgb_to_hsv(color[c_idx, 0],color[c_idx, 1],color[c_idx, 2])

# avoid dark black fibers that the cluster color is hard to see
# convert to hsv, and make all brightness values the same
# that way shadows show geometry, not cluster number
#v2 = wma.render.histeq(v, 200)[0] + 50
#v2 = histeq(v,200)[0] + 50
v2=numpy.ones(v.shape)*180
for c_idx in range(0, color.shape[0]):
    color[c_idx, 0], color[c_idx, 1], color[c_idx, 2] = colorsys.hsv_to_rgb(h[c_idx], s[c_idx], v2[c_idx])

plt.figure()
plt.plot(embed[:,0],color[:,0],'r.')
plt.plot(embed[:,1],color[:,1],'g.')
plt.plot(embed[:,2],color[:,2],'b.')

fiber_mask = numpy.ones(len(c))
pd4 = wma.filter.mask(pd3, fiber_mask, color)
# View embedding as RGB 
ren = wma.render.render(pd4)

wma.render.save_views(ren)

# view it cluster by cluster (run line 2 several times)
cidx = 1
fiber_mask = c == cidx; pd4 = wma.filter.mask(pd3, fiber_mask, color); wma.render.render(pd4); cidx = cidx+1

# evaluate

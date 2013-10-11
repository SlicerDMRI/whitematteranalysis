#!/usr/bin/env python

import whitematteranalysis as wma
import vtk
import numpy
import matplotlib.pyplot as plt
import multiprocessing
import glob
import pickle

number_of_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', number_of_jobs

print 'Read and preprocess'

#indir = '/Users/lauren/Dropbox/Coding/OUTPUTS/Dec2012/ipmi-reg-gender/iteration_6'
indir = '/Users/lauren/Dropbox/Coding/OUTPUTS/March2013/input_control_atlas/one_t_trad_binary'
outdir = '/Users/lauren/Dropbox/Coding/OUTPUTS/March2013'

landmark_file = '/Users/odonnell/Dropbox/Coding/OUTPUTS/March2013/input_control_atlas/fMRI_vtk/landmarks_all.p'
landmarks_from_file = pickle.load(open(landmark_file))

# parameters for clustering/atlas creation
number_of_clusters = 400
#number_of_clusters = 20
number_of_fibers_per_subject = 1000
#number_of_fibers_per_subject = 500
minimum_length = 40
number_of_sampled_fibers = 500
#number_of_sampled_fibers = 100
#number_of_eigenvectors = 15
#number_of_eigenvectors = 5
number_of_eigenvectors = 10
sigma = 50
#sigma = 30

# parameters for atlas labeling
number_of_fibers_to_label = 4000

# read all polydata in the directory
#input_mask = "{0}/*.vtk".format(indir)
input_mask = "{0}/*.vtp".format(indir)
input_poly_datas = glob.glob(input_mask)
#input_poly_datas = input_poly_datas[0:1]
#input_poly_datas = input_poly_datas[0:6]
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
    if (vtk.vtkVersion().GetVTKMajorVersion() > 5.0):
        appender.AddInputData(pd)
    else:
        appender.AddInput(pd)
 
appender.Update()
input_data = appender.GetOutput()

del input_pds

if 0:
    # synthetic landmarks for now
    number_of_landmarks = 10
    #number_of_landmarks = 3
    synthetic_list = list()
    for lidx in range(number_of_landmarks):
        synthetic_list.append(numpy.random.random_sample(3)*180 - 90)
    # landmarks per subject, repeated per-fiber (since fiber order 
    # changes with Nystrom)
    number_of_subjects = len(input_poly_datas)
    landmarks = numpy.zeros([number_of_subjects*number_of_fibers_per_subject, number_of_landmarks, 3])
    for sidx in range(number_of_subjects):
        for lidx in range(number_of_landmarks):
            subject_landmark = synthetic_list[lidx] + numpy.random.random_sample()*10
            for fidx in range(number_of_fibers_per_subject):
                landmarks[sidx*number_of_fibers_per_subject + fidx,lidx,:] = subject_landmark

# real landmarks
number_of_subjects = len(input_poly_datas)
landmarks_from_file = pickle.load(open(landmark_file))
if len(landmarks_from_file) != number_of_subjects:
    print "ERROR: number of subjects in landmark file does not match number of input brains."
number_of_landmarks = len(landmarks_from_file[0])    
landmarks = numpy.zeros([number_of_subjects*number_of_fibers_per_subject, number_of_landmarks, 3])
for sidx in range(len(landmarks_from_file)):
    for lidx in range(number_of_landmarks):
        for fidx in range(number_of_fibers_per_subject):
            landmarks[sidx*number_of_fibers_per_subject + fidx,lidx,:] = landmarks_from_file[sidx][lidx]


# CLUSTERING step
nystrom_mask = numpy.random.permutation(input_data.GetNumberOfLines()) < number_of_sampled_fibers

output_polydata_s, cluster_numbers_s, color, embed, distortion, atlas = \
    wma.cluster.spectral(input_data, number_of_clusters=number_of_clusters, \
                             number_of_jobs=number_of_jobs, use_nystrom=True, \
                             nystrom_mask = nystrom_mask, \
                             number_of_eigenvectors = number_of_eigenvectors, \
                             sigma = sigma, landmarks = landmarks)

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

ren2 = wma.cluster.view_cluster_number(output_polydata_s, cidx, cluster_numbers_s); print cidx; cidx = cidx+1

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
plt.hist(centroid_distance,1000); plt.savefig('centroid_distances.pdf')
fiber_mask = centroid_distance > 3.5
pd_dist = wma.filter.mask(output_polydata_s, fiber_mask, centroid_distance)
ren_dmax = wma.render.render(pd_dist, 1000)


fiber_mask = centroid_distance < .25
pd_dist = wma.filter.mask(output_polydata_s, fiber_mask, centroid_distance)
ren_dmin = wma.render.render(pd_dist, 1000)


# compute closest fiber to each centroid
# and furthest, within cluster
closest_fibers = list()
furthest_fibers = list()
for cidx in range(atlas.centroids.shape[0]):
    cluster = numpy.where(cluster_numbers_s==cidx)[0]
    cluster_diff = atlas.centroids[cidx,:] - embed[cluster,:]
    dist = numpy.sum(numpy.multiply(cluster_diff, cluster_diff), 1)
    closest_fibers.append(cluster[numpy.argmin(dist)])
    furthest_fibers.append(cluster[numpy.argmax(dist)])

closest_fibers = numpy.array(closest_fibers)
furthest_fibers = numpy.array(furthest_fibers)

close_mask = numpy.zeros(embed.shape[0])
close_mask[closest_fibers] = 1
far_mask = numpy.zeros(embed.shape[0])
far_mask[furthest_fibers] = 1

pd_close = wma.filter.mask(output_polydata_s, close_mask, cluster_numbers_s)
pd_far = wma.filter.mask(output_polydata_s, far_mask, cluster_numbers_s)
ren_close = wma.render.render(pd_close, 1000)
ren_far = wma.render.render(pd_far, 1000)

# save all polydatas individually.
# color by subject for now
subject_fiber_list = list()
for sidx in range(number_of_subjects):
    for fidx in range(number_of_fibers_per_subject):
        subject_fiber_list.append(sidx)
subject_fiber_list = numpy.array(subject_fiber_list)

#for cidx in range(atlas.centroids.shape[0]):
#    cluster_mask = (cluster_numbers_s==cidx)
#    pd_cluster = wma.filter.mask(output_polydata_s, cluster_mask, subject_fiber_list)
#    fname = 'cluster_{0}.vtk'.format(cidx)
#    wma.io.write_polydata(pd_cluster, fname)
#    #writer = vtk.vtkPolyDataWriter()
#    #writer.SetFileName(fname)
#    #writer.SetInput(pd_cluster)
#    #writer.Update()
#    #del writer
#    #print cidx


# to view the subjects in each cluster
cidx = 0
cluster_mask = (cluster_numbers_s==cidx); pd_cluster = wma.filter.mask(output_polydata_s, cluster_mask, subject_fiber_list); ren2 = wma.render.render(pd_cluster); print "C:", cidx, 'S:', len(set(subject_fiber_list[cluster_mask])), '/', number_of_subjects; cidx = cidx+1

# figure out how many subjects in most clusters
subjects_per_cluster = list()
for cidx in range(atlas.centroids.shape[0]):
    cluster_mask = (cluster_numbers_s==cidx) 
    subjects_per_cluster.append(len(set(subject_fiber_list[cluster_mask])))

plt.figure()
plt.hist(subjects_per_cluster, number_of_subjects)
plt.savefig('subjects_per_cluster_hist.pdf')
plt.close()


mask = (numpy.array(subjects_per_cluster) < 6)[cluster_numbers_s]
pd_cluster = wma.filter.mask(output_polydata_s, mask, subject_fiber_list)
ren2 = wma.render.render(pd_cluster)

keep = [76, 195, 233, 234, 242, 244, 270, 280, 290, 330, 348, 353, 357, 358, 365, 368, 371, 383, 387]

for cidx in keep:
    cluster_mask = (cluster_numbers_s==cidx)
    pd_cluster = wma.filter.mask(output_polydata_s, cluster_mask, subject_fiber_list);
    fname = 'cluster_{0}.vtk'.format(cidx)
    wma.io.write_polydata(pd_cluster, fname)

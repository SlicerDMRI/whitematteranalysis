#!/usr/bin/env python

import whitematteranalysis as wma
import vtk
import numpy
import matplotlib.pyplot as plt
import multiprocessing
import glob
import sys
import os

number_of_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', number_of_jobs

print 'Read and preprocess'

indir = '/Users/odonnell/Data/tbi_with_scalars'
outdir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/March2013/TBI/cluster_5'

# parameters for clustering/atlas creation
number_of_clusters = 400
#number_of_clusters = 80
number_of_fibers_per_subject = 1000
minimum_length = 40
#number_of_sampled_fibers = 500
number_of_sampled_fibers = 1000
number_of_eigenvectors = 15
sigma = 50

# parameters for atlas labeling
# not used now
#number_of_fibers_to_label = 4000

# read all polydata in the directory
input_mask = "{0}/*.vtp".format(indir)
input_poly_datas = glob.glob(input_mask)

#input_poly_datas = input_poly_datas[0:1]
#input_poly_datas = input_poly_datas[0:6]

if len(input_poly_datas) == 0:
    print ""
    print "ERROR: no polydatas found in input directory:"
    print input_mask
    print "<cluster_atlas> exiting."
    sys.exit(0) 

print input_poly_datas

group_indices =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,\
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
group_indices = numpy.array(group_indices)
input_poly_datas_all = input_poly_datas
#input_poly_datas_all = input_poly_datas[0:1]
# cluster controls only to define regions
#input_poly_datas_controls = input_poly_datas[0:1]
input_poly_datas_controls = input_poly_datas[0:13]

# below this line the pds are read and clustered
number_of_subjects = len(input_poly_datas)

# read in data
input_pds = list()
for fname in input_poly_datas_controls:
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
    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        appender.AddInputData(pd)
    else:
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
ren = wma.render.render(output_polydata_s, 1000)
directory = os.path.join(outdir, 'atlas_embed') 
if not os.path.isdir(directory):
    os.mkdir(directory)
ren.save_views(directory)
del ren

# view cluster numbers
output_polydata_s.GetCellData().SetActiveScalars('ClusterNumber')
ren = wma.render.render(output_polydata_s)
directory = os.path.join(outdir, 'atlas_clusters') 
if not os.path.isdir(directory):
    os.mkdir(directory)
ren.save_views(directory)
del ren

# compute distances to centroids
diff = atlas.centroids[cluster_numbers_s] - embed
centroid_distance = numpy.sum(numpy.multiply(diff,diff), 1)
fiber_mask = centroid_distance > 10.0
pd_dist = wma.filter.mask(output_polydata_s, fiber_mask, centroid_distance)
ren_dmax = wma.render.render(pd_dist, 1000)
plt.hist(centroid_distance,1000); plt.savefig('centroid_distances.pdf')
directory = os.path.join(outdir, 'max_dist') 
if not os.path.isdir(directory):
    os.mkdir(directory)
ren_dmax.save_views(directory)
del ren_dmax

fiber_mask = centroid_distance < 1.0
pd_dist = wma.filter.mask(output_polydata_s, fiber_mask, centroid_distance)
ren_dmin = wma.render.render(pd_dist, 1000)
directory = os.path.join(outdir, 'min_dist') 
if not os.path.isdir(directory):
    os.mkdir(directory)
ren_dmin.save_views(directory)
del ren_dmin

# figure out which subject each fiber was from 
subject_fiber_list = list()
for sidx in range(number_of_subjects):
    for fidx in range(number_of_fibers_per_subject):
        subject_fiber_list.append(sidx)
subject_fiber_list = numpy.array(subject_fiber_list)

# figure out how many subjects in most clusters
subjects_per_cluster = list()
for cidx in range(atlas.centroids.shape[0]):
    cluster_mask = (cluster_numbers_s==cidx) 
    subjects_per_cluster.append(len(set(subject_fiber_list[cluster_mask])))

plt.figure()
plt.hist(subjects_per_cluster, number_of_subjects)
plt.savefig('subjects_per_cluster_hist.pdf')
plt.close()

# Now we want to use this atlas to label all subjects, all fibers...
# read in data and label each subject separately
# separate short fibers
input_pds = list()
for fname in input_poly_datas_all:
    print fname
    pd = wma.io.read_polydata(fname)
    number_of_lines = pd.GetNumberOfLines()
    pd2, line_indices = wma.filter.preprocess(pd, minimum_length, return_indices=True)
    output_polydata_2, cluster_numbers_2, color2, embed2 = \
        wma.cluster.spectral_atlas_label(pd2, atlas)
    final_cluster_id = numpy.zeros(number_of_lines) - 1
    final_cluster_id[line_indices] = cluster_numbers_2
    # put cluster numbers into original input pd
    cluster_colors = vtk.vtkIntArray()
    cluster_colors.SetName('ClusterNumber_CellData')
    for lidx in range(number_of_lines):
        cluster_colors.InsertNextTuple1(final_cluster_id[lidx])
    pd.GetCellData().AddArray(cluster_colors)
    del cluster_colors
    # also put cluster numbers into the point data
    # for viewing in slicer and measuring from tensors
    cluster_colors = vtk.vtkIntArray()
    cluster_colors.SetName('ClusterNumber_PointData')
    ptids = vtk.vtkIdList()
    pd.GetLines().InitTraversal()
    for lidx in range(number_of_lines):
        pd.GetLines().GetNextCell(ptids)
        for pidx in range(ptids.GetNumberOfIds()):
            cluster_colors.InsertNextTuple1(final_cluster_id[lidx])
    pd.GetPointData().AddArray(cluster_colors)
    del cluster_colors

    subj_id = os.path.splitext(os.path.split(fname)[1])[0]
    fname = 'clusters_' + subj_id +'.vtp'
    fname = os.path.join(outdir, fname) 
    wma.io.write_polydata(pd, fname)
    del pd
    del pd2
    #del output_polydata_2

    #outdir2 = os.path.join(outdir, subj_id) 
    #if not os.path.isdir(outdir2):
    #    os.mkdir(outdir2)
    #ren = wma.render.render(pd2, 1000)
    #directory = os.path.join(outdir2, 'label_embed') 
    #if not os.path.isdir(directory):
    #    os.mkdir(directory)
    #ren.save_views(directory)
    #del ren

    #ren = wma.render.render(pd, 1000)
    #directory = os.path.join(outdir2, 'label_clusters') 
    #if not os.path.isdir(directory):
    #    os.mkdir(directory)
    #ren.save_views(directory)
    #del ren



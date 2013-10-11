
import os
import glob
import multiprocessing

import numpy
import vtk
import whitematteranalysis as wma

print 'LAUREN can similarity avoid explicit matrix replication'

print 'Read and preprocess'

minimum_length = 60

#outdir = os.path.join('.', 'test_outputs')
#os.mkdir(outdir)

number_of_fibers = 300
number_of_fibers = 1000

number_of_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', number_of_jobs
#number_of_jobs = 12
#number_of_jobs = 16

input_dir = '/Users/lauren/Dropbox/Coding/OUTPUTS/group_register_test_13s_brage_moredatathreads'
input_mask = "{0}/*.vtk".format(input_dir)
input_pd_fnames = glob.glob(input_mask)
input_pds = list()

# 9 to 11 are above the others. need to center?
#input_pd_fnames = input_pd_fnames[0:8] + input_pd_fnames[11:13]
#input_pd_fnames = input_pd_fnames[0:8] + input_pd_fnames[11:15] + input_pd_fnames[16:17] 
input_pd_fnames = input_pd_fnames[0:5]

print input_pd_fnames

for fname in input_pd_fnames:
    print fname
    pd = wma.io.read_polydata(fname)
    pd2 = wma.filter.preprocess(pd, minimum_length)
    pd3 = wma.filter.downsample(pd2, number_of_fibers)
    input_pds.append(pd3)

# test clustering the data
appender = vtk.vtkAppendPolyData()
for pd in input_pds:
    appender.AddInput(pd)

appender.Update()
pd_all_registered = appender.GetOutput()

# figure out which fibers are from which subject
subject_idx = list()
sidx = 0
for pd in input_pds:
    num_lines = pd.GetNumberOfLines()
    print sidx, num_lines
    for line in range(0, num_lines):
        subject_idx.append(sidx)
    sidx = sidx + 1

total_number_of_fibers = len(subject_idx)

# USE REAL DISTANCES
distances = wma.cluster._pairwise_distance_matrix(pd_all_registered, 0, number_of_jobs)

# try to measure potential cluster quality and find good clusters
# at various size scales
k = 25.0
k = 50.0
fiber_centroid_quality = numpy.zeros(total_number_of_fibers)
fiber_idxs_remaining = numpy.array(range(0, total_number_of_fibers))
fiber_idxs_all = numpy.array(range(0, total_number_of_fibers))
cluster_error = numpy.zeros(total_number_of_fibers)
cluster_error1 = numpy.zeros(total_number_of_fibers)
cluster_error2 = numpy.zeros(total_number_of_fibers)
n_subj = len(input_pds)
expected_percent = 1.0/n_subj

for fidx in fiber_idxs_all:
    # compute its quality as a centroid for this value of k
    print fidx
    # find k nearest neighbors
    dist = distances[fidx,:]
    dist_sorted = numpy.sort(dist)
    # find distances to these fibers. ignore first one, that is to self
    dist_chosen = dist_sorted[1 : k + 1]
    # find these k fibers
    k_fibers = list()
    for d in dist_chosen:
        k_fibers.append(numpy.nonzero(d == dist)[0][0])
    # see how good the cluster is
    # test if the subjects are present, in equal amounts
    subj_count = numpy.zeros(len(input_pds))
    for kidx in k_fibers:
        sidx = subject_idx[kidx]
        subj_count[sidx]  = subj_count[sidx] + 1
    # compare to expected percent of fibers from each subject
    subj_percent = numpy.divide(subj_count, k)
    error1 = expected_percent - subj_percent
    error1 = numpy.sqrt(numpy.sum(error1*error1))
    # also rate cluster by linkage properties
    # what is worst pairwise distance?
    # for test try max dist to centroid not final method, cant compare across k
    error3 = numpy.max(dist_chosen)
    local_dists = list()
    for kidx in range(1, int(k)):
        local_dists.append(dist_chosen[kidx])
    for kidx in k_fibers:
        # this does count the distances each twice, but we are using the max anyway
        for kidx2 in k_fibers:
            if not kidx == kidx2:
                local_dists.append(distances[kidx,kidx2])
    error2 = numpy.max(numpy.array(local_dists))
    # output value
    cluster_error[fidx] = error1 * error2 * error3
    cluster_error1[fidx] = error1
    cluster_error2[fidx] = error2


nf = pd_all_registered.GetNumberOfLines()
print nf
mask = numpy.ones(nf)
pd_quality = wma.filter.mask(pd_all_registered, mask, cluster_error)
#ren = wma.render.render(pd_quality, scalar_bar=True)

#mask = cluster_error < 0.5*numpy.median(cluster_error)
#pd_quality_masked = wma.filter.mask(pd_all_registered, mask, cluster_error)

# try to assign to clusters greedily for this k
cluster_idx = numpy.zeros(total_number_of_fibers) -1
cluster_unassigned = numpy.ones(total_number_of_fibers) 

# sort cluster_error and find indices into it
# this could cause a bug if there are duplicate entries fix
cq_sort = numpy.sort(cluster_error)
cq_indices = list()
for c in cq_sort:
    cq_indices.append(numpy.nonzero(c==cluster_error)[0][0])

#cq_indices = numpy.array(cq_indices)
centroid_idx = list()

# BUG?  what if k neighbors were already in another cluster?
# can this happen? check for fiber already assigned messages
for c in cq_indices:
    print c
    # if this centroid is still available make its cluster
    if (cluster_unassigned[c]):
        print 'CENTROID', c
        centroid_idx.append(c)
        # find members of the cluster
        # find k nearest neighbors
        dist = distances[c,:]
        dist_sorted = numpy.sort(dist)
        # find distances to these fibers. ignore first one, that is to self
        dist_chosen = dist_sorted[1 : k + 1]
        # find these k fibers
        k_fibers = list()
        for d in dist_chosen:
            k_fibers.append(numpy.nonzero(d == dist)[0][0])
        # label k cluster members
        # remove members from further search
        for kidx in k_fibers:
            if cluster_unassigned[kidx]:
                cluster_idx[kidx] = c
                cluster_unassigned[kidx] = 0
            else:
                print "fiber already assigned"
            

ren_list=list()
for cidx in centroid_idx[1:10]:
    print cidx
    ren = wma.cluster.view_cluster_number(pd_all_registered,cidx, cluster_idx)
    ren_list.append(ren)
    tmplist = list()
    for tmpidx in numpy.nonzero(cidx==cluster_idx)[0]:
        tmplist.append(subject_idx[tmpidx])
    print "CLUSTER", cidx, "SUBJECTS", tmplist
    del ren

#del ren_list
print "del ren_list to kill renderers"

appender = vtk.vtkAppendPolyData()
colors = numpy.ones(len(cluster_error))
color_idx = 0

#for cidx in centroid_idx[400:-1]:
#for cidx in centroid_idx[0:50]:
for cidx in centroid_idx[101:150]:
    print cidx
    pd = wma.filter.mask(pd_all_registered, cidx==cluster_idx, colors*color_idx)
    appender.AddInput(pd)
    color_idx = color_idx + 1

appender.Update()

top_clusters = appender.GetOutput()






### LAUREN you are here == make this work multiscale saturday
nf = pd_all_registered.GetNumberOfLines()
print nf
mask = numpy.ones(nf)
mce = numpy.mean(cluster_error,1)
pd_quality = wma.filter.mask(pd_all_registered, mask, mce)
#pd_quality = wma.filter.mask(pd_all_registered, mask, cluster_error)
#ren = wma.render.render(pd_quality, scalar_bar=True)

#mask = cluster_error < 0.5*numpy.median(cluster_error)
#pd_quality_masked = wma.filter.mask(pd_all_registered, mask, cluster_error)

cluster_error = (1 + cluster_error1 )* cluster_error2 * cluster_error3 / cluster_error4

maxerr = numpy.max(cluster_error)
badidx = numpy.nonzero(cluster_error4 < 0.80)
cluster_error[badidx] = cluster_error[badidx] + maxerr

badidx = cluster_error2 >= 5.0
cluster_error[badidx] = cluster_error[badidx] + maxerr

badidx = cluster_error3 >= 40.0
cluster_error[badidx] = cluster_error[badidx] + maxerr

k_indices = list()
min_ce = numpy.min(cluster_error, 1)
for idx in range(0, len(min_ce)):
    k_indices.append(numpy.nonzero(min_ce[idx]== cluster_error[idx,:])[0][0])

bigclusters = cluster_error2 < 15.0
num_rows, c = bigclusters.shape
for r in range(0, num_rows):
    min_ce[r] = min(cluster_error[r]) + maxerr
    # change k index?? huh.
    if bigclusters[r].any():
        idx = sum(bigclusters[r]) - 1
        min_ce[r] = cluster_error[r, idx]
        k_indices[r] = idx
    
#and (cluster_error4[c, k_indices[c]] >= 0.75) and (cluster_error3[c, k_indices[c]] <=10.0):
 
# sort cluster_error and find indices into it
# this could cause a bug if there are duplicate entries fix
#cq_sort = numpy.sort(cluster_error)
# LAUREN really look at scales separately
cq_sort = numpy.sort(min_ce)
cq_indices = list()
for c in cq_sort:
    cq_indices.append(numpy.nonzero(c==min_ce)[0][0])
    #cq_indices.append(numpy.nonzero(c==cluster_error)[0][0])

#cq_indices = numpy.array(cq_indices)
# try to assign to clusters greedily for this k
cluster_idx = numpy.zeros(total_number_of_fibers) -1
cluster_unassigned = numpy.ones(total_number_of_fibers)
centroid_idx = list()
cluster_sizes = list()
cluster_subjects = list()

# BUG?  what if k neighbors were already in another cluster?
# can this happen? check for fiber already assigned messages
for c in cq_indices:
    print c
    # if this centroid is still available make its cluster
    # only do this if more than x% of subjects are present
    # otherwise very very tight clusters from 1 subject can be found
    if (cluster_unassigned[c]):
        print 'CENTROID', c
        cluster_unassigned[c] = 0
        cluster_idx[c] = c
        centroid_idx.append(c)
        cluster_size = 1
        subjects_in_cluster = numpy.zeros(len(input_pds))
        subjects_in_cluster[subject_idx[c]]  = 1
                
        # find members of the cluster, k nearest neighbors
        dist = distances[c,:]
        dist_sorted = numpy.sort(dist)
        # find distances to these fibers. ignore first one, that is to self
        #k = k_list[k_indices[c]]
        k = k_later_sampling[k_indices[c]]
        #k = k_list[k_indices[c]] + 15
        dist_chosen = dist_sorted[1 : k + 1]
        # find these k fibers
        k_fibers = list()
        for d in dist_chosen:
            k_fibers.append(numpy.nonzero(d == dist)[0][0])
        # label k cluster members and remove members from further search
        for kidx in k_fibers:
            if cluster_unassigned[kidx]:
                cluster_idx[kidx] = c
                cluster_unassigned[kidx] = 0
                cluster_size = cluster_size + 1
                sidx = subject_idx[kidx]
                subjects_in_cluster[sidx]  = 1
            else:
                print "fiber already assigned"
        cluster_sizes.append(cluster_size)
        cluster_subjects.append(numpy.sum(subjects_in_cluster))

print "UNASSIGNED FIBERS: ", numpy.sum(cluster_unassigned)

ren_list=list()
for cidx in centroid_idx[1:10]:
    print cidx
    ren = wma.cluster.view_cluster_number(pd_all_registered, cidx, cluster_idx)
    ren_list.append(ren)
    tmplist = list()
    for tmpidx in numpy.nonzero(cidx==cluster_idx)[0]:
        tmplist.append(subject_idx[tmpidx])
    print "CLUSTER", cidx, "SUBJECTS", tmplist
    del ren

#del ren_list
print "del ren_list to kill renderers"

# find non-singleton, non-tiny clusters
keepers_size = (numpy.array(cluster_sizes) > 5)
keepers_subjects = (numpy.array(cluster_subjects) > 10)
keepers = keepers_size.astype(int) * keepers_subjects.astype(int)

#keepers = (numpy.array(cluster_sizes) > 5)

centroid_idx_keepers = numpy.array(centroid_idx)[numpy.nonzero(keepers)[0]]
cluster_sizes_keepers = numpy.array(cluster_sizes)[numpy.nonzero(keepers)[0]]
cluster_subjects_keepers = numpy.array(cluster_subjects)[numpy.nonzero(keepers)[0]]

# see how many subjects are in those good ones
#for cidx in centroid_idx[400:-1]:
#for cidx in centroid_idx[101:150]:
#
#for cidx in centroid_idx[1:20]:
#for cidx in centroid_idx[21:100]:
#for cidx in centroid_idx[101:200]:
#for cidx in centroid_idx_keepers[61:90]:
ren_list=list()
colors = numpy.ones(len(cluster_error))
idx = 0
while ( idx < len(centroid_idx_keepers) -10):
    appender = vtk.vtkAppendPolyData()
    color_idx = 0
    for cidx in centroid_idx_keepers[idx:idx + 10]:
        print cidx
        pd = wma.filter.mask(pd_all_registered, cidx==cluster_idx, colors*color_idx)
        appender.AddInput(pd)
        color_idx = color_idx + 1
    appender.Update()
    top_clusters = appender.GetOutput()
    ren = wma.render.render(top_clusters)
    ren_list.append(ren)
    idx = idx + 10

appender = vtk.vtkAppendPolyData()
color_idx = 0
for cidx in centroid_idx_keepers:
    pd = wma.filter.mask(pd_all_registered, cidx==cluster_idx, colors*color_idx)
    appender.AddInput(pd)
    color_idx = color_idx + 1
appender.Update()
ren = wma.render.render(appender.GetOutput())

print "Good clusters:", len(centroid_idx_keepers)
print "Subjects per cluster:", numpy.min(cluster_subjects_keepers), numpy.max(cluster_subjects_keepers), numpy.mean(cluster_subjects_keepers)
print "Sizes per cluster:", numpy.min(cluster_sizes_keepers), numpy.max(cluster_sizes_keepers), numpy.mean(cluster_sizes_keepers)
print "Fibers in good clusters:", numpy.sum(cluster_sizes_keepers), '/', numpy.sum(cluster_sizes), '=', \
    numpy.sum(cluster_sizes_keepers)/float(numpy.sum(cluster_sizes))

kidx = 0
centroid_ks = list()
centroid_keeper_ks = list()
cluster_ks = list()
for k in k_list:
    centroid_keeper_ks.append(sum(numpy.array(k_indices)[centroid_idx_keepers] == kidx))
    centroid_ks.append(sum(numpy.array(k_indices)[centroid_idx] == kidx))
    cluster_ks.append(sum(numpy.array(k_indices) == kidx))
    kidx = kidx + 1

print k_list
print cluster_ks
print centroid_ks
print centroid_keeper_ks



# try to measure potential cluster quality and find good clusters
# at various size scales
#k = 25.0
#k = 50.0
# scales to test
k_list = [10, 25, 50, 100]
number_k = len(k_list)

number_of_subjects = len(input_pds)
expected_percent = 1.0/number_of_subjects

# outputs measuring quality of putative cluster at each centroid
cluster_error = numpy.zeros([total_number_of_fibers, number_k])
cluster_error1 = numpy.zeros([total_number_of_fibers, number_k])
cluster_error2 = numpy.zeros([total_number_of_fibers, number_k])
cluster_error3 = numpy.zeros([total_number_of_fibers, number_k])
current_k_idx = 0

#LAUREN should loop over fidx, then sort will happen in outer loop

for k in k_list:
    fiber_idxs_remaining = numpy.array(range(0, total_number_of_fibers))
    fiber_idxs_all = numpy.array(range(0, total_number_of_fibers))
    for fidx in fiber_idxs_all:
        # compute its quality as a centroid for this value of k
        if numpy.mod(fidx, 250) == 0:
            print 'K:', k, '(', current_k_idx + 1, '/', number_k, ') F:', fidx
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
        for kidx1 in range(1, int(k)):
            local_dists.append(dist_chosen[kidx1])
        for kidx1 in k_fibers:
            # this does count the distances each twice, but we are using the max anyway
            for kidx2 in k_fibers:
                if not kidx1 == kidx2:
                    local_dists.append(distances[kidx1,kidx2])
        error2 = numpy.max(numpy.array(local_dists))
        # output value
        #cluster_error[fidx, current_k_idx] = error1 * error2 * error3
        cluster_error1[fidx, current_k_idx] = error1
        cluster_error2[fidx, current_k_idx] = error2
        cluster_error3[fidx, current_k_idx] = error3
    current_k_idx = current_k_idx + 1

cluster_error = cluster_error1 * cluster_error2 * cluster_error3

# this below works terribly. fibers are unlike centroid.
#cluster_error = cluster_error1 * cluster_error2

print "Lauren if ERROR 1 == 0 error is 0, not good."

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

k_indices = list()
min_ce = numpy.min(cluster_error, 1)
for idx in range(0, len(min_ce)):
    k_indices.append(numpy.nonzero(min_ce[idx]== cluster_error[idx,:])[0][0])


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
    if (cluster_unassigned[c]):
        print 'CENTROID', c
        cluster_unassigned[c] = 0
        cluster_idx[c] = c
        centroid_idx.append(c)
        cluster_size = 1
        # find members of the cluster, k nearest neighbors
        dist = distances[c,:]
        dist_sorted = numpy.sort(dist)
        # find distances to these fibers. ignore first one, that is to self
        k = k_list[k_indices[c]]
        #k = k_list[k_indices[c]] + 15
        dist_chosen = dist_sorted[1 : k + 1]
        # find these k fibers
        k_fibers = list()
        for d in dist_chosen:
            k_fibers.append(numpy.nonzero(d == dist)[0][0])
        # label k cluster members and remove members from further search
        subj_present = numpy.zeros(len(input_pds))
        for kidx in k_fibers:
            if cluster_unassigned[kidx]:
                cluster_idx[kidx] = c
                cluster_unassigned[kidx] = 0
                cluster_size = cluster_size + 1
                sidx = subject_idx[kidx]
                subj_present[sidx]  = 1
            else:
                print "fiber already assigned"
        cluster_sizes.append(cluster_size)
        cluster_subjects.append(numpy.sum(subj_present))

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

# find non-singleton, non-tiny clusters
keepers = numpy.array(cluster_sizes)>10
centroid_idx_keepers = numpy.array(centroid_idx)[numpy.nonzero(keepers)[0]]

appender = vtk.vtkAppendPolyData()
colors = numpy.ones(len(cluster_error))
color_idx = 0
#for cidx in centroid_idx[400:-1]:
#for cidx in centroid_idx[101:150]:
#
#for cidx in centroid_idx[1:20]:
#for cidx in centroid_idx[21:100]:
#for cidx in centroid_idx[101:200]:
#for cidx in centroid_idx_keepers[61:90]:
for cidx in centroid_idx_keepers[61:70]:
    print cidx
    pd = wma.filter.mask(pd_all_registered, cidx==cluster_idx, colors*color_idx)
    appender.AddInput(pd)
    color_idx = color_idx + 1

appender.Update()

top_clusters = appender.GetOutput()

ren = wma.render.render(top_clusters)


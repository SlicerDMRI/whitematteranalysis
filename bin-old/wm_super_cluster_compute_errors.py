
# try to measure potential cluster quality and find good clusters
# at various size scales
#k = 25.0
#k = 50.0
# scales to test
# these values work for 5 subjects, 1000 fibers each
#k_list = [10, 25, 50, 100]
#k_later_sampling = [17, 37, 75, 125]
k_list = [10, 12, 15, 17, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100]
#k_later_sampling = [11, 13, 18, 23, 35, 75, 125]
k_later_sampling = k_list
number_k = len(k_list)

number_of_subjects = len(input_pds)
expected_percent = 1.0/number_of_subjects

# outputs measuring quality of putative cluster at each centroid
cluster_error = numpy.zeros([total_number_of_fibers, number_k])
cluster_error1 = numpy.zeros([total_number_of_fibers, number_k])
cluster_error2 = numpy.zeros([total_number_of_fibers, number_k])
cluster_error3 = numpy.zeros([total_number_of_fibers, number_k])
cluster_error4 = numpy.zeros([total_number_of_fibers, number_k])
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
        error4 = numpy.sum(subj_count > 0) / float(number_of_subjects)
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
        cluster_error4[fidx, current_k_idx] = error4
    current_k_idx = current_k_idx + 1

cluster_error = (1 + cluster_error1 )* cluster_error2 * cluster_error3 / cluster_error4

# this also works poorly
#cluster_error = cluster_error1 * cluster_error3

# this below works terribly. fibers are unlike centroid.
#cluster_error = cluster_error1 * cluster_error2

print "Lauren if ERROR 1 == 0 error is 0, not good."

import glob
import os
import numpy
import scipy.stats
import vtk
import whitematteranalysis as wma


# ================= FUNCTIONS (move to module code) =======================
# create some polydata objects to view the results
def fiber_list_to_fiber_array(fiber_list):
    fiber_array = wma.fibers.FiberArray()    
    fiber_array.number_of_fibers = len(fiber_list)
    fiber_array.points_per_fiber = len(fiber_list[0].r)
    dims = [fiber_array.number_of_fibers, fiber_array.points_per_fiber]
    # fiber data
    fiber_array.fiber_array_r = numpy.zeros(dims)
    fiber_array.fiber_array_a = numpy.zeros(dims)
    fiber_array.fiber_array_s = numpy.zeros(dims)
    curr_fidx = 0
    for curr_fib in fiber_list:
        fiber_array.fiber_array_r[curr_fidx] = curr_fib.r
        fiber_array.fiber_array_a[curr_fidx] = curr_fib.a
        fiber_array.fiber_array_s[curr_fidx] = curr_fib.s
        curr_fidx += 1
    return fiber_array

def add_array_to_polydata(pd, array, array_name='Test', array_type='Cell'):
    out_array = vtk.vtkFloatArray()
    for idx in range(len(array)):
        out_array.InsertNextTuple1(array[idx])
    out_array.SetName(array_name)
    ret = pd.GetCellData().AddArray(out_array)
    print ret
    pd.GetCellData().SetActiveScalars(array_name)
    return(pd)

def find_islands_and_sizes(distances, distance_threshold, statistic_mask):
    # input distances between sample and all other fibers
    [number_of_fiber_centroids, total_number_of_fibers] = distances.shape
    # initialize data structures
    island = numpy.zeros([number_of_fiber_centroids,1]) - 100
    curr_island = 0
    # loop over centroids, agglomerating nearest neighborhoods
    for cidx in range(number_of_fiber_centroids):
        #print cidx
        # was if statistic < threshold.
        # but in some cases (KL) want if statistic > threshold.
        # so just input a mask of whether or not the fiber passes the test
        if statistic_mask[cidx] == True:
            # start a new island if needed
            if island[cidx] == -100:
                island[cidx] = curr_island
                curr_island += 1
            #print cidx, island[cidx]
            # find nearest neighbors that are centroids. skip self (first one)
            centroids_sorted = numpy.argsort(distances[cidx,:])[1:]
            growing = 1
            candidate = 0
            number_of_candidates = len(centroids_sorted)
            while growing & (candidate < number_of_candidates):
                cidx2 = centroids_sorted[candidate]
                # both statistic and distance must be over threshold to stop agglomerating
                if distances[cidx, cidx2] > distance_threshold:
                    growing = 0
                # if this one can be agglomerated with the current island do it
                elif statistic_mask[cidx2] == True:
                    # if this one is unassigned so far assign to the current island
                    if island[cidx2] == -100:
                        #print "join: ", island[cidx2], island[cidx], cidx, cidx2
                        island[cidx2] = island[cidx]
                    elif island[cidx2] != island[cidx]:
                        # otherwise we have to agglomerate the existing island with this one.
                        merge_idx = island[cidx2]
                        #print "agglomerate: ", merge_idx, island[cidx], cidx, cidx2
                        island[numpy.nonzero(island == merge_idx)] = island[cidx]
                        #else:
                        #if distances[cidx, cidx2] > distance_threshold:
                        #growing = 0
                candidate += 1
    # compute sizes
    island_ids, island_indices = numpy.unique(island, return_index=True)
    print "Number of islands: ", len(island_ids)
    island_sizes = list()
    island_ids = island[island_indices]
    for island_id in island_ids:
        island_sizes.append(len(numpy.nonzero(island == island_id)[0]))
    return island, numpy.array(island_ids), numpy.array(island_sizes)

# compute average brain fibers list and polydata
# note sigma around 10.0 preserves structure.
def compute_avg_brain(fiber_array, distances, sigma):
    sigma_sq = sigma*sigma
    # ignore fibers whose weight in the sum is less than 0.018
    neighborhood_threshold = sigma*2.0

    average_fibers_list = list()
    weights_list = list()
    for idx in range(number_of_fiber_centroids):
        if idx % 100 == 0:
            print idx, '/', number_of_fiber_centroids
        neighborhood_indices = numpy.nonzero(distances[idx,:] < neighborhood_threshold)[0]
        d_sq = numpy.multiply(distances[idx,neighborhood_indices], distances[idx,neighborhood_indices])
        neighborhood_weights = numpy.exp(numpy.divide(-d_sq, sigma_sq))
        hood_count = len(neighborhood_indices)
        # find average fiber in neighborhood
        avg_fiber = fiber_array.get_fiber(neighborhood_indices[0])
        avg_fiber *= neighborhood_weights[0]
        for idx in range(1, hood_count):
            curr_fiber = fiber_array.get_fiber(neighborhood_indices[idx])
            curr_fiber *= neighborhood_weights[idx]
            avg_fiber += curr_fiber
        total_weights = numpy.sum(neighborhood_weights)
        avg_fiber /= total_weights
        average_fibers_list.append(avg_fiber)
        weights_list.append(total_weights)

    # output as polydata
    outpd = fiber_list_to_fiber_array(average_fibers_list).convert_to_polydata()
    outpd = add_array_to_polydata(outpd, weights_list, array_name='LocalWeights')
    return outpd, numpy.array(weights_list)

# compute statistics about each local fiber neighborhood
def compute_local_stats(distances, distances_sq, sigma, group_indices, input_group_idx_per_fiber, input_subject_idx_per_fiber):
    sigma_sq = sigma*sigma
    # ignore fibers whose weight in the sum is less than 0.018
    neighborhood_threshold = sigma*2.0
    # groups
    g0 = numpy.nonzero(group_indices == 0)
    g1 = numpy.nonzero(group_indices == 1)
    number_of_subjects = len(group_indices)
    fibers_per_subject = numpy.divide(len(input_group_idx_per_fiber), number_of_subjects)
    fibers_per_group = numpy.divide(len(input_group_idx_per_fiber), 2)
    
    # compute stats
    fiber_count_p = numpy.ones(number_of_fiber_centroids)
    fiber_count_mean_g0 = numpy.zeros(number_of_fiber_centroids)
    fiber_count_std_g0 = numpy.zeros(number_of_fiber_centroids)
    fiber_count_mean_g1 = numpy.zeros(number_of_fiber_centroids)
    fiber_count_std_g1 = numpy.zeros(number_of_fiber_centroids)
    
    fiber_probability_p  = numpy.ones(number_of_fiber_centroids)
    fiber_probability_mean_g0 = numpy.zeros(number_of_fiber_centroids)
    fiber_probability_std_g0 = numpy.zeros(number_of_fiber_centroids)
    fiber_probability_mean_g1 = numpy.zeros(number_of_fiber_centroids)
    fiber_probability_std_g1 = numpy.zeros(number_of_fiber_centroids)

    local_number_subjects = numpy.zeros(number_of_fiber_centroids)

    KL = numpy.ones(number_of_fiber_centroids)
    BC = numpy.ones(number_of_fiber_centroids)
    DIFF2 = numpy.ones(number_of_fiber_centroids)
    DIFF = numpy.ones(number_of_fiber_centroids)
     
    for idx in range(number_of_fiber_centroids):
        if idx % 100 == 0:
            print idx, '/', number_of_fiber_centroids
        neighborhood_indices = numpy.nonzero(distances[idx,:] < neighborhood_threshold)[0]
        # if remove self, error if neighborhood empty otherwise. think about this.
        # use near enough only (for speed) and not self (to remove bias) fibers
        #mask = (distances[idx,:] < neighborhood_threshold) & (distances[idx,:] > 0.0)
        #neighborhood_indices = numpy.nonzero(mask)[0]
        d_sq = distances_sq[idx,neighborhood_indices]
    
        # find neighborhood info per subject and group
        group_list = list()
        subject_list = list()
        for hood_idx in neighborhood_indices:
            #data_list.append(input_data_per_fiber[hood_idx])
            group_list.append(input_group_idx_per_fiber[hood_idx])
            subject_list.append(input_subject_idx_per_fiber[hood_idx])

        subject_list = numpy.array(subject_list)
        group_list = numpy.array(group_list)

        # calculate hard fiber counts
        subject_counts = numpy.bincount(subject_list, minlength=number_of_subjects)
        [t, p] = scipy.stats.ttest_ind(subject_counts[g0], subject_counts[g1])
        fiber_count_p[idx] = p
        fiber_count_mean_g0[idx] = numpy.mean(subject_counts[g0])
        fiber_count_std_g0[idx] =  numpy.std(subject_counts[g0])
        fiber_count_mean_g1[idx] = numpy.mean(subject_counts[g1])
        fiber_count_std_g1[idx] = numpy.std(subject_counts[g1])
        fiber_count_retval = [fiber_count_p, fiber_count_mean_g0, fiber_count_std_g0, fiber_count_mean_g1, fiber_count_std_g1]
        # calculate number of subjects with fibers present (enable removal of outliers)
        # this works for outlier removal only with a relatively lower sigma
        # otherwise can remove low probability fibers, for example
        local_number_subjects[idx] = numpy.sum(subject_counts > 0)
        
        # calculate soft fiber counts (probability density at this centroid)
        fiber_probabilities = numpy.exp(numpy.divide(-d_sq, sigma_sq))
        subject_probabilities = numpy.zeros(number_of_subjects)
        logit_subject_probabilities = numpy.zeros(number_of_subjects)
        for sidx in range(number_of_subjects):
            # sum fiber probabilities from this subject
            sfibers = numpy.nonzero(subject_list == sidx)[0]
            subject_probabilities[sidx] = numpy.sum(fiber_probabilities[sfibers])

        # normalization (size of the set of reference fibers) is the number of fibers per subject
        subject_probabilities = numpy.divide(subject_probabilities, fibers_per_subject)
        # logit transform to make 0 to 1 numbers more appropriate for t-test
        eps = 0.000000001
        logit_subject_probabilities = numpy.log(numpy.divide(subject_probabilities, 1 - subject_probabilities + eps) + eps)
        # multiply by 100 to use probabilities, a bit easier numbers to view in output
        subject_probabilities = subject_probabilities * 100.0

        [t, p] = scipy.stats.ttest_ind(logit_subject_probabilities[g0], logit_subject_probabilities[g1])
        fiber_probability_p[idx] = p
        fiber_probability_mean_g0[idx] = numpy.mean(subject_probabilities[g0])
        fiber_probability_std_g0[idx] =  numpy.std(subject_probabilities[g0])
        fiber_probability_mean_g1[idx] = numpy.mean(subject_probabilities[g1])
        fiber_probability_std_g1[idx] = numpy.std(subject_probabilities[g1])
        fiber_probability_retval = [fiber_probability_p, fiber_probability_mean_g0, fiber_probability_std_g0, fiber_probability_mean_g1, fiber_probability_std_g1]

        # now KL divergence: sum probability in each group and compute: (p-q)(lnp-lnq)
        # normalization (size of the set of reference fibers) is the number of fibers per group 
        eps = 0.000000001
        # don't normalize twice, already normalized by fibers per subject
        #p0 = numpy.divide(numpy.sum(subject_probabilities[g0]), fibers_per_group)
        #p1 = numpy.divide(numpy.sum(subject_probabilities[g1]), fibers_per_group)
        p0 = numpy.sum(subject_probabilities[g0])
        p1 = numpy.sum(subject_probabilities[g1])        
        KL[idx] = numpy.multiply((p0 - p1 + eps),(numpy.log(p0 + eps) - numpy.log(p1 + eps)))
        
        # Bhattacharyya coefficient
        # sum of sqrt of p*q
        BC[idx] = numpy.sqrt(numpy.multiply(p0, p1))


        # simple subtraction (p -q)^2
        tmp = p0 - p1
        DIFF[idx] = numpy.multiply(tmp, tmp)
        # does each group represent 50% here??
        diff_50 = numpy.divide(p0, p0 + p1) - 0.5
        # are subjects evenly dividing the group density?
        #eps = 0
        g0_err = numpy.divide(subject_probabilities[g0], numpy.sum(subject_probabilities[g0]) + eps)
        g1_err = numpy.divide(subject_probabilities[g1], numpy.sum(subject_probabilities[g1]) + eps)
        # subtract expected value if evenly divided among subjects
        expected_val = numpy.divide(2.0, float(number_of_subjects))
        g0_err = g0_err - expected_val
        g1_err = g1_err - expected_val
        #print g0_err
        #print g1_err
        # numerator is large when groups differ from each representing half
        #numerator = numpy.multiply(diff_50, diff_50)
        numerator = numpy.abs(diff_50)

        # square and sum errors, put in the denominator to decrease test statistic if not trustworthy
        #denominator = numpy.sum(numpy.multiply(g0_err, g0_err)) + numpy.sum(numpy.multiply(g1_err, g1_err)) 
        denominator = numpy.sqrt(numpy.sum(numpy.multiply(g0_err, g0_err)) + numpy.sum(numpy.multiply(g1_err, g1_err)) )

        #DIFF2[idx] = numpy.divide(numpy.multiply(diff_50, diff_50), denominator + eps)
        #DIFF2[idx] = numpy.divide(numerator, numpy.multiply(number_of_subjects*number_of_subjects, denominator))
        DIFF2[idx] = numpy.divide(numerator, denominator)
        
    return [fiber_count_retval, fiber_probability_retval, KL, BC, DIFF2, DIFF, local_number_subjects]


def suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic_mask):
    # now find max suprathreshold cluster size
    island_labels, island_ids, island_sizes = find_islands_and_sizes(distances_centroids, distance_threshold, statistic_mask)
    non_background_island_sizes = island_sizes[numpy.nonzero(island_ids != -100)[0]]
    # if there were any significant islands append the maximum size
    if len(non_background_island_sizes) > 0:
        return numpy.max(non_background_island_sizes)
    else:
        # otherwise append 0 for no significant islands
        return 0

def island_correct(distances_centroids, distance_threshold, statistic_mask, null):
    # now find max suprathreshold cluster size
    island_labels, island_ids, island_sizes = find_islands_and_sizes(distances_centroids, distance_threshold, statistic_mask)
    # find significance values for each island size
    background_island = numpy.nonzero(island_ids == -100)[0][0]
    island_p = list()
    for sz in island_sizes:
        island_p.append(len(numpy.nonzero(null >= sz)[0])/float(number_of_permutations))
    island_p = numpy.array(island_p)
    island_p[background_island] = 1.0
    print "P::::", island_p
    # find significance values for each fiber too
    corrected_p = numpy.ones(len(island_labels))
    for is_idx in range(len(island_ids)):
        fiber_idx = numpy.nonzero(island_ids[is_idx] == island_labels)[0]
        corrected_p[fiber_idx] = island_p[is_idx]

    return island_p, corrected_p, island_labels, island_ids, island_sizes

def visualize_data_array_to_disk(polydata, data_array, array_name, outdir, scalar_range=None, p_values=False, p_threshold=0.05):
    add_array_to_polydata(polydata, data_array, array_name=array_name)
    mean_wm.GetCellData().SetActiveScalars(array_name)
    if not os.path.exists(outdir):
            os.makedirs(outdir)
    if p_values == True:
        mask_colors = numpy.array(data_array)
        mask_colors[mask_colors >= p_threshold] = numpy.nan
        polydata = wma.filter.mask(polydata, numpy.ones(len(mask_colors)), mask_colors)
        ren = wma.render.render(polydata, scalar_bar=True, scalar_range=[0.0, p_threshold], colormap='hot')
        actors = ren.renderer.GetActors()
        actors.GetLastActor().GetMapper().GetLookupTable().SetNanColor(0.8,0.8,0.8,0.2)
    elif scalar_range == None:
        ren = wma.render.render(polydata, scalar_bar=True)
    else:
        ren = wma.render.render(polydata, scalar_bar=True, scalar_range=scalar_range)

    ren.save_views(outdir)
    del ren

def do_permutation_NOOP(pidx, number_of_permutations, group_indices_real, number_of_subjects, number_of_fibers_per_subject, distances, distances_sq, neighborhood_sigma, statistic_threshold, distance_threshold, kl_threshold, bc_threshold):
    print "permutation: ", pidx, "/", number_of_permutations
    return 0

def do_permutation(pidx, number_of_permutations, group_indices_real, number_of_subjects, number_of_fibers_per_subject, distances, distances_sq, neighborhood_sigma, statistic_threshold, distance_threshold, kl_threshold, bc_threshold):

    print "==========================================================================="
    print "permutation: ", pidx, "/", number_of_permutations
    print "==========================================================================="    
    # assign groups for each subject and fiber
    group_indices = numpy.random.permutation(group_indices_real)
    input_subject_idx_per_fiber = list()
    input_group_idx_per_fiber = list()
    for subj_idx in range(number_of_subjects):
        for fiber_idx in range(number_of_fibers_per_subject):
            input_subject_idx_per_fiber.append(subj_idx)
            input_group_idx_per_fiber.append(group_indices[subj_idx])
    g0 = numpy.nonzero(group_indices == 0)
    g1 = numpy.nonzero(group_indices == 1)

    [count_info, prob_info, KL, BC, DIFF2, DIFF, subject_count] = compute_local_stats(distances, distances_sq, neighborhood_sigma, group_indices, input_group_idx_per_fiber, input_subject_idx_per_fiber)

    statistic_mask = count_info[0] < statistic_threshold
    null_count = suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic_mask)
    statistic_mask = prob_info[0] < statistic_threshold
    null_prob = suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic_mask)
    statistic_mask = KL > kl_threshold
    null_kl = suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic_mask)
    null_kl_total = numpy.sum(KL)

    #statistic_mask = BC < bc_threshold    
    #null_bc = suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic_mask)
    # do not use, does not make sense, is just high when high p... not appropriate really
    null_bc = 100
    null_bc_total = numpy.sum(BC)

    statistic_mask = DIFF2 > diff2_threshold    
    null_diff2 = suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic_mask)
    null_diff2_total = numpy.sum(DIFF2)

    statistic_mask = DIFF > diff_threshold    
    null_diff = suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic_mask)
    null_diff_total = numpy.sum(DIFF)

    return [null_count, null_prob, null_kl, null_kl_total, null_bc, null_bc_total, null_diff2, null_diff2_total, null_diff, null_diff_total]

# ================= END FUNCTIONS (move to module code) =======================

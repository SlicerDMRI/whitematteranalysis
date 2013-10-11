execfile('/Users/odonnell/Dropbox/Coding/Python/whitematteranalysis/bin/publications/ipmi-stats-paper/stats-computation-code.py')

outdir = '.'
neighborhood_sigma = 20.0
group_indices_real = numpy.array([0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1])
# 20 choose 10 is: 184756
#number_of_permutations = 1000
number_of_permutations = 200
#number_of_permutations = 600
# parameters for island computation
statistic_threshold = 0.05
kl_threshold = 2e-05
distance_threshold = 10.0
# ================= END OF SETTINGS =================

# compute mean brain fibers
# ================= 
print "Computing mean white matter"
#mean_wm, mean_wm_weights = compute_avg_brain(fiber_array, distances, neighborhood_sigma)

# ================= neighborhood parameters
sigma_sq = neighborhood_sigma*neighborhood_sigma

# ================= 
# COMPUTE NULL DISTRIBUTION
# =================
print "Computing null distribution"
null_count = numpy.zeros(number_of_permutations)
null_prob = numpy.zeros(number_of_permutations)
null_kl = numpy.zeros(number_of_permutations)
null_kl_total = numpy.zeros(number_of_permutations)

for pidx in range(number_of_permutations):
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

    [count_info, prob_info, KL, subject_count] = compute_local_stats(distances, distances_sq, neighborhood_sigma, group_indices, input_group_idx_per_fiber, input_subject_idx_per_fiber)

    statistic_mask = count_info[0] < statistic_threshold
    null_count[pidx] = suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic_mask)
    statistic_mask = prob_info[0] < statistic_threshold
    null_prob[pidx] = suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic_mask)
    statistic_mask = KL > kl_threshold
    null_kl[pidx] = suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic_mask)
    null_kl_total[pidx] = numpy.sum(KL)

# ================= 
# COMPUTE SIGNIFICANCES
# ================= 
# assign groups for each subject and fiber
group_indices = group_indices_real
input_subject_idx_per_fiber = list()
input_group_idx_per_fiber = list()
for subj_idx in range(number_of_subjects):
    for fiber_idx in range(number_of_fibers_per_subject):
        input_subject_idx_per_fiber.append(subj_idx)
        input_group_idx_per_fiber.append(group_indices[subj_idx])
g0 = numpy.nonzero(group_indices == 0)
g1 = numpy.nonzero(group_indices == 1)
    
[count_info, prob_info, KL, subject_count] = compute_local_stats(distances, distances_sq, neighborhood_sigma, group_indices, input_group_idx_per_fiber, input_subject_idx_per_fiber)

kl_total =  numpy.sum(KL)
kl_total_p = len(numpy.nonzero(null_kl_total <= kl_total)[0])/float(number_of_permutations)

statistic_mask = count_info[0] < statistic_threshold
island_p_c, corrected_p_count, island_labels_c, island_ids_c, island_sizes_c = island_correct(distances_centroids, distance_threshold, statistic_mask, null_count)

statistic_mask = prob_info[0] < statistic_threshold
island_p_p, corrected_p_prob, island_labels_p, island_ids_p, island_sizes_p = island_correct(distances_centroids, distance_threshold, statistic_mask, null_prob)

statistic_mask = KL > kl_threshold
island_p_kl, corrected_p_kl, island_labels_kl, island_ids_kl, island_sizes_kl = island_correct(distances_centroids, distance_threshold, statistic_mask, null_kl)

# outlier removal. 1% threshold gets rid of 7 noise fibers outside of brain. 1.1% 27 fibers
prob_fiber_thresh = 1.1
mask = (prob_info[1] > prob_fiber_thresh) & (prob_info[3] > prob_fiber_thresh)
mean_wm_clean = wma.filter.mask(mean_wm, mask, prob_info[1] + prob_info[3])

outdir_current = os.path.join(outdir, 'count_p_corrected')
visualize_data_array_to_disk(mean_wm_clean, corrected_p_count[mask], 'count_p_corrected', outdir_current, p_values=True, p_threshold=statistic_threshold)
outdir_current = os.path.join(outdir, 'prob_p_corrected')
visualize_data_array_to_disk(mean_wm_clean, corrected_p_prob[mask], 'prob_p_corrected', outdir_current, p_values=True, p_threshold=statistic_threshold)
outdir_current = os.path.join(outdir, 'kl_p_corrected')
visualize_data_array_to_disk(mean_wm_clean, corrected_p_kl[mask], 'kl_p_corrected', outdir_current, p_values=True, p_threshold=statistic_threshold)

outfile = os.path.join(outdir, "output.txt")
text_file = open(outfile, "w")
text_file.write("COUNT p: {0}\n".format(island_p_c))
text_file.write("PROB p: {0}\n".format(island_p_p))
text_file.write("KL p: {0}\n".format(island_p_kl))
text_file.write("KL total val: {0}\n".format(kl_total))
text_file.write("KL total p: {0}\n".format(kl_total_p))
text_file.close()

outdir_current = os.path.join(outdir, 'k_l')
visualize_data_array_to_disk(mean_wm_clean, KL[mask], 'k_l', outdir_current, scalar_range=[0,numpy.max(KL)])

outdir_current = os.path.join(outdir, 'count_p')
visualize_data_array_to_disk(mean_wm_clean, count_info[0][mask], 'count_p', outdir_current, scalar_range=[0,1])
outdir_current = os.path.join(outdir, 'count_mean_g0')
visualize_data_array_to_disk(mean_wm_clean, count_info[1][mask], 'count_mean_g0', outdir_current)
outdir_current = os.path.join(outdir, 'count_std_g0')
visualize_data_array_to_disk(mean_wm_clean, count_info[2][mask], 'count_std_g0', outdir_current)
outdir_current = os.path.join(outdir, 'count_perc_err_g0')
visualize_data_array_to_disk(mean_wm_clean, numpy.divide(count_info[2], count_info[1])[mask], 'count_perc_err_g0', outdir_current)
outdir_current = os.path.join(outdir, 'count_mean_g1')
visualize_data_array_to_disk(mean_wm_clean, count_info[3][mask], 'count_mean_g1', outdir_current)
outdir_current = os.path.join(outdir, 'count_std_g1')
visualize_data_array_to_disk(mean_wm_clean, count_info[4][mask], 'count_std_g1', outdir_current)
outdir_current = os.path.join(outdir, 'count_perc_err_g1')
visualize_data_array_to_disk(mean_wm_clean, numpy.divide(count_info[4], count_info[3])[mask], 'count_perc_err_g1', outdir_current)


outdir_current = os.path.join(outdir, 'prob_p')
visualize_data_array_to_disk(mean_wm_clean, prob_info[0][mask], 'prob_p', outdir_current, scalar_range=[0,1])
outdir_current = os.path.join(outdir, 'prob_mean_g0')
visualize_data_array_to_disk(mean_wm_clean, prob_info[1][mask], 'prob_mean_g0', outdir_current)
outdir_current = os.path.join(outdir, 'prob_std_g0')
visualize_data_array_to_disk(mean_wm_clean, prob_info[2][mask], 'prob_std_g0', outdir_current)
outdir_current = os.path.join(outdir, 'prob_perc_err_g0')
visualize_data_array_to_disk(mean_wm_clean, numpy.divide(prob_info[2], prob_info[1])[mask], 'prob_perc_err_g0', outdir_current)
outdir_current = os.path.join(outdir, 'prob_mean_g1')
visualize_data_array_to_disk(mean_wm_clean, prob_info[3][mask], 'prob_mean_g1', outdir_current)
outdir_current = os.path.join(outdir, 'prob_std_g1')
visualize_data_array_to_disk(mean_wm_clean, prob_info[4][mask], 'prob_std_g1', outdir_current)
outdir_current = os.path.join(outdir, 'prob_perc_err_g1')
visualize_data_array_to_disk(mean_wm_clean, numpy.divide(prob_info[4], prob_info[3])[mask], 'prob_perc_err_g1', outdir_current)


try:
    from joblib import Parallel, delayed
    USE_PARALLEL = 1
except ImportError:
    USE_PARALLEL = 0
    print "<compute_statistics.py> Failed to import joblib, cannot multiprocess."
    print " Please install joblib for this functionality."


execfile('/Users/odonnell/Dropbox/Coding/Python/whitematteranalysis/bin/publications/ipmi-stats-paper/stats-computation-code.py')

# ================= PARAMETER SETTINGS =================

#exp1 = False
#exp2 = False
#exp3 = False

# experiment 1
# ==============
if exp1:
    outdir = 'sig25_dthresh10'
    # spatial scale
    neighborhood_sigma = 25.0
    # spatial thresholds for cluster computation
    distance_threshold = 10.0
    compute_mean_wm = True
    # 20 choose 10 is: 184756
    number_of_permutations = 1000
    #number_of_permutations = 10

# experiment 2
# ==============
if exp2:
    outdir = 'sig15_dthresh10'
    # spatial scale
    neighborhood_sigma = 15.0
    # spatial thresholds for cluster computation
    distance_threshold = 10.0
    compute_mean_wm = True
    # 20 choose 10 is: 184756
    number_of_permutations = 1000
    #number_of_permutations = 10    

# experiment 3
# ==============
if exp3:
    outdir = 'sig20_dthresh10'
    # spatial scale
    neighborhood_sigma = 20.0
    # spatial thresholds for cluster computation
    distance_threshold = 10.0
    compute_mean_wm = True
    # 20 choose 10 is: 184756
    number_of_permutations = 1000
    #number_of_permutations = 10
    
#outdir = 'test_sig15_dthresh15'
# spatial scale
#neighborhood_sigma = 15.0
# spatial thresholds for cluster computation
#distance_threshold = 15.0

# thresholds for cluster computation
statistic_threshold = 0.05
#statistic_threshold = 0.1
#kl_threshold = 0.001
kl_threshold = 0.02
#kl_threshold = 0.0004
# what should this be??
bc_threshold = 0.1
diff2_threshold = 0.4
diff_threshold = 0.01

# final thresholds for corrected p values
final_statistic_threshold = 0.05
n_jobs = 8

# ================= END OF SETTINGS =================
# define dataset (gender data N=20, 10+10)
group_indices_real = numpy.array([0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1])

# compute mean brain fibers
# ================= 
print "Computing mean white matter"
if compute_mean_wm:
    mean_wm, mean_wm_weights = compute_avg_brain(fiber_array, distances, neighborhood_sigma)

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
null_bc = numpy.zeros(number_of_permutations)
null_bc_total = numpy.zeros(number_of_permutations)
null_diff2 = numpy.zeros(number_of_permutations)
null_diff2_total = numpy.zeros(number_of_permutations)
null_diff = numpy.zeros(number_of_permutations)
null_diff_total = numpy.zeros(number_of_permutations)

null_diff2, null_diff2_total, null_diff, null_diff_total

# parallel computation of the null distribution
USE_PARALLEL = 0
if USE_PARALLEL:
    print "Using parallel computation DON'T THIS IS SLOWER THAN REGULAR! BUG/INEFFICIENT"
    retval = Parallel(n_jobs=n_jobs, verbose=1)(
        delayed(do_permutation)(
            pidx, number_of_permutations, group_indices_real, number_of_subjects, number_of_fibers_per_subject, distances, distances_sq, neighborhood_sigma, statistic_threshold, distance_threshold, kl_threshold, bc_threshold
            ) for pidx in range(number_of_permutations))

    print "Done parallel"
    retval = numpy.array(retval)
    null_count = retval[:,0]
    null_prob = retval[:,1]
    null_kl = retval[:,2]
    null_kl_total = retval[:,3]
    null_bc = retval[:,4]
    null_bc_total = retval[:,5]

    null_diff2 = retval[:,6]
    null_diff2_total = retval[:,7]

    null_diff = retval[:,8]
    null_diff_total = retval[:,9]

else:
    for pidx in range(number_of_permutations):
        [null_count[pidx], null_prob[pidx], null_kl[pidx], null_kl_total[pidx], null_bc[pidx], null_bc_total[pidx], null_diff2[pidx], null_diff2_total[pidx], null_diff[pidx], null_diff_total[pidx]] = do_permutation(pidx, number_of_permutations, group_indices_real, number_of_subjects, number_of_fibers_per_subject, distances, distances_sq, neighborhood_sigma, statistic_threshold, distance_threshold, kl_threshold, bc_threshold)

print "Done computing null dists"

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
    
[count_info, prob_info, KL, BC, DIFF2, DIFF, subject_count] = compute_local_stats(distances, distances_sq, neighborhood_sigma, group_indices, input_group_idx_per_fiber, input_subject_idx_per_fiber)

kl_total =  numpy.sum(KL)
kl_total_p = len(numpy.nonzero(null_kl_total <= kl_total)[0])/float(number_of_permutations)

bc_total =  numpy.sum(BC)
bc_total_p = len(numpy.nonzero(null_bc_total <= bc_total)[0])/float(number_of_permutations)

diff2_total =  numpy.sum(DIFF2)
diff2_total_p = len(numpy.nonzero(null_diff2_total <= diff2_total)[0])/float(number_of_permutations)

diff_total =  numpy.sum(DIFF)
diff_total_p = len(numpy.nonzero(null_diff_total <= diff_total)[0])/float(number_of_permutations)

statistic_mask = count_info[0] < statistic_threshold
island_p_c, corrected_p_count, island_labels_c, island_ids_c, island_sizes_c = island_correct(distances_centroids, distance_threshold, statistic_mask, null_count)

statistic_mask = prob_info[0] < statistic_threshold
island_p_p, corrected_p_prob, island_labels_p, island_ids_p, island_sizes_p = island_correct(distances_centroids, distance_threshold, statistic_mask, null_prob)

statistic_mask = KL > kl_threshold
island_p_kl, corrected_p_kl, island_labels_kl, island_ids_kl, island_sizes_kl = island_correct(distances_centroids, distance_threshold, statistic_mask, null_kl)

#statistic_mask = BC < bc_threshold
#island_p_bc, corrected_p_bc, island_labels_bc, island_ids_bc, island_sizes_bc = island_correct(distances_centroids, distance_threshold, statistic_mask, null_bc)

statistic_mask = DIFF2 > diff2_threshold
island_p_diff2, corrected_p_diff2, island_labels_diff2, island_ids_diff2, island_sizes_diff2 = island_correct(distances_centroids, distance_threshold, statistic_mask, null_diff2)

statistic_mask = DIFF > diff_threshold
island_p_diff, corrected_p_diff, island_labels_diff, island_ids_diff, island_sizes_diff = island_correct(distances_centroids, distance_threshold, statistic_mask, null_diff)

# outlier removal. 1% threshold gets rid of 7 noise fibers outside of brain. 1.1% 27 fibers
#prob_fiber_thresh = .005 * numpy.max(prob_info[1])
prob_fiber_thresh = 0.075
mask = (prob_info[1] > prob_fiber_thresh) & (prob_info[3] > prob_fiber_thresh)
mean_wm_clean = wma.filter.mask(mean_wm, mask, prob_info[1] + prob_info[3])

outdir_current = os.path.join(outdir, 'count_p_corrected')
visualize_data_array_to_disk(mean_wm_clean, corrected_p_count[mask], 'count_p_corrected', outdir_current, p_values=True, p_threshold=final_statistic_threshold)
outdir_current = os.path.join(outdir, 'prob_p_corrected')
visualize_data_array_to_disk(mean_wm_clean, corrected_p_prob[mask], 'prob_p_corrected', outdir_current, p_values=True, p_threshold=final_statistic_threshold)
outdir_current = os.path.join(outdir, 'kl_p_corrected')
visualize_data_array_to_disk(mean_wm_clean, corrected_p_kl[mask], 'kl_p_corrected', outdir_current, p_values=True, p_threshold=final_statistic_threshold)
#outdir_current = os.path.join(outdir, 'bc_p_corrected')
#visualize_data_array_to_disk(mean_wm_clean, corrected_p_bc[mask], 'bc_p_corrected', outdir_current, p_values=True, p_threshold=final_statistic_threshold)
outdir_current = os.path.join(outdir, 'diff_p_corrected')
visualize_data_array_to_disk(mean_wm_clean, corrected_p_diff[mask], 'diff_p_corrected', outdir_current, p_values=True, p_threshold=final_statistic_threshold)
outdir_current = os.path.join(outdir, 'diff_50_p_corrected')
visualize_data_array_to_disk(mean_wm_clean, corrected_p_diff2[mask], 'diff_50_corrected', outdir_current, p_values=True, p_threshold=final_statistic_threshold)

outfile = os.path.join(outdir, "output.txt")
text_file = open(outfile, "w")
text_file.write("COUNT p: {0}\n".format(island_p_c))
text_file.write("PROB p: {0}\n".format(island_p_p))
text_file.write("KL p: {0}\n".format(island_p_kl))
text_file.write("KL total val: {0}\n".format(kl_total))
text_file.write("KL total p: {0}\n".format(kl_total_p))
#text_file.write("BC p: {0}\n".format(island_p_bc))
text_file.write("BC total val: {0}\n".format(bc_total))
text_file.write("BC total p: {0}\n".format(bc_total_p))
text_file.write("DIFF2 (50) p: {0}\n".format(island_p_diff2))
text_file.write("DIFF2 total val: {0}\n".format(diff2_total))
text_file.write("DIFF2 total p: {0}\n".format(diff2_total_p))
text_file.write("DIFF p: {0}\n".format(island_p_diff))
text_file.write("DIFF total val: {0}\n".format(diff_total))
text_file.write("DIFF total p: {0}\n".format(diff_total_p))
text_file.close()

outdir_current = os.path.join(outdir, 'k_l')
visualize_data_array_to_disk(mean_wm_clean, KL[mask], 'k_l', outdir_current, scalar_range=[0,numpy.max(KL[mask])])

outdir_current = os.path.join(outdir, 'b_c')
visualize_data_array_to_disk(mean_wm_clean, BC[mask], 'b_c', outdir_current, scalar_range=[0,numpy.max(BC[mask])])

outdir_current = os.path.join(outdir, 'diff')
visualize_data_array_to_disk(mean_wm_clean, DIFF[mask], 'diff', outdir_current, scalar_range=[0,numpy.max(DIFF[mask])])
outdir_current = os.path.join(outdir, 'diff_50')
visualize_data_array_to_disk(mean_wm_clean, DIFF2[mask], 'diff2', outdir_current, scalar_range=[0,numpy.max(DIFF2[mask])])

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

outdir_current = os.path.join(outdir, 'differences')
visualize_data_array_to_disk(mean_wm_clean, prob_info[3][mask] - prob_info[1][mask], 'mean_g1_minus_g0', outdir_current)

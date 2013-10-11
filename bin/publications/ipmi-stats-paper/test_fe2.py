#import glob
#import os
#import numpy
#import scipy.stats
#import vtk
#import whitematteranalysis as wma

execfile('/Users/odonnell/Dropbox/Coding/Python/whitematteranalysis/bin/publications/ipmi-stats-paper/stats-computation-code.py')
# ================= SETTINGS HERE  =================
indir = '/Volumes/brage-raid/odonnell/TBI_FE_CONTROLS_PNL/Tracts/reg-fe-sz/iteration_5'
outdir = '.'
#number_of_fibers_per_subject = 300
#number_of_fibers_per_subject = 600
number_of_fibers_per_subject = 3000
#number_of_fiber_centroids = 100
#number_of_fiber_centroids = 500
number_of_fiber_centroids = 1500
points_per_fiber = 30
neighborhood_sigma = 20.0
#neighborhood_sigma = 10.0

minimum_fiber_length = 30
group_indices_real = numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
number_of_permutations = 1000
# NOTE ONLY 1024 possible, could enumerate them exactly...
#number_of_permutations = 600
#number_of_permutations = 200
#number_of_permutations = 50
# parameters for island computation
statistic_threshold = 0.05
kl_threshold = 10.0
distance_threshold = 10.0
# ================= END OF SETTINGS =================

# ================= read data
input_mask = "{0}/*.vtk".format(indir)
input_poly_datas = glob.glob(input_mask)
#input_poly_datas = input_poly_datas[0:3]
number_of_subjects = len(input_poly_datas)
print input_poly_datas
input_pds = list()
for fname in input_poly_datas:
    print fname
    pd = wma.io.read_polydata(fname)
    pd = wma.filter.preprocess(pd, minimum_fiber_length) 
    input_pds.append(pd)
    #print pd
# get requested number of lines per subject
input_pds_downsampled = list()
for pd in input_pds:
    input_pds_downsampled.append(wma.filter.downsample(pd, number_of_fibers_per_subject))

   
# ================= compute distances =================
# convert to array representation
print 'Appending inputs into one polydata'
appender = vtk.vtkAppendPolyData()
for pd in input_pds_downsampled:
    appender.AddInputData(pd)
appender.Update()
print 'Converting fibers to array representation for dist and averaging'
fiber_array = wma.fibers.FiberArray()
fiber_array.convert_from_polydata(appender.GetOutput(), points_per_fiber)
print 'Done converting fibers to array representation for dist and averaging'

# random sample of "center/centroid" fibers for stats computation
total_number_of_fibers = number_of_fibers_per_subject*number_of_subjects
fiber_sample = numpy.random.permutation(total_number_of_fibers - 1)
fiber_sample = fiber_sample[0:number_of_fiber_centroids]

# compute distances between sample and all other fibers
print "Computing distances"
distances = numpy.zeros([number_of_fiber_centroids, total_number_of_fibers])

for idx in range(number_of_fiber_centroids):
    if idx % 100 == 0:
            print 'distances: ', idx, '/', number_of_fiber_centroids
    fiber = fiber_array.get_fiber(fiber_sample[idx])
    distances[idx,:] = wma.similarity.fiber_distance(fiber, fiber_array, threshold=0, distance_method='Hausdorff')

distances_centroids = distances[:,fiber_sample]
distances_sq = numpy.multiply(distances, distances)
print "Done computing distances"

# compute mean brain fibers
# ================= 
print "Computing mean white matter"
mean_wm, mean_wm_weights = compute_avg_brain(fiber_array, distances, neighborhood_sigma)

# ================= neighborhood parameters
sigma_sq = neighborhood_sigma*neighborhood_sigma

# ================= 
# COMPUTE NULL DISTRIBUTION
# =================
print "Computing null distribution"
null_count = numpy.zeros([number_of_permutations, 1])
null_prob = numpy.zeros([number_of_permutations, 1])
null_kl = numpy.zeros([number_of_permutations, 1])

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

    [count_info, prob_info, KL] = compute_local_stats(distances, distances_sq, neighborhood_sigma, group_indices, input_group_idx_per_fiber, input_subject_idx_per_fiber)

    statistic_mask = count_info[0] < statistic_threshold
    null_count[pidx] = suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic_mask)
    statistic_mask = prob_info[0] < statistic_threshold
    null_prob[pidx] = suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic_mask)
    statistic_mask = KL > kl_threshold
    null_kl[pidx] = suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic_mask)

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

[count_info, prob_info, KL] = compute_local_stats(distances, distances_sq, neighborhood_sigma, group_indices, input_group_idx_per_fiber, input_subject_idx_per_fiber)

statistic_mask = count_info[0] < statistic_threshold
island_p_c, corrected_p_count, island_labels_c, island_ids_c, island_sizes_c = island_correct(distances_centroids, distance_threshold, statistic_mask, null_count)

statistic_mask = prob_info[0] < statistic_threshold
island_p_p, corrected_p_prob, island_labels_p, island_ids_p, island_sizes_p = island_correct(distances_centroids, distance_threshold, statistic_mask, null_prob)

statistic_mask = KL > kl_threshold
island_p_kl, corrected_p_kl, island_labels_kl, island_ids_kl, island_sizes_kl = island_correct(distances_centroids, distance_threshold, statistic_mask, null_kl)

outdir_current = os.path.join(outdir, 'count_p_corrected')
visualize_data_array_to_disk(mean_wm, corrected_p_count, 'count_p_corrected', outdir_current, p_values=True, p_threshold=statistic_threshold)
outdir_current = os.path.join(outdir, 'prob_p_corrected')
visualize_data_array_to_disk(mean_wm, corrected_p_prob, 'prob_p_corrected', outdir_current, p_values=True, p_threshold=statistic_threshold)
outdir_current = os.path.join(outdir, 'kl_p_corrected')
visualize_data_array_to_disk(mean_wm, corrected_p_kl, 'kl_p_corrected', outdir_current, p_values=True, p_threshold=statistic_threshold)

outfile = os.path.join(outdir, "output.txt")
text_file = open(outfile, "w")
text_file.write("COUNT p: {0}\n".format(island_p_c))
text_file.write("PROB p: {0}\n".format(island_p_p))
text_file.write("KL p: {0}\n".format(island_p_kl))
text_file.close()

outdir_current = os.path.join(outdir, 'k_l')
visualize_data_array_to_disk(mean_wm, KL, 'k_l', outdir_current, scalar_range=[0,30])

outdir_current = os.path.join(outdir, 'count_p')
visualize_data_array_to_disk(mean_wm, count_info[0], 'count_p', outdir_current, scalar_range=[0,1])
outdir_current = os.path.join(outdir, 'count_mean_g0')
visualize_data_array_to_disk(mean_wm, count_info[1], 'count_mean_g0', outdir_current)
outdir_current = os.path.join(outdir, 'count_std_g0')
visualize_data_array_to_disk(mean_wm, count_info[2], 'count_std_g0', outdir_current)
outdir_current = os.path.join(outdir, 'count_perc_err_g0')
visualize_data_array_to_disk(mean_wm, numpy.divide(count_info[2], count_info[1]), 'count_perc_err_g0', outdir_current)
outdir_current = os.path.join(outdir, 'count_mean_g1')
visualize_data_array_to_disk(mean_wm, count_info[3], 'count_mean_g1', outdir_current)
outdir_current = os.path.join(outdir, 'count_std_g1')
visualize_data_array_to_disk(mean_wm, count_info[4], 'count_std_g1', outdir_current)
outdir_current = os.path.join(outdir, 'count_perc_err_g1')
visualize_data_array_to_disk(mean_wm, numpy.divide(count_info[4], count_info[3]), 'count_perc_err_g1', outdir_current)


outdir_current = os.path.join(outdir, 'prob_p')
visualize_data_array_to_disk(mean_wm, prob_info[0], 'prob_p', outdir_current, scalar_range=[0,1])
outdir_current = os.path.join(outdir, 'prob_mean_g0')
visualize_data_array_to_disk(mean_wm, prob_info[1], 'prob_mean_g0', outdir_current)
outdir_current = os.path.join(outdir, 'prob_std_g0')
visualize_data_array_to_disk(mean_wm, prob_info[2], 'prob_std_g0', outdir_current)
outdir_current = os.path.join(outdir, 'prob_perc_err_g0')
visualize_data_array_to_disk(mean_wm, numpy.divide(prob_info[2], prob_info[1]), 'prob_perc_err_g0', outdir_current)
outdir_current = os.path.join(outdir, 'prob_mean_g1')
visualize_data_array_to_disk(mean_wm, prob_info[3], 'prob_mean_g1', outdir_current)
outdir_current = os.path.join(outdir, 'prob_std_g1')
visualize_data_array_to_disk(mean_wm, prob_info[4], 'prob_std_g1', outdir_current)
outdir_current = os.path.join(outdir, 'prob_perc_err_g1')
visualize_data_array_to_disk(mean_wm, numpy.divide(prob_info[4], prob_info[3]), 'prob_perc_err_g1', outdir_current)


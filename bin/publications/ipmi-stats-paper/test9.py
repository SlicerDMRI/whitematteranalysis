import glob
import os
import numpy
import scipy.stats
import vtk
import statsmodels.sandbox.stats.multicomp as mcc
import whitematteranalysis as wma

# ================= SETTINGS HERE  =================
indir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/Dec2012/ipmi-reg-gender/iteration_6'
outdir = '.'
#number_of_fibers_per_subject = 300
#number_of_fibers_per_subject = 600
number_of_fibers_per_subject = 3000
#number_of_fiber_centroids = 100
#number_of_fiber_centroids = 500
number_of_fiber_centroids = 1500
points_per_fiber = 30
neighborhood_sigma = 20.0

minimum_fiber_length = 30
#group_indices = numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
group_indices_real = numpy.array([0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1])
#g0 = numpy.nonzero(group_indices == 0)
#g1 = numpy.nonzero(group_indices == 1)
#number_of_permutations = 1000
# NOTE ONLY 1024 possible, could enumerate them exactly...
number_of_permutations = 600
# parameters for island computation
statistic_threshold = 0.05
kl_threshold = 10.0
distance_threshold = 10.0
# ================= END OF SETTINGS =================

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

def find_islands_and_sizes(distances, distance_threshold, statistic, threshold):
    # input distances between sample and all other fibers
    [number_of_fiber_centroids, total_number_of_fibers] = distances.shape
    # initialize data structures
    island = numpy.zeros([number_of_fiber_centroids,1]) - 100
    curr_island = 0
    # loop over centroids, agglomerating nearest neighborhoods
    for cidx in range(number_of_fiber_centroids):
        #print cidx
        if statistic[cidx] < threshold:
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
                # if this one can be agglomerated with the current island do it
                if statistic[cidx2] < threshold:
                    # if this one is unassigned so far assign to the current island
                    if island[cidx2] == -100:
                        #print "join: ", island[cidx2], island[cidx], cidx, cidx2
                        island[cidx2] = island[cidx]
                    elif island[cidx2] != island[cidx]:
                        # otherwise we have to agglomerate the existing island with this one.
                        merge_idx = island[cidx2]
                        #print "agglomerate: ", merge_idx, island[cidx], cidx, cidx2
                        island[numpy.nonzero(island == merge_idx)] = island[cidx]
                else:
                    # both statistic and distance must be over threshold to stop agglomerating
                    if distances[cidx, cidx2] > distance_threshold:
                        growing = 0
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
    
    # compute stats
    fiber_count_p = numpy.ones([number_of_fiber_centroids,1])
    fiber_count_mean_g0 = numpy.zeros([number_of_fiber_centroids,1])
    fiber_count_std_g0 = numpy.zeros([number_of_fiber_centroids,1])
    fiber_count_mean_g1 = numpy.zeros([number_of_fiber_centroids,1])
    fiber_count_std_g1 = numpy.zeros([number_of_fiber_centroids,1])
    
    fiber_probability_p  = numpy.ones([number_of_fiber_centroids,1])
    fiber_probability_mean_g0 = numpy.zeros([number_of_fiber_centroids,1])
    fiber_probability_std_g0 = numpy.zeros([number_of_fiber_centroids,1])
    fiber_probability_mean_g1 = numpy.zeros([number_of_fiber_centroids,1])
    fiber_probability_std_g1 = numpy.zeros([number_of_fiber_centroids,1])

    KL = numpy.ones([number_of_fiber_centroids,1])
     
    for idx in range(number_of_fiber_centroids):
        if idx % 100 == 0:
            print idx, '/', number_of_fiber_centroids
        neighborhood_indices = numpy.nonzero(distances[idx,:] < neighborhood_threshold)[0]
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
        # calculate soft fiber counts (probability density at this centroid)
        fiber_probabilities = numpy.exp(numpy.divide(-d_sq, sigma_sq))
        subject_probabilities = numpy.zeros(number_of_subjects)
        for sidx in range(number_of_subjects):
            # sum fiber probabilities from this subject
            sfibers = numpy.nonzero(subject_list == sidx)[0]
            subject_probabilities[sidx] = numpy.sum(fiber_probabilities[sfibers])
        [t, p] = scipy.stats.ttest_ind(subject_probabilities[g0], subject_probabilities[g1])
        fiber_probability_p[idx] = p
        fiber_probability_mean_g0[idx] = numpy.mean(subject_probabilities[g0])
        fiber_probability_std_g0[idx] =  numpy.std(subject_probabilities[g0])
        fiber_probability_mean_g1[idx] = numpy.mean(subject_probabilities[g1])
        fiber_probability_std_g1[idx] = numpy.std(subject_probabilities[g1])
        fiber_probability_retval = [fiber_probability_p, fiber_probability_mean_g0, fiber_probability_std_g0, fiber_probability_mean_g1, fiber_probability_std_g1]

        # now KL divergence: sum probability in each group and compute: (p-q)(lnp-lnq)
        eps = 0.000000001
        p0 = numpy.sum(subject_probabilities[g0])
        p1 = numpy.sum(subject_probabilities[g1])
        KL[idx] = numpy.multiply((p0 - p1 + eps),(numpy.log(p0 + eps) - numpy.log(p1 + eps)))
        
        
    return fiber_count_retval, fiber_probability_retval, KL


def suprathreshold_cluster_max(distances_centroids, distance_threshold, statistic, statistic_threshold):
        # now find max suprathreshold cluster size
    island_labels, island_ids, island_sizes = find_islands_and_sizes(distances_centroids, distance_threshold, statistic, statistic_threshold)
    non_background_island_sizes = island_sizes[numpy.nonzero(island_ids != -100)[0]]
    # if there were any significant islands append the maximum size
    if len(non_background_island_sizes) > 0:
        return numpy.max(non_background_island_sizes)
    else:
        # otherwise append 0 for no significant islands
        return 0

def island_correct(distances_centroids, distance_threshold, statistic, statistic_threshold, null):
    # now find max suprathreshold cluster size
    island_labels, island_ids, island_sizes = find_islands_and_sizes(distances_centroids, distance_threshold, statistic, statistic_threshold)
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

def visualize_data_array_to_disk(polydata, data_array, array_name, outdir, scalar_range=None):
    add_array_to_polydata(polydata, data_array, array_name=array_name)
    mean_wm.GetCellData().SetActiveScalars(array_name)
    if not os.path.exists(outdir):
            os.makedirs(outdir)
    if scalar_range == None:
        ren = wma.render.render(polydata, scalar_bar=True)
    else:
        ren = wma.render.render(polydata, scalar_bar=True, scalar_range=scalar_range)
    ren.save_views(outdir)
    del ren
    
# ================= END FUNCTIONS (move to module code) =======================

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
group_indices_to_permute = numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
null_count = numpy.zeros([number_of_permutations, 1])
null_prob = numpy.zeros([number_of_permutations, 1])
null_kl = numpy.zeros([number_of_permutations, 1])

for pidx in range(number_of_permutations):
    print "==========================================================================="
    print "permutation: ", pidx, "/", number_of_permutations
    print "==========================================================================="    
    # assign groups for each subject and fiber
    group_indices = numpy.random.permutation(group_indices_to_permute)
    input_subject_idx_per_fiber = list()
    input_group_idx_per_fiber = list()
    for subj_idx in range(number_of_subjects):
        for fiber_idx in range(number_of_fibers_per_subject):
            input_subject_idx_per_fiber.append(subj_idx)
            input_group_idx_per_fiber.append(group_indices[subj_idx])
    g0 = numpy.nonzero(group_indices == 0)
    g1 = numpy.nonzero(group_indices == 1)

    [count_info, prob_info, KL] = compute_local_stats(distances, distances_sq, neighborhood_sigma, group_indices, input_group_idx_per_fiber, input_subject_idx_per_fiber)

    null_count[pidx] = suprathreshold_cluster_max(distances_centroids, distance_threshold, count_info[0], statistic_threshold)
    null_prob[pidx] = suprathreshold_cluster_max(distances_centroids, distance_threshold, prob_info[0], statistic_threshold)
    null_kl[pidx] = suprathreshold_cluster_max(distances_centroids, distance_threshold, KL, kl_threshold)

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

island_p_c, corrected_p_count, island_labels_c, island_ids_c, island_sizes_c = island_correct(distances_centroids, distance_threshold, count_info[0], statistic_threshold, null_count)

island_p_p, corrected_p_prob, island_labels_p, island_ids_p, island_sizes_p = island_correct(distances_centroids, distance_threshold, prob_info[0], statistic_threshold, null_prob)

island_p_kl, corrected_p_kl, island_labels_kl, island_ids_kl, island_sizes_kl = island_correct(distances_centroids, distance_threshold, KL, kl_threshold, null_kl)

outdir_current = os.path.join(outdir, 'count_p_corrected')
visualize_data_array_to_disk(mean_wm, corrected_p_count, 'count_p_corrected', outdir_current, scalar_range=[0,0.05])
outdir_current = os.path.join(outdir, 'prob_p_corrected')
visualize_data_array_to_disk(mean_wm, corrected_p_prob, 'prob_p_corrected', outdir_current, scalar_range=[0,0.05])
outdir_current = os.path.join(outdir, 'kl_p_corrected')
visualize_data_array_to_disk(mean_wm, corrected_p_kl, 'kl_p_corrected', outdir_current, scalar_range=[0,0.05])

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


add_array_to_polydata(mean_wm, island_labels_c, array_name='Islands_count')
mean_wm.GetCellData().SetActiveScalars('Islands_count')
ren2 = wma.render.render(mean_wm, scalar_bar=True, scalar_range=[-10,len(island_ids_c)])

# mask pd if outlier hood, less than half of subjects
#mask = numpy.array(number_subjects_list) > 9
#outpd_masked = wma.filter.mask(outpd, mask, numpy.array(p_num_fibers_list))
#ren2 = wma.render.render(outpd_masked, scalar_bar=True, scalar_range=[0,1])



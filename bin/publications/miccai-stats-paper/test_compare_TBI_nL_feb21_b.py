import os
import glob
import matplotlib.pyplot as plt
import numpy
import scipy.stats

import vtk

import whitematteranalysis as wma

import multiprocessing


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


parallel_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', parallel_jobs
#parallel_jobs *= 3
#parallel_jobs = 101
parallel_jobs = 15
#parallel_jobs = 10
print 'Using N jobs:', parallel_jobs

#group_indices = [1, 0, 1, 0, 0, 1, 1, 0]
# 1 T, 2 C, 3 T, 4 C, 5 C, 6 T, 7 T, 8 C

execfile('/Users/odonnell/Dropbox/Coding/Python/WhiteMatterAnalysis/bin/test_compute_FA.py')

indir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/MICCAI2012/tbi_with_scalars'

input_mask = "{0}/*.vtk".format(indir)
input_poly_datas = glob.glob(input_mask)

print input_poly_datas

input_pds = list()
input_pds_downsampled = list()

#number_of_fibers_per_subject = 3000
#number_of_fiber_centroids = 1000
# this is about 2.4 GB of memory for the distances...
number_of_fibers_per_subject = 6000
number_of_fiber_centroids = 2000
number_of_subjects = len(input_poly_datas)
points_per_fiber = 30

#input_poly_datas = input_poly_datas[0:1]

# this is to produce the files with scalars
if 0:
    for fname in input_poly_datas:
        print fname
        pd = wma.io.read_polydata(fname)
        pd, fa_lines_list, fa_avg_list = compute_scalar_measures(pd)
        fname2 =  'scalars_' + os.path.basename(fname)
        wma.io.write_polydata(pd, fname2)

# read in ones with scalars already
for fname in input_poly_datas:
    print fname
    pd = wma.io.read_polydata(fname)
    input_pds.append(pd)

# grab scalars of interest
input_mean_fas_per_subject = list()
input_pds_downsampled = list()
downsample_indices = list()
for pd in input_pds:
    pd2, fiber_indices = wma.filter.downsample(pd, number_of_fibers_per_subject,return_indices=True)
    # get FA only at fibers of interest
    pd.GetCellData().RemoveArray('mean_FA')
    # the files on disk only have median FA, get mean instead
    pd = compute_mean_measures(pd)
    fa = pd.GetCellData().GetArray('mean_FA')
    fa_subj = list()
    for idx in fiber_indices:
        fa_subj.append(fa.GetTuple1(idx))
    #fa_subj = numpy.array(fa_subj)
    #fa_subj = numpy.array(fa_avg_list)[fiber_indices]
    input_mean_fas_per_subject.append(fa_subj)    
    input_pds_downsampled.append(pd2)
    downsample_indices.append(fiber_indices)
    

# convert to arrays for dist and averaging
# use entire appended polydata (perhaps in future compute per-subject)
print 'Appending inputs into one polydata'
appender = vtk.vtkAppendPolyData()
for pd in input_pds_downsampled:
    appender.AddInput(pd)

appender.Update()
print 'Done appending inputs into one polydata'

# convert to array representation
print 'Converting fibers to array representation for dist and averaging'
fiber_array = wma.fibers.FiberArray()
fiber_array.convert_from_polydata(appender.GetOutput(), points_per_fiber)
print 'Done converting fibers to array representation for dist and averaging'

# try to do some statistics
# random sample of fibers for stats
total_number_of_fibers = number_of_fibers_per_subject*number_of_subjects
fiber_sample = numpy.random.permutation(total_number_of_fibers - 1)
fiber_sample = fiber_sample[0:number_of_fiber_centroids]

# compute dists
# find the sample's distances to all other fibers
distances = numpy.zeros([number_of_fiber_centroids, total_number_of_fibers])

for idx in range(number_of_fiber_centroids):
    print idx, '/', number_of_fiber_centroids
    fiber = fiber_array.get_fiber(fiber_sample[idx])
    distances[idx,:] = wma.similarity.fiber_distance(fiber, fiber_array, threshold=0, distance_method='Hausdorff')



# ------------------------------------------------------
# Can re-run the below after changing neighborhood_threshold
# Or after changing group membership
# All slow processing happens above this line
# ------------------------------------------------------
neighborhood_threshold = 20

group_indices =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,\
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

# test for artifacts
#group_indices =  [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,\
#                  1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
#group_indices =  [1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0,\
#                  1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0]
#group_indices = numpy.random.randint(0,2,26)
group_indices = numpy.array(group_indices)



# assign to flat lists with subject idx
subj_idx = 0
input_data_per_fiber = list()
input_subject_idx_per_fiber = list()
input_group_idx_per_fiber = list()
for fa_subj in input_mean_fas_per_subject:
    input_data_per_fiber += fa_subj
    for data_point in fa_subj:
        input_subject_idx_per_fiber.append(subj_idx)
        input_group_idx_per_fiber.append(group_indices[subj_idx])
    subj_idx +=1
    
# according to neighborhood definition:
# compute avg fibers and FA stats in neighborhoods
average_fibers_list = list()
statistic_list = list()
significance_list = list()
density_statistic_list = list()
density_significance_list = list()
density_0_list = list()
density_1_list = list()
data_0_list = list()
data_1_list = list()
data_0_std_list = list()
data_1_std_list = list()

for idx in range(number_of_fiber_centroids):
    print idx, '/', number_of_fiber_centroids
    neighborhood_indices = numpy.nonzero(distances[idx,:] < neighborhood_threshold)[0]
    hood_count = len(neighborhood_indices)
    # find average fiber in neighborhood
    avg_fiber = fiber_array.get_fiber(neighborhood_indices[0])
    for hood_idx in neighborhood_indices[1:]:
        avg_fiber += fiber_array.get_fiber(hood_idx)
    avg_fiber /= hood_count
    average_fibers_list.append(avg_fiber)
    # compute statistic(s) in neighborhood
    data_list = list()
    group_list = list()
    subject_list = list()
    for hood_idx in neighborhood_indices:
        data_list.append(input_data_per_fiber[hood_idx])
        group_list.append(input_group_idx_per_fiber[hood_idx])
        subject_list.append(input_subject_idx_per_fiber[hood_idx])
    # statistic: "fiber density" in group 1 vs group 2
    subject_stats_list = list()
    subject_list = numpy.array(subject_list)
    for sidx in range(number_of_subjects):
        subject_stats_list.append(len(numpy.nonzero(subject_list == sidx)[0]))
    subject_stats_list = numpy.array(subject_stats_list)
    # figure out which subject is from which group
    g0 = numpy.nonzero(group_indices == 0)[0]
    g1 = numpy.nonzero(group_indices == 1)[0]
    t, p = scipy.stats.ttest_ind(subject_stats_list[g0], subject_stats_list[g1])
    density_0_list.append(numpy.sum(subject_stats_list[g0]))
    density_1_list.append(numpy.sum(subject_stats_list[g1]))
    density_statistic_list.append(t)
    density_significance_list.append(p)    
    # statistic: FA in group 1 vs FA in group 2
    subject_stats_list = list()
    data_list = numpy.array(data_list)
    for sidx in range(number_of_subjects):
        subject_fibers = numpy.nonzero(subject_list == sidx)[0]
        if len(subject_fibers):
            subject_stats_list.append(numpy.mean(data_list[subject_fibers]))
        else:
            subject_stats_list.append(numpy.nan)
    subject_stats_list = numpy.array(subject_stats_list)
    g0 = numpy.nonzero((group_indices == 0) & ~numpy.isnan(subject_stats_list))[0]
    g1 = numpy.nonzero((group_indices == 1) & ~numpy.isnan(subject_stats_list))[0]
    if len(g0) and len(g1):
        # non parametric
        #t1=numpy.random.randint(1,10,10)
        #t2=numpy.random.randint(1,10,10)
        #t, p = scipy.stats.ks_2samp(t1,t2)
        #t, p = scipy.stats.ks_2samp(subject_stats_list[g0], subject_stats_list[g1])
        t, p = scipy.stats.ttest_ind(subject_stats_list[g0], subject_stats_list[g1])
        data_0_list.append(numpy.mean(subject_stats_list[g0]))
        data_1_list.append(numpy.mean(subject_stats_list[g1]))
        data_0_std_list.append(numpy.std(subject_stats_list[g0]))
        data_1_std_list.append(numpy.std(subject_stats_list[g1]))
    else:
        t = p = numpy.nan
        data_0_list.append(numpy.nan)
        data_1_list.append(numpy.nan)
        data_0_std_list.append(numpy.nan)
        data_1_std_list.append(numpy.nan)
    statistic_list.append(t)
    significance_list.append(p)
    
f = open('pvals_Density.txt','w')
for lidx in range(len(density_significance_list)):
    f.write(str(density_significance_list[lidx]))
    f.write('\n')
f.close()

f = open('pvals_FA.txt','w')
for lidx in range(len(significance_list)):
    f.write(str(significance_list[lidx]))
    f.write('\n')
f.close()

print 'load pvals.txt into matlab to test fdr for now'

plt.figure()
plt.plot(data_0_list, data_1_list, 'o')
plt.title('FA in neighborhoods')
plt.xlabel('Group 0')
plt.ylabel('Group 1')
plt.savefig('FA-groups-neighborhoods.pdf')
plt.close()

plt.figure()
plt.plot(density_0_list, density_1_list, 'o')
plt.title('Number of trajectories in neighborhoods')
plt.xlabel('Group 0')
plt.ylabel('Group 1')
plt.savefig('density-groups-neighborhoods.pdf')
plt.close()

plt.figure()
pvals = numpy.array(significance_list)
plt.hist(pvals[~numpy.isnan(pvals)],300)
plt.savefig('pvals-data-neighborhoods.pdf')
plt.close()

plt.figure()
svals = numpy.array(statistic_list)
plt.hist(svals[~numpy.isnan(svals)],300)
plt.savefig('statistic-data-neighborhoods.pdf')
plt.close()

plt.figure()
plt.plot(numpy.array(data_1_list) - numpy.array(data_0_list), pvals, 'o')
plt.savefig('difference_vs_p.pdf')
plt.close()

plt.figure()
plt.plot(numpy.array(data_1_std_list) - numpy.array(data_0_std_list), pvals, 'o')
plt.savefig('std_difference_vs_p.pdf')
plt.close()

plt.figure()
plt.plot(numpy.array(data_1_std_list), numpy.array(data_0_std_list), 'o')
plt.title('Data std in neighborhoods')
plt.xlabel('Group 0')
plt.ylabel('Group 1')
plt.savefig('std-groups-neighborhoods.pdf')
plt.close()

plt.figure()
plt.plot(numpy.array(data_1_std_list), pvals, 'o')
plt.title('Data std g1 vs pval')
plt.savefig('std-pval-neighborhoods.pdf')
plt.close()


plt.figure()
plt.plot(numpy.array(density_1_list), pvals, 'o')
plt.title('Data density g1 vs pval')
plt.savefig('density-pval-neighborhoods.pdf')
plt.close()

plt.figure()
test = numpy.array(data_1_list) - numpy.array(data_0_list)
mask = ~numpy.isnan(test)
#test = numpy.divide(test, numpy.array(data_1_std_list) + numpy.array(data_0_std_list))
plt.hist(test[mask], 800)
plt.title('Approx test in neighborhoods')
plt.savefig('test-groups-neighborhoods.pdf')
plt.close()

# output as pd
outpd = fiber_list_to_fiber_array(average_fibers_list).convert_to_polydata()
outpd = add_array_to_polydata(outpd, significance_list, array_name='P')
outpd = add_array_to_polydata(outpd, data_0_list, array_name='Data0')
outpd = add_array_to_polydata(outpd, data_1_list, array_name='Data1')
outpd = add_array_to_polydata(outpd, data_0_std_list, array_name='Data0Std')
outpd = add_array_to_polydata(outpd, data_1_std_list, array_name='Data1Std')
outpd = add_array_to_polydata(outpd, density_significance_list, array_name='P-density')
outpd = add_array_to_polydata(outpd, density_0_list, array_name='Density0')
outpd = add_array_to_polydata(outpd, density_1_list, array_name='Density1')
outpd = add_array_to_polydata(outpd, numpy.array(data_1_list)-numpy.array(data_0_list), array_name='Difference')
outpd.GetCellData().SetActiveScalars('Data0')

mask = ~numpy.isnan(data_1_list)
mask_idx = numpy.nonzero(mask)[0]
mask_pd = wma.filter.mask(outpd, mask)
mask_pd = add_array_to_polydata(mask_pd, numpy.array(data_1_list)[mask_idx], array_name='Data1')
mask_pd = add_array_to_polydata(mask_pd, numpy.array(data_0_list)[mask_idx], array_name='Data0')
mask_pd = add_array_to_polydata(mask_pd, numpy.array(data_1_list)[mask_idx]-numpy.array(data_0_list)[mask_idx], array_name='Difference')
mask_pd = add_array_to_polydata(mask_pd, numpy.array(significance_list)[mask_idx], array_name='P')
mask_pd.GetCellData().SetActiveScalars('Data0')
mask_pd.GetCellData().SetActiveScalars('Data1')
#mask_pd.GetCellData().SetActiveScalars('Difference')
#mask_pd.GetCellData().SetActiveScalars('P')
mask_pd = add_array_to_polydata(mask_pd, test[mask_idx], array_name='test')

ren = wma.render.render(outpd, scalar_bar=True, scalar_range=[0.35,0.65])

mask_pd.GetCellData().SetActiveScalars('Difference')
ren2 = wma.render.render(mask_pd, scalar_bar=True, scalar_range=[-.05,.05])

# from matlab fdr
# hood size 15 mm
# [pID,pN] = fdr(pval_tbi, .05)
# this is the pID (pN:2.5276e-04)
thresh = 0.0054
# hood size 20mm
# [pID,pN] = fdr(pval_tbi, .05)
# this is the pID (pN:4.4417e-04)
thresh = 0.0078
# for .1
thresh = 0.0243
# hood size 30 (too large, arcuates disappear)
thresh = 0.0414

# FA, nhood 15, fdr .1
thresh = 0.0153
# hood 20
thresh = 0.0278
#ren2 = wma.render.render(mask_pd, scalar_bar=True, scalar_range=[0,thresh])

# testing other groups, also ks
thresh = 0.001

mask = (numpy.array(significance_list) < thresh) & (numpy.array(density_significance_list) > 0.05)
mask_idx = numpy.nonzero(mask)[0]
mask_pd2 = wma.filter.mask(outpd, mask)
mask_pd2 = add_array_to_polydata(mask_pd2, numpy.array(significance_list)[mask_idx], array_name='P')
mask_pd2 = add_array_to_polydata(mask_pd2, numpy.array(data_1_list)[mask_idx]-numpy.array(data_0_list)[mask_idx], array_name='Difference')
#mask_pd2 = add_array_to_polydata(mask_pd2, test[mask_idx], array_name='test')

ren3 = wma.render.render(mask_pd2, scalar_bar=True)


# suprathreshold clusters??
# say p .1
# say want clusters of 5 or more fibers.
# must be "neighbors"
# do we have any of these?
centroid_dists = distances[:,fiber_sample]
cluster_neighbors_list = list()
cluster_size_list = list()
pval_threshold = 0.05
cluster_size_threshold = 5

pval_threshold = 0.1
pval_threshold = 0.15
pval_threshold = 0.2
cluster_size_threshold = 3

for idx in range(number_of_fiber_centroids):
    print idx, '/', number_of_fiber_centroids
    if pvals[idx] < pval_threshold:
        neighbors = numpy.nonzero(centroid_dists[idx] < neighborhood_threshold)[0]
        good_neighbors = numpy.nonzero(pvals[neighbors] < pval_threshold)[0]
        cluster_neighbors_list.append(neighbors[good_neighbors])
        cluster_size_list.append(len(good_neighbors))
        #print numpy.min(pvals[neighbors]), '...', numpy.max(pvals[neighbors])
        print len(good_neighbors)
    else:
        cluster_neighbors_list.append([])
        cluster_size_list.append(0)
        
# potentially significant clusters
cluster_size_list2 = numpy.array(cluster_size_list)

# include all fibers from each cluster in the size calculation
for idx in range(number_of_fiber_centroids):
    neighbors = cluster_neighbors_list[idx]
    for n_idx in neighbors:
        if cluster_size_list[n_idx] < cluster_size_list[idx]:
            cluster_size_list2[n_idx] = cluster_size_list[idx]
            
mask = cluster_size_list2 >= cluster_size_threshold

#mask_pd3 = wma.filter.mask(outpd, mask, numpy.array(significance_list))
mask_pd3 = wma.filter.mask(outpd, mask, pvals)
mask_pd3 = wma.filter.mask(outpd, mask, numpy.array(data_1_list)-numpy.array(data_0_list))
#mask_pd3 = add_array_to_polydata(mask_pd3, numpy.array(data_1_list)[mask_idx]-numpy.array(data_0_list)[mask_idx], array_name='Difference')
#mask_pd2 = add_array_to_polydata(mask_pd2, test[mask_idx], array_name='test')

ren4 = wma.render.render(mask_pd3, scalar_bar=True)

#======================
# test individuals
# =====
## NOTE: this needs to be re-done with fiber centroids
## from normal subjects only

age_control = [29, 43, 38, 31, 23, 29, 40, 24, 42, 26, 47, 23, 40]
age_tbi = [44, 37, 43, 27, 29, 42, 27, 24, 25, 29, 24, 39, 44]
age_control= numpy.array(age_control)
age_tbi= numpy.array(age_tbi)

# according to neighborhood definition:
# compute avg fibers and FA stats in neighborhoods
average_fibers_list = list()
statistic_list = list()
significance_list = list()
density_statistic_list = list()
density_significance_list = list()
density_0_list = list()
density_1_list = list()
data_0_list = list()
data_1_list = list()
data_0_std_list = list()
data_1_std_list = list()

for idx in range(number_of_fiber_centroids):
    print idx, '/', number_of_fiber_centroids
    neighborhood_indices = numpy.nonzero(distances[idx,:] < neighborhood_threshold)[0]
    hood_count = len(neighborhood_indices)
    # find average fiber in neighborhood
    avg_fiber = fiber_array.get_fiber(neighborhood_indices[0])
    for hood_idx in neighborhood_indices[1:]:
        avg_fiber += fiber_array.get_fiber(hood_idx)
    avg_fiber /= hood_count
    average_fibers_list.append(avg_fiber)
    # compute statistic(s) in neighborhood
    data_list = list()
    group_list = list()
    subject_list = list()
    for hood_idx in neighborhood_indices:
        data_list.append(input_data_per_fiber[hood_idx])
        group_list.append(input_group_idx_per_fiber[hood_idx])
        subject_list.append(input_subject_idx_per_fiber[hood_idx])
    # statistic: "fiber density" in group 1 vs group 2
    # fiber density in each subject in group 1, vs group 0 model
    subject_stats_list = list()
    subject_list = numpy.array(subject_list)
    for sidx in range(number_of_subjects):
        subject_stats_list.append(len(numpy.nonzero(subject_list == sidx)[0]))
    subject_stats_list = numpy.array(subject_stats_list)
    # figure out which subject is from which group
    g0 = numpy.nonzero(group_indices == 0)[0]
    g1 = numpy.nonzero(group_indices == 1)[0]
    #t, p = scipy.stats.ttest_ind(subject_stats_list[g0], subject_stats_list[g1])
    subject_density_significance_list = list()
    for sidx in g1:
        #t, p = scipy.stats.ttest_ind(subject_stats_list[g0], subject_stats_list[sidx])
        p = scipy.stats.zmap(subject_stats_list[sidx], subject_stats_list[g0])
        subject_density_significance_list.append(p)
        
    density_0_list.append(numpy.sum(subject_stats_list[g0]))
    density_1_list.append(numpy.sum(subject_stats_list[g1]))
    #density_statistic_list.append(t)
    density_significance_list.append(subject_density_significance_list)    
    # statistic: FA in group 1 vs FA in group 2
    subject_stats_list = list()
    data_list = numpy.array(data_list)
    for sidx in range(number_of_subjects):
        subject_fibers = numpy.nonzero(subject_list == sidx)[0]
        if len(subject_fibers):
            subject_stats_list.append(numpy.mean(data_list[subject_fibers]))
        else:
            subject_stats_list.append(numpy.nan)
    subject_stats_list = numpy.array(subject_stats_list)
    g0 = numpy.nonzero((group_indices == 0) & ~numpy.isnan(subject_stats_list))[0]
    #g1 = numpy.nonzero((group_indices == 1) & ~numpy.isnan(subject_stats_list))[0]
    g1 = numpy.nonzero(group_indices == 1)[0]
    subject_data_significance_list = list()
    if len(g0):
        for sidx in g1:
            slope, intercept, r, p, err = scipy.stats.linregress(subject_stats_list[g0], age_control[g0])
            if (p < 0.05):
                print "fitting line"
                subject_stats_list[g0] += -intercept - numpy.multiply(slope,age_control[g0])
                subject_stats_list[sidx] += -intercept - numpy.multiply(slope,age_tbi[sidx-13])
            p = scipy.stats.zmap(subject_stats_list[sidx], subject_stats_list[g0])
            #t, p = scipy.stats.ttest_ind(subject_stats_list[g0], subject_stats_list[sidx])
            subject_data_significance_list.append(p)
        #t, p = scipy.stats.ttest_ind(subject_stats_list[g0], subject_stats_list[g1])
        data_0_list.append(numpy.mean(subject_stats_list[g0]))
        data_1_list.append(numpy.mean(subject_stats_list[g1]))
        data_0_std_list.append(numpy.std(subject_stats_list[g0]))
        data_1_std_list.append(numpy.std(subject_stats_list[g1]))
    else:
        t = p = numpy.nan
        for sidx in g1:
            subject_data_significance_list.append(p)
        data_0_list.append(numpy.nan)
        data_1_list.append(numpy.nan)
        data_0_std_list.append(numpy.nan)
        data_1_std_list.append(numpy.nan)
    statistic_list.append(t)
    significance_list.append(subject_data_significance_list)
    

# view per-subject
g1 = numpy.nonzero(group_indices == 1)[0]
out_masked_pds_zscore_FA = list()
number_abnormal_high_subj = list()
number_abnormal_low_subj = list()
number_centroids_subj = list()
#abnormal_threshold = 3
abnormal_threshold = 2

for sidx in range(len(g1)):
    zscore_subj = list()
    for sigs in significance_list:
        zscore_subj.append(sigs[sidx])
    zscore_subj = numpy.array(zscore_subj)
    #mask = ~numpy.isnan(zscore_subj) & ~numpy.isinf(zscore_subj)
    number_abnormal_high_subj.append(numpy.sum(zscore_subj[mask] > abnormal_threshold))
    number_abnormal_low_subj.append(numpy.sum(zscore_subj[mask] < -abnormal_threshold))
    # test for just decreases
    #number_abnormal_subj.append(numpy.sum(zscore_subj[mask] < -abnormal_threshold))
    number_centroids_subj.append(numpy.sum(mask))
    out_masked_pds_zscore_FA.append(wma.filter.mask(outpd, mask, zscore_subj))

number_abnormal_high_subj = numpy.array(number_abnormal_high_subj)
number_abnormal_low_subj = numpy.array(number_abnormal_low_subj)
number_abnormal_subj = number_abnormal_high_subj + number_abnormal_low_subj

percent_abnormal_subj = numpy.divide(number_abnormal_subj.astype(float), numpy.array(number_centroids_subj))
percent_abnormal_high_subj = numpy.divide(number_abnormal_high_subj.astype(float), numpy.array(number_centroids_subj))
percent_abnormal_low_subj = numpy.divide(number_abnormal_low_subj.astype(float), numpy.array(number_centroids_subj))
# what percentage of zscores are more than 2 std away from mean?
# really should correct for age also!


print "LAUREN why are there -inf z scores??? are these bad samples from tbi only?"

#good 4 5
#moderate 2 3
#severe 0 1
#disability = [5, 5, 3, 2, 0, 5, 2, 5, 5, 4, 3, 3, 4]
# too many groups. combine into good high/low, moderate, severe
#disability = [5, 5, 2, 2, 0, 5, 2, 5, 5, 5, 2, 2, 5]

# last patient has no score yet.
indices = range(12)
disability = [5, 3, 2, 0, 5, 2, 5, 5, 4, 3, 3, 4]
stroop = numpy.array([112, 76, 90, 56, 112, 103, 112, 112, 111, 112, 99])
plt.figure()
plt.plot(age_tbi[indices], percent_abnormal_subj[indices], 'ko')
plt.plot(age_tbi[indices], percent_abnormal_high_subj[indices], 'ro')
plt.plot(age_tbi[indices], percent_abnormal_low_subj[indices], 'bo')
plt.savefig('age_vs_perc_abnormal.pdf')
print percent_abnormal_subj

plt.figure()
plt.plot(stroop[0:10], percent_abnormal_subj[0:10], 'ko')
plt.savefig('stroop_vs_perc_abnormal.pdf')
plt.close('all')

plt.figure()
#plt.plot(disability, percent_abnormal_subj, 'o')
plt.plot(disability, percent_abnormal_subj[indices], 'ko')
plt.plot(disability, percent_abnormal_high_subj[indices], 'ro')
plt.plot(disability, percent_abnormal_low_subj[indices], 'bo')
plt.savefig('disability_vs_perc_abnormal.pdf')
plt.close('all')

plt.figure()
plt.plot(disability, percent_abnormal_high_subj[indices] - percent_abnormal_low_subj[indices], 'bo')
plt.savefig('disability_vs_perc_highminuslowabnormal.pdf')
plt.close('all')


# groups
dis_groups = list()
dis_groups_2 = list()
dis_groups_3 = list()
disability = numpy.array(disability)
percent_abnormal_subj = numpy.array(percent_abnormal_subj)
for dis in range(0,6):
    subj = numpy.nonzero(disability == dis)[0]
    dis_groups.append(percent_abnormal_subj[subj])
    dis_groups_2.append(percent_abnormal_high_subj[subj])
    dis_groups_3.append(percent_abnormal_low_subj[subj])

plt.figure()
plt.boxplot(dis_groups)
plt.savefig('disability_vs_perc_abnormal_box.pdf')
plt.close()
plt.figure()
plt.boxplot(dis_groups_2)
plt.savefig('disability_vs_perc_high_box.pdf')
plt.close()
plt.figure()
plt.boxplot(dis_groups_3)
plt.savefig('disability_vs_perc_low_box.pdf')
plt.close()

# to remove outlier 8. WHY is FA elevated everywhere?
#indices = [0, 1, 2, 3, 4, 5, 6, 7,9, 10, 11, 12]
ret = scipy.stats.pearsonr(disability[indices], percent_abnormal_subj[indices])
print ret

#scan = [0.5, 0, 0.5, 0.5, 1, 0, nan, 0, 0, 0, 1, 0.5, 0]

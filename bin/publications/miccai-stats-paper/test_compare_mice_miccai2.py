import os
import glob
import matplotlib.pyplot as plt
import numpy
import scipy.stats
import statsmodels.sandbox.stats.multicomp as mcc
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
#parallel_jobs = 15
parallel_jobs = 10
print 'Using N jobs:', parallel_jobs

#group_indices = [1, 0, 1, 0, 0, 1, 1, 0]
# 1 T, 2 C, 3 T, 4 C, 5 C, 6 T, 7 T, 8 C

indir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/MICCAI2012/mice_with_scalars'


input_mask = "{0}/*.vtk".format(indir)
input_poly_datas = glob.glob(input_mask)

print input_poly_datas

input_pds = list()
input_pds_downsampled = list()

#number_of_fibers_per_subject = 3000
#number_of_fibers_per_subject = 1000
# this is max possible with length thresh of 15
number_of_fibers_per_subject = 3500
number_of_fiber_centroids = 1000
number_of_subjects = len(input_poly_datas)
points_per_fiber = 30
neighborhood_threshold = 15
fiber_length = 15

#input_poly_datas = input_poly_datas[0:1]
# read in ones with scalars already
for fname in input_poly_datas:
    print fname
    pd = wma.io.read_polydata(fname)
    input_pds.append(pd)

# grab scalars of interest
input_mean_mds_per_subject = list()
input_pds_downsampled = list()
for pd in input_pds:
    pd2, fiber_indices1 = wma.filter.preprocess(pd, fiber_length,return_indices=True)
    pd2, fiber_indices2 = wma.filter.downsample(pd2, number_of_fibers_per_subject,return_indices=True)
    fiber_indices = fiber_indices1[fiber_indices2]
    # get MD only at fibers of interest
    md = pd.GetCellData().GetArray('mean_MD')
    md_subj = list()
    for idx in fiber_indices:
        md_subj.append(md.GetTuple1(idx))
    #md_subj = numpy.array(md_subj)
    #md_subj = numpy.array(md_avg_list)[fiber_indices]
    input_mean_mds_per_subject.append(md_subj)    
    input_pds_downsampled.append(pd2)

# convert to arrays for dist and averaging
# use entire appended polydata (perhaps in future compute per-subject)
print 'Appending inputs into one polydata'
appender = vtk.vtkAppendPolyData()
for pd in input_pds_downsampled:
    appender.AddInput(pd)

appender.Update()

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
neighborhood_threshold = 15.0
#neighborhood_threshold = 10.0

group_indices = [1, 0, 1, 0, 0, 1, 1, 0]
group_indices = numpy.array(group_indices)

data_to_analyze = input_mean_mds_per_subject
 
# assign to flat lists with subject idx
subj_idx = 0
input_data_per_fiber = list()
input_subject_idx_per_fiber = list()
input_group_idx_per_fiber = list()
for fa_subj in data_to_analyze:
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

sigma = neighborhood_threshold/2.0
sigma_sq = sigma*sigma

# TO LOOK AT AVERAGE BRAIN, AND GROUP STATS
for idx in range(number_of_fiber_centroids):
    print idx, '/', number_of_fiber_centroids
    neighborhood_indices = numpy.nonzero(distances[idx,:] < neighborhood_threshold)[0]
    d_sq = numpy.multiply(distances[idx,neighborhood_indices],distances[idx,neighborhood_indices])
    neighborhood_weights = numpy.exp(numpy.divide(-d_sq, sigma_sq))
    
    hood_count = len(neighborhood_indices)
    # find average fiber in neighborhood
    avg_fiber = fiber_array.get_fiber(neighborhood_indices[0])
    avg_fiber *= neighborhood_weights[0]
    for idx in range(1, hood_count):
        curr_fiber = fiber_array.get_fiber(neighborhood_indices[idx])
        curr_fiber *= neighborhood_weights[idx]
        avg_fiber += curr_fiber
    #avg_fiber /= hood_count
    avg_fiber /= numpy.sum(neighborhood_weights)
    
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
        # to replace with weighted total
        subject_stats_list.append(numpy.sum(neighborhood_weights[numpy.nonzero(subject_list == sidx)[0]]))
        # this is just number of fibers
        #subject_stats_list.append(len(numpy.nonzero(subject_list == sidx)[0]))
    subject_stats_list = numpy.array(subject_stats_list)
    # figure out which subject is from which group
    g0 = numpy.nonzero(group_indices == 0)[0]
    g1 = numpy.nonzero(group_indices == 1)[0]
    t, p = scipy.stats.ttest_ind(subject_stats_list[g0], subject_stats_list[g1])
    # totals within group for visualization, error checking
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
            # weighted mean
            w_mean = numpy.multiply(neighborhood_weights[subject_fibers], data_list[subject_fibers])
            w_mean = numpy.sum(w_mean)
            w_mean = numpy.divide(w_mean, numpy.sum(neighborhood_weights[subject_fibers]))
            subject_stats_list.append(w_mean)
            # regular mean
            #subject_stats_list.append(numpy.mean(data_list[subject_fibers]))
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

f = open('pvals_Data.txt','w')
for lidx in range(len(significance_list)):
    f.write(str(significance_list[lidx]))
    f.write('\n')
f.close()

print 'load pvals.txt into matlab to test fdr for now'

plt.figure()
plt.plot(data_0_list, data_1_list, 'o')
plt.title('Data in neighborhoods')
plt.xlabel('Group 0')
plt.ylabel('Group 1')
plt.savefig('Data-groups-neighborhoods.pdf')
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
wma.io.write_polydata(outpd,'atlas_info.vtp')
# mean MD
ren = wma.render.render(outpd, scalar_bar=True, scalar_range=[0.0001,0.000497])

mask = ~numpy.isnan(data_1_list)
mask_idx = numpy.nonzero(mask)[0]
mask_pd = wma.filter.mask(outpd, mask)
mask_pd = add_array_to_polydata(mask_pd, numpy.array(data_1_list)[mask_idx], array_name='Data1')
mask_pd = add_array_to_polydata(mask_pd, numpy.array(data_0_list)[mask_idx], array_name='Data0')
mask_pd = add_array_to_polydata(mask_pd, numpy.array(data_1_list)[mask_idx]-numpy.array(data_0_list)[mask_idx], array_name='Difference')
mask_pd = add_array_to_polydata(mask_pd, numpy.array(significance_list)[mask_idx], array_name='P')
mask_pd.GetCellData().SetActiveScalars('Data0')
mask_pd.GetCellData().SetActiveScalars('Data1')

mask_pd.GetCellData().SetActiveScalars('Data0')
ren2 = wma.render.render(mask_pd, scalar_bar=True, scalar_range=[0.00019,0.00036])

mask_pd.GetCellData().SetActiveScalars('Data1')
ren2 = wma.render.render(mask_pd, scalar_bar=True, scalar_range=[0.00019,0.00036])

mask_pd.GetCellData().SetActiveScalars('Difference')
ren2 = wma.render.render(mask_pd, scalar_bar=True, scalar_range=[-.05,.05])

mask_pd.GetCellData().SetActiveScalars('P')
ren2 = wma.render.render(mask_pd, scalar_bar=True, scalar_range=[0,0.005])

# pvals = load('pvals_Data.txt');
#pvals = pvals(~isnan(pvals))
#[pID,pN] = fdr(pvals, 0.05)
#[h, crit_p, adj_p]=fdr_bh(pvals, 0.05, 'pdep', 'yes');
# hood 15, 3500 fibers/subject
thresh = 0.0167
mask_pd.GetCellData().SetActiveScalars('P')
ren2 = wma.render.render(mask_pd, scalar_bar=True, scalar_range=[0.0, thresh])

mask_colors = numpy.array(significance_list)[mask_idx]
mask_colors[mask_colors >= thresh] = numpy.nan
mask_pd2 = wma.filter.mask(mask_pd, numpy.ones(len(mask_colors)), mask_colors)
ren2 = wma.render.render(mask_pd2, scalar_bar=True, scalar_range=[0.0, thresh], colormap='hot')
actors = ren2.renderer.GetActors()
actors.GetLastActor().GetMapper().GetLookupTable().SetNanColor(0.8,0.8,0.8,0.3)


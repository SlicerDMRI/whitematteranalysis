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


# ------------------------------------------------------
# Can re-run the below after changing neighborhood_threshold
# Or after changing group membership
# All slow processing happens above this line
# ------------------------------------------------------
sigma = 25
#sigma = 15
sigma_sq = sigma*sigma
neighborhood_threshold = sigma*2

group_indices =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,\
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

group_indices = numpy.array(group_indices)


data_to_analyze = input_mean_fas_per_subject
#data_to_analyze = input_mean_perp_per_subject

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
plt.grid(1)
plt.savefig('FA-groups-neighborhoods.pdf')
plt.close()

# is change related to original FA??
plt.figure()
plt.plot(data_0_list, numpy.array(data_1_list) - numpy.array(data_0_list),'o')
plt.grid(1)
plt.savefig('FA-neighborhood-difference_vs_control_FA.pdf')
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
# mean FA
ren = wma.render.render(outpd, scalar_bar=True, scalar_range=[0.35,0.65])
# max FA
#ren = wma.render.render(outpd, scalar_bar=True, scalar_range=[0.5,0.95])
# lambda parallel
#ren = wma.render.render(outpd, scalar_bar=True, scalar_range=[0.0003,0.0011])

mask = ~numpy.isnan(data_1_list)
mask_idx = numpy.nonzero(mask)[0]
mask_pd = wma.filter.mask(outpd, mask)
mask_pd = add_array_to_polydata(mask_pd, numpy.array(data_1_list)[mask_idx], array_name='Data1')
mask_pd = add_array_to_polydata(mask_pd, numpy.array(data_0_list)[mask_idx], array_name='Data0')
mask_pd = add_array_to_polydata(mask_pd, numpy.array(data_1_list)[mask_idx]-numpy.array(data_0_list)[mask_idx], array_name='Difference')
mask_pd.GetCellData().SetActiveScalars('Data0')
mask_pd.GetCellData().SetActiveScalars('Data1')
mask_pd.GetCellData().SetActiveScalars('Difference')
ren2 = wma.render.render(mask_pd, scalar_bar=True, scalar_range=[-.05,.05])
ren2.save_views()

#======================
# test individuals
# ===== NOTE: this needs to be re-done with registration of normals in group and
## others to it, rather than registration of all datasets together for initial testing
age_control = [29, 43, 38, 31, 23, 29, 40, 24, 42, 26, 47, 23, 40]
age_tbi = [44, 37, 43, 27, 29, 42, 27, 24, 25, 29, 24, 39, 44]
age_control= numpy.array(age_control)
age_tbi= numpy.array(age_tbi)
age =  numpy.array(list(age_control) + list(age_tbi))

# according to neighborhood definition:
# compute avg fibers and FA stats in neighborhoods
average_fibers_list = list()
zscore_list = list()
density_zscore_list = list()


for idx in range(number_of_fiber_centroids):
    print idx, '/', number_of_fiber_centroids
    print idx, '/', number_of_fiber_centroids
    neighborhood_indices = numpy.nonzero(distances[idx,:] < neighborhood_threshold)[0]
    d_sq = numpy.multiply(distances[idx,neighborhood_indices],distances[idx,neighborhood_indices])
    neighborhood_weights = numpy.exp(numpy.divide(-d_sq, sigma_sq))
    hood_count = len(neighborhood_indices)

    # find average fiber in neighborhood
    # THIS SHOULD BE PER SUBJECT??
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
    # fiber density in each subject in group 1, vs group 0 model
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
    subject_density_zscore_list = list()
    for sidx in range(number_of_subjects):
        #t, p = scipy.stats.ttest_ind(subject_stats_list[g0], subject_stats_list[sidx])
        # leave one out if this is in g0
        leave_out = numpy.nonzero(g0==sidx)[0]
        if len(leave_out):
            gmodel = numpy.array(list(g0[1:leave_out]) + list(g0[leave_out+1:]))
            #print "LOO: ", sidx, "////", gmodel
        else:
            gmodel = g0
        z = scipy.stats.zmap(subject_stats_list[sidx], subject_stats_list[gmodel])
        subject_density_zscore_list.append(z)

    density_zscore_list.append(subject_density_zscore_list)    
    # statistic: FA in subject vs atlas
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
    subject_data_zscore_list = list()
    for sidx in range(number_of_subjects):
        leave_out = numpy.nonzero(g0==sidx)[0]
        if len(leave_out):
            gmodel = numpy.array(list(g0[1:leave_out]) + list(g0[leave_out+1:]))
            #print "LOO: ", sidx, "////", gmodel
        else:
            gmodel = g0
        if (len(gmodel) > 1) & ~numpy.isnan(subject_stats_list[sidx]):
            slope, intercept, r, p, err = scipy.stats.linregress(subject_stats_list[gmodel], age[gmodel])
            if (p < 0.05):
                #print "fitting line"
                #subject_stats_list[gmodel] += -intercept - numpy.multiply(slope,age[gmodel])
                #subject_stats_list[sidx] += -intercept - numpy.multiply(slope,age[sidx])
                subject_stats_list_model =  subject_stats_list[gmodel]- numpy.multiply(slope,age[gmodel])
                subject_stats_list_sidx =  subject_stats_list[sidx] - numpy.multiply(slope,age[sidx])
                z = scipy.stats.zmap(subject_stats_list_sidx, subject_stats_list_model)                
            else:
                z = scipy.stats.zmap(subject_stats_list[sidx], subject_stats_list[gmodel])
        else:
            z = numpy.nan
        subject_data_zscore_list.append(z)

    zscore_list.append(subject_data_zscore_list)
    

# view per-subject
#g1 = numpy.nonzero(group_indices == 1)[0]
#g1 = range(len(group_indices))
out_masked_pds_zscore_FA = list()
number_abnormal_high_subj = list()
number_abnormal_low_subj = list()
number_centroids_subj = list()
#abnormal_threshold = 3
abnormal_threshold = 1
mean_zscore_subj = list()

for sidx in range(number_of_subjects):
    zscore_subj = list()
    #for zscores in density_zscore_list
    for zscores in zscore_list:
        zscore_subj.append(zscores[sidx])
    zscore_subj = numpy.array(zscore_subj)
    mask = ~numpy.isnan(zscore_subj) & ~numpy.isinf(zscore_subj)
    number_abnormal_high_subj.append(numpy.sum(zscore_subj[mask] > abnormal_threshold))
    number_abnormal_low_subj.append(numpy.sum(zscore_subj[mask] < -abnormal_threshold))
    # test for just decreases
    #number_abnormal_subj.append(numpy.sum(zscore_subj[mask] < -abnormal_threshold))
    number_centroids_subj.append(numpy.sum(mask))
    out_masked_pds_zscore_FA.append(wma.filter.mask(outpd, mask, zscore_subj))
    mean_zscore_subj.append(numpy.mean(numpy.abs(zscore_subj[mask])))
    # below is very unstable
    #mean_zscore_subj.append(numpy.mean(numpy.multiply(zscore_subj[mask],zscore_subj[mask])))
    #mean_zscore_subj.append(numpy.max(numpy.abs(zscore_subj[mask])))
    #mean_zscore_subj.append(scipy.stats.tmax(numpy.abs(zscore_subj[mask]), 10))
                            
number_abnormal_high_subj = numpy.array(number_abnormal_high_subj)
number_abnormal_low_subj = numpy.array(number_abnormal_low_subj)
number_abnormal_subj = number_abnormal_high_subj + number_abnormal_low_subj

percent_abnormal_subj = numpy.divide(number_abnormal_subj.astype(float), numpy.array(number_centroids_subj))
percent_abnormal_high_subj = numpy.divide(number_abnormal_high_subj.astype(float), numpy.array(number_centroids_subj))
percent_abnormal_low_subj = numpy.divide(number_abnormal_low_subj.astype(float), numpy.array(number_centroids_subj))

plt.figure()
plt.plot(age,mean_zscore_subj,'o')
plt.plot(age_tbi,mean_zscore_subj[13:],'ro')
plt.xlabel('age')
plt.ylabel('mean abs z score tract FA')
#plt.ylabel('max abs z score tract FA')
plt.title('Controls (blue) vs. mTBI (red)')
plt.savefig('age_vs_z_nLandPt.pdf')
plt.close('all')


plt.figure()
plt.plot(age,percent_abnormal_subj,'o')
plt.plot(age_tbi,percent_abnormal_subj[13:],'ro')
plt.xlabel('age')
plt.ylabel('percent abnormal tract FA')
plt.title('Controls (blue) vs. mTBI (red)')
plt.savefig('age_vs_perc_ab_nLandPt.pdf')


plt.figure()
plt.plot(age,percent_abnormal_low_subj,'o')
plt.plot(age_tbi,percent_abnormal_low_subj[13:],'ro')
plt.xlabel('age')
plt.ylabel('percent decrease tract FA')
plt.title('Controls (blue) vs. mTBI (red)')
plt.savefig('age_vs_perc_lowab_nLandPt.pdf')


plt.figure()
plt.plot(age,percent_abnormal_high_subj,'o')
plt.plot(age_tbi,percent_abnormal_high_subj[13:],'ro')
plt.xlabel('age')
plt.ylabel('percent increase tract FA')
plt.title('Controls (blue) vs. mTBI (red)')
plt.savefig('age_vs_perc_highab_nLandPt.pdf')
plt.close('all')
print "LAUREN why are there -inf z scores??? are these bad samples from tbi only?"

ret = scipy.stats.ttest_ind(percent_abnormal_subj[0:13],percent_abnormal_subj[13:])
print ret

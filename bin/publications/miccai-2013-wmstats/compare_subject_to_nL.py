# assume we just ran setup and compare to normal

# for each subject, does its FA increase/decrease where it was low/high?
# i.e. does the group pattern seem to hold true in individuals?

# data structures for storing FA per subject
subject_fa = numpy.zeros([number_of_subjects, number_of_fiber_centroids])
subject_fa_age = numpy.zeros([number_of_subjects, number_of_fiber_centroids])
# mean to compare to. different per subject because of leave-one-out for controls
control_mean_fa = numpy.zeros([number_of_subjects, number_of_fiber_centroids])
control_mean_fa_age = numpy.zeros([number_of_subjects, number_of_fiber_centroids])

# now loop over neighborhoods and compute FA per subject
# also compute difference from normal model, taking age into account
for idx in range(number_of_fiber_centroids):

    neighborhood_indices = numpy.nonzero(distances[idx,:] < neighborhood_threshold)[0]
    d_sq = numpy.multiply(distances[idx,neighborhood_indices],distances[idx,neighborhood_indices])
    neighborhood_weights = numpy.exp(numpy.divide(-d_sq, sigma_sq))
    hood_count = len(neighborhood_indices)

    # compute statistic(s) in neighborhood
    data_list = list()
    #group_list = list()
    subject_list = list()
    for hood_idx in neighborhood_indices:
        data_list.append(input_data_per_fiber[hood_idx])
        #group_list.append(input_group_idx_per_fiber[hood_idx])
        subject_list.append(input_subject_idx_per_fiber[hood_idx])

    # statistic: FA in subject vs atlas
    data_list = numpy.array(data_list)
    subject_list = numpy.array(subject_list)
    
    for sidx in range(number_of_subjects):
        print idx, '/', number_of_fiber_centroids, '.....', sidx, '/', number_of_subjects
        subject_fibers = numpy.nonzero(subject_list == sidx)[0]
        if len(subject_fibers):
            # weighted mean
            w_mean = numpy.multiply(neighborhood_weights[subject_fibers], data_list[subject_fibers])
            w_mean = numpy.sum(w_mean)
            w_mean = numpy.divide(w_mean, numpy.sum(neighborhood_weights[subject_fibers]))
            subject_fa[sidx, idx] = w_mean
            # regular mean
            #subject_fa[sidx].append(numpy.mean(data_list[subject_fibers]))
        else:
           subject_fa[sidx, idx] = numpy.nan

    # correct for age
    # build model in this neighborhood ignoring any nan values
    g0 = numpy.nonzero((group_indices == 0) & ~numpy.isnan(subject_stats_list))[0]
    subject_data_zscore_list = list()
    # find control group leaving out current subject if control
    for sidx in range(number_of_subjects):
        leave_out = numpy.nonzero(g0==sidx)[0]
        if len(leave_out):
            gmodel = numpy.array(list(g0[1:leave_out]) + list(g0[leave_out+1:]))
            #print "LOO: ", sidx, "////", gmodel
        else:
            gmodel = g0

        if (len(gmodel) > 1) & ~numpy.isnan(subject_fa[sidx, idx]):
            control_mean_fa[sidx, idx] =  numpy.mean(subject_fa[gmodel, idx])
            slope, intercept, r, p, err = scipy.stats.linregress(age[gmodel], subject_fa[gmodel, idx])
            if (p < 0.05):
                print "fitting line"
                #subject_stats_list[gmodel] += -intercept - numpy.multiply(slope,age[gmodel])
                #subject_stats_list[sidx] += -intercept - numpy.multiply(slope,age[sidx])
                
                control_mean_fa_age[sidx, idx] =  numpy.mean(subject_fa[gmodel, idx] - numpy.multiply(slope,age[gmodel]))
                subject_fa_age[sidx, idx] =  subject_fa[sidx, idx] - numpy.multiply(slope,age[sidx])
                #z = scipy.stats.zmap(subject_stats_list_sidx, subject_stats_list_model)                
            else:
                control_mean_fa_age[sidx, idx] = control_mean_fa[sidx, idx]
                subject_fa_age[sidx, idx] =  subject_fa[sidx, idx]
                #z = scipy.stats.zmap(subject_stats_list[sidx], subject_stats_list[gmodel])
        else:
            #z = numpy.nan
            control_mean_fa[sidx, idx] = numpy.nan
            control_mean_fa_age[sidx, idx] = numpy.nan
            subject_fa_age[sidx, idx] = numpy.nan

            
#  what are changes in FA?
delta_fa = subject_fa_age - control_mean_fa_age

test = numpy.nanmax(delta_fa, axis=1) - numpy.nanmin(delta_fa, axis=1) 

for sidx in range(number_of_subjects):
    plt.figure()
    plt.plot(control_mean_fa_age[sidx,:],delta_fa[sidx,:],'o')
    plt.savefig('test{0}.pdf'.format(sidx))

g0 = numpy.nonzero(group_indices == 0)[0]
g1 = numpy.nonzero(group_indices == 1)[0]

plt.figure()
for sidx in g0:
    plt.plot(control_mean_fa_age[sidx,:],delta_fa[sidx,:],'bo')
for sidx in g1:
    plt.plot(control_mean_fa_age[sidx,:],delta_fa[sidx,:],'o')
plt.savefig('test_all.pdf'.format(sidx))


#  what are changes in FA?
delta_fa_uncorrected = subject_fa - control_mean_fa

for sidx in range(number_of_subjects):
    plt.figure()
    plt.plot(control_mean_fa[sidx,:],delta_fa_uncorrected[sidx,:],'o')
    plt.savefig('test_uncorr_{0}.pdf'.format(sidx))

g0 = numpy.nonzero(group_indices == 0)[0]
g1 = numpy.nonzero(group_indices == 1)[0]

plt.figure()
for sidx in g0:
    plt.plot(control_mean_fa[sidx,:],delta_fa_uncorrected[sidx,:],'bo')
for sidx in g1:
    plt.plot(control_mean_fa[sidx,:],delta_fa_uncorrected[sidx,:],'o')
plt.savefig('test_all_uncorr.pdf'.format(sidx))


# looks like should plot the line delta fa vs fa, build model of that in normals, i.e. what
# does the smoothing in the group do vs individuals, then look for who is away from that a
# lot. can have a z-score based on normal model for each FA value.
# this however will ignore local variability in the model in the neighborhood. may not be needed?
# for now, simpler, look for regions more than 0.09 away from mean, and count how many.

test = numpy.abs(delta_fa_uncorrected) > 0.08
outcome = numpy.sum(test, axis=1)

plt.figure()
plt.boxplot([outcome[g0], outcome[g1]])
plt.savefig('diagnose_maybe.pdf'.format(sidx))


# now look at the best idea so far... what is the expected mean and standard deviation of the
# difference from 'normal', versus FA.  model that and measure deviations from it.
#bin_width = 0.025
bin_width = 0.01
half_bin_width = bin_width / 2.0
fa_axis = numpy.arange(0.25, 0.75, bin_width)

difference_mean = list()
difference_std = list()

for fa_val in fa_axis:
    # model control_mean_fa data "belonging" to this bin
    indices = numpy.where(numpy.abs(control_mean_fa[g0,:] - fa_val) <= half_bin_width)
    difference_mean.append(numpy.mean(delta_fa_uncorrected[g0][indices]))
    difference_std.append(numpy.std(delta_fa_uncorrected[g0][indices]))

difference_mean = numpy.array(difference_mean)
difference_std = numpy.array(difference_std)

plt.figure()
plt.errorbar(fa_axis, difference_mean, difference_std)
plt.savefig('difference_model2.pdf')

# for subject data, compare every delta to the model. find a z-score per delta
zscore_delta_fa = numpy.zeros(delta_fa_uncorrected.shape)
diff_delta_fa = numpy.zeros(delta_fa_uncorrected.shape)
mean_z = numpy.zeros([number_of_subjects, len(fa_axis)])

for fa_index in range(0,len(fa_axis)):
    indices = numpy.where(numpy.abs(control_mean_fa - fa_axis[fa_index]) <= half_bin_width)
    zscore_delta_fa[indices] = numpy.divide((delta_fa_uncorrected[indices] - difference_mean[fa_index]), difference_std[fa_index])
    diff_delta_fa[indices] = delta_fa_uncorrected[indices] - difference_mean[fa_index]
    # per-subject mean z for each fa bin of interest
    for sidx in range(number_of_subjects):
        idxcol = numpy.where(indices[0] == sidx)
        mean_z[sidx,fa_index] = numpy.mean(zscore_delta_fa[sidx][indices[1][idxcol]]) 


plt.figure()
plt.plot(fa_axis, mean_z[g0].T, 'bo')
plt.plot(fa_axis, mean_z[g1].T, 'o')
plt.savefig('meanz_vsfa.pdf')

#test = numpy.abs(mean_z) > 0.5
indices = fa_axis < 0.4
test1 = mean_z[:, indices] > 1
indices = fa_axis > 0.55
test2 = mean_z[:,indices] < -0.3
#outcome = numpy.sum(test1, axis=1) + numpy.sum(test2, axis=1)
#outcome = numpy.sum(test, axis=1)
outcome = numpy.sum(test1, axis=1)
outcome[g0]
outcome[g1]



test = numpy.abs(zscore_delta_fa) > 3.0
outcome = numpy.sum(test, axis=1)
#outcome = numpy.sum(numpy.multiply(test,numpy.abs(zscore_delta_fa)), axis=1)
outcome[g0]
outcome[g1]

plt.figure()
plt.boxplot([outcome[g0], outcome[g1]])
plt.savefig('diagnose_maybe.pdf'.format(sidx))

plt.figure()
plt.plot(age[g0], outcome[g0],'bo')
plt.plot(age[g1], outcome[g1],'ro')
plt.savefig('diagnose_maybe2.pdf'.format(sidx))


test2 = numpy.abs(diff_delta_fa) > 0.08
outcome = numpy.sum(test2, axis=1)

plt.figure()
plt.boxplot([outcome[g0], outcome[g1]])
plt.savefig('diagnose_maybe3.pdf'.format(sidx))


mask2 = ~numpy.isnan(data_1_list) & (numpy.abs(numpy.array(data_1_list)-numpy.array(data_0_list))> 0.015)
mask_diff_pd = mask_pd = wma.filter.mask(outpd, mask2)
mask_pd = add_array_to_polydata(mask_pd, numpy.array(data_1_list)[mask2]-numpy.array(data_0_list)[mask2], array_name='Difference')
ren3 = wma.render.render(mask_diff_pd, scalar_bar=True, scalar_range=[-.05,.05])

# look at data in mask

plt.figure()
test_data = (delta_fa_uncorrected[:,mask2])
plt.plot(test_data[g0,:].T, 'bo')
plt.plot(test_data[g1,:].T, 'ro')
plt.savefig('feature_ident.pdf')

test = numpy.abs(test_data) > 0.15
outcome = numpy.sum(test, axis=1)
outcome[g0]
outcome[g1]


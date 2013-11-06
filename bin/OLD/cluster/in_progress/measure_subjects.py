import whitematteranalysis as wma
import vtk
import numpy
import matplotlib.pyplot as plt
#import multiprocessing
import glob
#import sys
#import os
import scipy.stats

number_of_clusters = 394
cluster_ids = range(0, number_of_clusters)
#indir = '/Users/lauren/Dropbox/Coding/OUTPUTS/March2013/TBI/cluster_4'
indir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/March2013/TBI/cluster_5'
input_mask = "{0}/*.vtp".format(indir)
input_poly_datas = glob.glob(input_mask)

#input_poly_datas = input_poly_datas[0:4]
#group_indices = [0,0,1,1]

group_indices =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,\
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
group_indices = numpy.array(group_indices)
g0 = numpy.where(group_indices == 0)
g1 = numpy.where(group_indices == 1)
# ===== NOTE: this needs to be re-done with registration of normals in group and
## others to it, rather than registration of all datasets together for initial testing
age_control = [29, 43, 38, 31, 23, 29, 40, 24, 42, 26, 47, 23, 40]
age_tbi = [44, 37, 43, 27, 29, 42, 27, 24, 25, 29, 24, 39, 44]
age_control= numpy.array(age_control)
age_tbi= numpy.array(age_tbi)
age =  numpy.array(list(age_control) + list(age_tbi))

# read them in, measure stuff in each cluster
subject_cluster_mean_fa = list()
subject_cluster_mean_para = list()
subject_cluster_mean_perp = list()
subject_cluster_std_fa = list()
subject_cluster_std_para = list()
subject_cluster_std_perp = list()

for fname in input_poly_datas:
    print fname
    pd = wma.io.read_polydata(fname)
    
    # this code assumes for now we have the following arrays
    #Array 0 name = FA
    #Array 1 name = Tensors_
    #Array 2 name = parallel_diffusivity
    #Array 3 name = perpendicular_diffusivity
    #Array 4 name = ClusterNumber_PointData
    
    fa = pd.GetPointData().GetArray('FA')
    para = pd.GetPointData().GetArray('parallel_diffusivity')
    perp = pd.GetPointData().GetArray('perpendicular_diffusivity')
    clusters = pd.GetPointData().GetArray('ClusterNumber_PointData')
    cluster_data_fa = list()
    cluster_data_para = list()
    cluster_data_perp = list()

    cl_range = clusters.GetRange()
    num_points = clusters.GetNumberOfTuples()
    #cluster_ids = range(int(cl_range[0]), int(cl_range[1]) + 1)
    # ignore negative ids (junk cluster is -1)
    # this can give an error if a subject is missing the last cluster.
    #cluster_ids = range(0, int(cl_range[1]) + 1)
    
    # make a list to hold all values in each cluster
    for cidx in cluster_ids:
        cluster_data_fa.append(list())
        cluster_data_para.append(list())
        cluster_data_perp.append(list())

    # put data from each cluster into the lists
    for pidx in range(num_points):
        cidx = int(clusters.GetTuple1(pidx))
        # ignore negative cluster ids
        if cidx >= 0:
            cluster_data_fa[cidx].append(fa.GetTuple1(pidx))
            cluster_data_para[cidx].append(para.GetTuple1(pidx))
            cluster_data_perp[cidx].append(perp.GetTuple1(pidx))

    # compute mean and standard deviation for now
    cluster_mean_fa = list()
    cluster_mean_para = list()
    cluster_mean_perp = list()
    cluster_std_fa = list()
    cluster_std_para = list()
    cluster_std_perp = list()

    for cidx in cluster_ids:
        cluster_mean_fa.append(numpy.mean(numpy.array(cluster_data_fa[cidx])))
        cluster_std_fa.append(numpy.std(numpy.array(cluster_data_fa[cidx])))
        cluster_mean_para.append(numpy.mean(numpy.array(cluster_data_para[cidx])))
        cluster_std_para.append(numpy.std(numpy.array(cluster_data_para[cidx])))
        cluster_mean_perp.append(numpy.mean(numpy.array(cluster_data_perp[cidx])))
        cluster_std_perp.append(numpy.std(numpy.array(cluster_data_perp[cidx])))

    # record the values
    subject_cluster_mean_fa.append(numpy.array(cluster_mean_fa))
    subject_cluster_mean_para.append(numpy.array(cluster_mean_para))
    subject_cluster_mean_perp.append(numpy.array(cluster_mean_perp))
    subject_cluster_std_fa.append(numpy.array(cluster_std_fa))
    subject_cluster_std_para.append(numpy.array(cluster_std_para))
    subject_cluster_std_perp.append(numpy.array(cluster_std_perp))


# plot something
subject_ids = range(len(input_poly_datas))
for sidx in subject_ids:
    if group_indices[sidx] == 0:
        plt.plot(cluster_ids, subject_cluster_mean_fa[sidx], 'bo')
    else:
        plt.plot(cluster_ids, subject_cluster_mean_fa[sidx], 'ro')

plt.savefig('mean_fa.pdf')
plt.close()

# test stuff
subject_cluster_mean_fa = numpy.array(subject_cluster_mean_fa)
cluster_stat_fa = list()
for cidx in cluster_ids:
    print cidx
    t,p = scipy.stats.ttest_ind(subject_cluster_mean_fa[g0,cidx][0], subject_cluster_mean_fa[g1,cidx][0])
    cluster_stat_fa.append(p)

# plot something
subject_ids = range(len(input_poly_datas))
for sidx in subject_ids:
    if group_indices[sidx] == 0:
        plt.plot(cluster_ids, subject_cluster_mean_para[sidx], 'bo')
    else:
        plt.plot(cluster_ids, subject_cluster_mean_para[sidx], 'ro')

plt.savefig('mean_para.pdf')
plt.close()

# test stuff
subject_cluster_mean_para = numpy.array(subject_cluster_mean_para)
cluster_stat_para = list()
for cidx in cluster_ids:
    print cidx
    t,p = scipy.stats.ttest_ind(subject_cluster_mean_para[g0,cidx][0], subject_cluster_mean_para[g1,cidx][0])
    cluster_stat_para.append(p)

# plot something
subject_ids = range(len(input_poly_datas))
for sidx in subject_ids:
    if group_indices[sidx] == 0:
        plt.plot(cluster_ids, subject_cluster_mean_perp[sidx], 'bo')
    else:
        plt.plot(cluster_ids, subject_cluster_mean_perp[sidx], 'ro')

plt.savefig('mean_perp.pdf')
plt.close()

# test stuff
subject_cluster_mean_perp = numpy.array(subject_cluster_mean_perp)
cluster_stat_perp = list()
for cidx in cluster_ids:
    print cidx
    t,p = scipy.stats.ttest_ind(subject_cluster_mean_perp[g0,cidx][0], subject_cluster_mean_perp[g1,cidx][0])
    cluster_stat_perp.append(p)


# --var

# plot something
subject_ids = range(len(input_poly_datas))
for sidx in subject_ids:
    if group_indices[sidx] == 0:
        plt.plot(cluster_ids, subject_cluster_std_fa[sidx], 'bo')
    else:
        plt.plot(cluster_ids, subject_cluster_std_fa[sidx], 'ro')

plt.savefig('std_fa.pdf')
plt.close()

# test stuff
subject_cluster_std_fa = numpy.array(subject_cluster_std_fa)
cluster_stat_fa = list()
for cidx in cluster_ids:
    print cidx
    t,p = scipy.stats.ttest_ind(subject_cluster_std_fa[g0,cidx][0], subject_cluster_std_fa[g1,cidx][0])
    cluster_stat_fa.append(p)

# plot something
subject_ids = range(len(input_poly_datas))
for sidx in subject_ids:
    if group_indices[sidx] == 0:
        plt.plot(cluster_ids, subject_cluster_std_para[sidx], 'bo')
    else:
        plt.plot(cluster_ids, subject_cluster_std_para[sidx], 'ro')

plt.savefig('std_para.pdf')
plt.close()

# test stuff
subject_cluster_std_para = numpy.array(subject_cluster_std_para)
cluster_stat_para = list()
for cidx in cluster_ids:
    print cidx
    t,p = scipy.stats.ttest_ind(subject_cluster_std_para[g0,cidx][0], subject_cluster_std_para[g1,cidx][0])
    cluster_stat_para.append(p)

# plot something
subject_ids = range(len(input_poly_datas))
for sidx in subject_ids:
    if group_indices[sidx] == 0:
        plt.plot(cluster_ids, subject_cluster_std_perp[sidx], 'bo')
    else:
        plt.plot(cluster_ids, subject_cluster_std_perp[sidx], 'ro')

plt.savefig('std_perp.pdf')
plt.close()

# test stuff
subject_cluster_std_perp = numpy.array(subject_cluster_std_perp)
cluster_stat_perp = list()
for cidx in cluster_ids:
    print cidx
    t,p = scipy.stats.ttest_ind(subject_cluster_std_perp[g0,cidx][0], subject_cluster_std_perp[g1,cidx][0])
    cluster_stat_perp.append(p)

# -- var

### COMPUTE Z SCORES ######
threshold = 2.0
# PERP ###
subject_z_scores_perp = list()
for sidx in subject_ids:
    subject_z_scores_perp.append(list())
    controls = set(g0[0])
    # leave current one out if it's a control
    controls.discard(sidx)
    controls = numpy.array(list(controls))
    for cidx in cluster_ids:
        c_mean_perp = numpy.mean(subject_cluster_mean_perp[controls,cidx])
        c_std_perp = numpy.std(subject_cluster_mean_perp[controls,cidx])
        subject_z_scores_perp[sidx].append(numpy.divide(subject_cluster_mean_perp[sidx,cidx] - c_mean_perp, c_std_perp))

mean_z_perp = list()
num_z_perp = list()
for sidx in subject_ids:
    data = numpy.array(subject_z_scores_perp[sidx])
    data = numpy.ma.masked_array(data, numpy.isnan(data))
    mean_z_perp.append(numpy.mean(numpy.multiply(data, data)))
    num_z_perp.append(len(numpy.where(data >= threshold)[0]))
    
mean_z_perp = numpy.array(mean_z_perp)
num_z_perp = numpy.array(num_z_perp)

plt.figure()
plt.boxplot((mean_z_perp[g0], mean_z_perp[g1]))
plt.savefig('mean-z2-perp.pdf')
plt.close()

plt.figure()
plt.boxplot((num_z_perp[g0], num_z_perp[g1]))
plt.savefig('num-z2-perp.pdf')
plt.close()


## FA
subject_z_scores_fa = list()
for sidx in subject_ids:
    subject_z_scores_fa.append(list())
    controls = set(g0[0])
    # leave current one out if it's a control
    controls.discard(sidx)
    controls = numpy.array(list(controls))
    for cidx in cluster_ids:
        c_mean_fa = numpy.mean(subject_cluster_mean_fa[controls,cidx])
        c_std_fa = numpy.std(subject_cluster_mean_fa[controls,cidx])
        subject_z_scores_fa[sidx].append(numpy.divide(subject_cluster_mean_fa[sidx,cidx] - c_mean_fa, c_std_fa))

mean_z_fa = list()
num_z_fa = list()
for sidx in subject_ids:
    data = numpy.array(subject_z_scores_fa[sidx])
    data = numpy.ma.masked_array(data, numpy.isnan(data))
    mean_z_fa.append(numpy.mean(numpy.multiply(data, data)))
    num_z_fa.append(len(numpy.where(data >= threshold)[0]))
    
mean_z_fa = numpy.array(mean_z_fa)
num_z_fa = numpy.array(num_z_fa)

plt.figure()
plt.boxplot((mean_z_fa[g0], mean_z_fa[g1]))
plt.savefig('mean-z2-fa.pdf')
plt.close()

plt.figure()
plt.boxplot((num_z_fa[g0], num_z_fa[g1]))
plt.savefig('num-z2-fa.pdf')
plt.close()


### PARA
subject_z_scores_para = list()
for sidx in subject_ids:
    subject_z_scores_para.append(list())
    controls = set(g0[0])
    # leave current one out if it's a control
    controls.discard(sidx)
    controls = numpy.array(list(controls))
    for cidx in cluster_ids:
        c_mean_para = numpy.mean(subject_cluster_mean_para[controls,cidx])
        c_std_para = numpy.std(subject_cluster_mean_para[controls,cidx])
        subject_z_scores_para[sidx].append(numpy.divide(subject_cluster_mean_para[sidx,cidx] - c_mean_para, c_std_para))

mean_z_para = list()
num_z_para = list()
for sidx in subject_ids:
    data = numpy.array(subject_z_scores_para[sidx])
    data = numpy.ma.masked_array(data, numpy.isnan(data))
    mean_z_para.append(numpy.mean(numpy.multiply(data, data)))
    num_z_para.append(len(numpy.where(data >= threshold)[0]))
    
mean_z_para = numpy.array(mean_z_para)
num_z_para = numpy.array(num_z_para)

plt.figure()
plt.boxplot((mean_z_para[g0], mean_z_para[g1]))
plt.savefig('mean-z2-para.pdf')
plt.close()

plt.figure()
plt.boxplot((num_z_para[g0], num_z_para[g1]))
plt.savefig('num-z2-para.pdf')
plt.close()


# what if we sum various measures
sum_mean_z = mean_z_para + mean_z_perp + mean_z_fa
sum_num_z = num_z_para + num_z_perp + num_z_fa

plt.figure()
plt.boxplot((sum_mean_z[g0], sum_mean_z[g1]))
plt.savefig('sum_mean_z.pdf')
plt.close()

plt.figure()
plt.boxplot((sum_num_z[g0], sum_num_z[g1]))
plt.savefig('sum_num_z.pdf')
plt.close()

scipy.stats.ttest_ind(sum_num_z[g0], sum_num_z[g1])

plt.figure()
plt.plot(age[g0], sum_num_z[g0], 'o')
plt.plot(age[g1], sum_num_z[g1], 'ro')
plt.savefig('sum_num_z_vs_age.pdf')
plt.close()

# -----------------------------------
# fa difference vs original fa???
mean_fa_cluster_control = list()
for cidx in cluster_ids:
    print cidx
    mean_fa_cluster_control.append(numpy.mean(subject_cluster_mean_fa[g0,cidx][0]))

#, subject_cluster_mean_fa[g1,cidx][0])
plt.figure()
for sidx in subject_ids:
    if group_indices[sidx] == 0:
        plt.plot(mean_fa_cluster_control, subject_cluster_mean_fa[sidx,:], 'bo')
    else:
        plt.plot(mean_fa_cluster_control, subject_cluster_mean_fa[sidx,:], 'or')
plt.savefig('fa_group_vs_mean.pdf')
plt.close()


# Perhaps looking at missing clusters is the way to go!
# Perhaps more sensitive with higher number of clusters, or not.
# nans?
nans = list()
for sidx in subject_ids:
    nans.append(numpy.sum(numpy.isnan(subject_cluster_mean_fa[sidx,:])))
nans=numpy.array(nans)
print scipy.stats.ttest_ind(nans[g0], nans[g1])


#

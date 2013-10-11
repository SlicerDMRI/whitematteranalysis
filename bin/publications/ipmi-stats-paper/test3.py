import glob
import numpy
import scipy.stats
import vtk
import statsmodels.sandbox.stats.multicomp as mcc
import whitematteranalysis as wma

# ================= SETTINGS HERE  =================
indir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/Dec2012/ipmi-reg-gender/iteration_6'
#outdir = 
number_of_fibers_per_subject = 300
#number_of_fiber_centroids = 100
number_of_fiber_centroids = 500
points_per_fiber = 30
neighborhood_sigma = 20.0
#neighborhood_threshold = 20.0
# too small.(gave sigma 5)
#neighborhood_threshold = 10.0

minimum_fiber_length = 30
#group_indices = numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
group_indices = numpy.array([0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1])
g0 = numpy.nonzero(group_indices == 0)
g1 = numpy.nonzero(group_indices == 1)
# ================= END OF SETTINGS =================


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
    island = numpy.zeros([number_of_fiber_centroids, 1]) - 100
    curr_island = 0
    # loop over centroids, agglomerating nearest neighborhoods
    for cidx in range(number_of_fiber_centroids):
        #print cidx
        if statistic[cidx] < threshold:
            # start a new island if needed
            if island[cidx] == -100:
                island[cidx] = curr_island
                curr_island += 1
            print cidx, island[cidx]
            # find nearest neighbors that are centroids. skip self (first one)
            centroids_sorted = numpy.argsort(distances[cidx,:])[1:]
            growing = 1
            candidate = 0
            number_of_candidates = len(centroids_sorted)
            while growing & (candidate < number_of_candidates):
                cidx2 = centroids_sorted[candidate]
                # if this one can be agglomerated with the current island do it
                if statistic[cidx2] < threshold:
                    #print cidx, island[cidx], ': ?', island[cidx2]
                    # if this one is unassigned so far assign to the current island
                    if island[cidx2] == -100:
                        print "join: ", island[cidx2], island[cidx], cidx, cidx2
                        island[cidx2] = island[cidx]
                    elif island[cidx2] != island[cidx]:
                        # otherwise we have to agglomerate the existing island with this one.
                        merge_idx = island[cidx2]
                        print "agglomerate: ", merge_idx, island[cidx], cidx, cidx2
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
    return island, island_ids, island_sizes

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

# assign to flat lists with subject idx
subj_idx = 0
input_subject_idx_per_fiber = list()
input_group_idx_per_fiber = list()
for subj_idx in range(number_of_subjects):
    for fiber_idx in range(number_of_fibers_per_subject):
        input_subject_idx_per_fiber.append(subj_idx)
        input_group_idx_per_fiber.append(group_indices[subj_idx])
    subj_idx +=1
    
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
distances = numpy.zeros([number_of_fiber_centroids, total_number_of_fibers])

for idx in range(number_of_fiber_centroids):
    print idx, '/', number_of_fiber_centroids
    fiber = fiber_array.get_fiber(fiber_sample[idx])
    distances[idx,:] = wma.similarity.fiber_distance(fiber, fiber_array, threshold=0, distance_method='Hausdorff')

distances_centroids = distances[:,fiber_sample]

# ================= compute something simple to test, # fibers in each 'hood
#sigma = neighborhood_threshold/2.0
sigma = neighborhood_sigma
sigma_sq = sigma*sigma
neighborhood_threshold = sigma*2
average_fibers_list = list()
data_0_list = list()
data_1_list = list()
data_2_list = list()
data_3_list = list()

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

    # count neighborhood size per subject
    #data_list = list()
    group_list = list()
    subject_list = list()
    for hood_idx in neighborhood_indices:
        #data_list.append(input_data_per_fiber[hood_idx])
        group_list.append(input_group_idx_per_fiber[hood_idx])
        subject_list.append(input_subject_idx_per_fiber[hood_idx])
    data_0_list.append(numpy.divide(numpy.sum(numpy.array(group_list)), float(len(group_list))))
    data_1_list.append(len(group_list))
    data_2_list.append(len(numpy.unique(subject_list)))
    subject_counts = numpy.bincount(numpy.array(subject_list), minlength=number_of_subjects)
    [t, p] = scipy.stats.ttest_ind(subject_counts[g0], subject_counts[g1])
    data_3_list.append(p)

    #
# correct p values
[reject, data_4, alpha1, alpha2] = mcc.multipletests(numpy.array(data_3_list), 0.1)

# output as pd
outpd = fiber_list_to_fiber_array(average_fibers_list).convert_to_polydata()
outpd = add_array_to_polydata(outpd, data_1_list, array_name='Data1')
outpd = add_array_to_polydata(outpd, data_2_list, array_name='Data2')
outpd = add_array_to_polydata(outpd, data_0_list, array_name='Data0')
outpd = add_array_to_polydata(outpd, data_3_list, array_name='Data3')
outpd = add_array_to_polydata(outpd, data_4, array_name='Data4')


ren = wma.render.render(outpd, scalar_bar=True, scalar_range=[0,1])

# mask pd if outlier hood, less than half of subjects
#mask = numpy.array(data_2_list) > 9
mask = numpy.array(data_2_list) > 3
#mask = ( numpy.array(data_2_list) > 9 ) & (abs(numpy.array(data_0_list)-0.5)>0.05)
#mask = ( numpy.array(data_2_list) > 9 ) & (numpy.array(data_3_list) < 0.05)
outpd_masked = wma.filter.mask(outpd, mask, numpy.array(data_0_list))
#outpd_masked = wma.filter.mask(outpd, mask, numpy.array(data_3_list))
#outpd_masked = wma.filter.mask(outpd, mask, numpy.array(data_4))

ren2 = wma.render.render(outpd_masked, scalar_bar=True, scalar_range=[0,1])

# ================= compute KL divergence
# sum ln (p(f)/q(f))
# sum ln (q(f)/p(f))


#--- other approach
# embed using landmark MDS formulae


# plot the first two dimensions, color by group?

# discretize space, make N-D histogram per subject

# compute KL divergence

#========
# test islands
statistic = numpy.array(data_3_list)
statistic_threshold = 0.15
distance_threshold = 6
island_labels, island_ids, island_sizes = find_islands_and_sizes(distances_centroids, distance_threshold, statistic, statistic_threshold)

mask=numpy.ones(len(island_labels))
pd_island = wma.filter.mask(outpd, mask, island_labels)
ren2 = wma.render.render(pd_island, scalar_bar=True, scalar_range=[-10, max(island_ids)])

numpy.argsort(island_sizes)
print numpy.sort(island_sizes)
# maximum island size, not background island (-100)
lg_island = island_ids[numpy.argsort(island_sizes)[-2]][0]
mask = numpy.array(island_labels==lg_island)
pd_island = wma.filter.mask(outpd, mask, island_labels)


import glob
import os
import numpy
import vtk
import whitematteranalysis as wma
try:
    from joblib import Parallel, delayed
    USE_PARALLEL = 1
except ImportError:
    USE_PARALLEL = 0
    print "<compute_statistics.py> Failed to import joblib, cannot multiprocess."
    print " Please install joblib for this functionality."

# ================= SETTINGS HERE  =================
indir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/Dec2012/ipmi-reg-fe2/iteration_5'
outdir = './fe_downsampled'
#number_of_fibers_per_subject = 300
#number_of_fibers_per_subject = 600
number_of_fibers_per_subject = 4000
#number_of_fiber_centroids = 100
#number_of_fiber_centroids = 500
#number_of_fiber_centroids = 1500
number_of_fiber_centroids = 2000
points_per_fiber = 30

minimum_fiber_length = 30
n_jobs = 20
# ================= END OF SETTINGS =================

def read_and_process_pd(fname, outdir):
    print fname
    pd = wma.io.read_polydata(fname)
    pd = wma.filter.preprocess(pd, minimum_fiber_length)
    pd = wma.filter.downsample(pd, number_of_fibers_per_subject)
    # can't return pd (can't pickle it) so write to disk
    outfile = os.path.splitext(os.path.basename(fname))[0] + '_downsampled.vtk'
    outfile = os.path.join(outdir, outfile)
    wma.io.write_polydata(pd, outfile)
    del pd
    return outfile

# ================= read data
input_mask = "{0}/*.vtk".format(indir)
input_poly_datas = glob.glob(input_mask)
#input_poly_datas = input_poly_datas[0:3]
number_of_subjects = len(input_poly_datas)
print input_poly_datas

input_pds = list()
# parallel input preprocessing
if USE_PARALLEL:
    input_processed_poly_datas = Parallel(n_jobs=n_jobs, verbose=1)(
        delayed(read_and_process_pd)(fname, outdir)
        for fname in input_poly_datas)
    # now read them in (smaller files, faster). why is this still slow??
    for fname in input_processed_poly_datas:
        print fname
        pd = wma.io.read_polydata(fname)
        input_pds.append(pd)
else:
    for fname in input_poly_datas:
        print fname
        pd = wma.io.read_polydata(fname)
        pd = wma.filter.preprocess(pd, minimum_fiber_length)
        pd = wma.filter.downsample(pd, number_of_fibers_per_subject)
        input_pds.append(pd)
    
# get requested number of lines per subject
input_pds_downsampled = input_pds
#input_pds_downsampled = list()
#for pd in input_pds:
#    input_pds_downsampled.append(wma.filter.downsample(pd, number_of_fibers_per_subject))

   
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
#distances = numpy.zeros([number_of_fiber_centroids, total_number_of_fibers])
#for idx in range(number_of_fiber_centroids):
#    if idx % 100 == 0:
#            print 'distances: ', idx, '/', number_of_fiber_centroids
#    fiber = fiber_array.get_fiber(fiber_sample[idx])
#    distances[idx,:] = wma.similarity.fiber_distance(fiber, fiber_array, threshold=0, distance_method='Hausdorff')

distances = Parallel(n_jobs=n_jobs, verbose=1)(
    delayed(wma.similarity.fiber_distance)(
        fiber_array.get_fiber(fiber_sample[idx]),
        fiber_array,
        0, distance_method='Hausdorff')
    for idx in range(number_of_fiber_centroids))
distances = numpy.array(distances)

        
distances_centroids = distances[:,fiber_sample]
distances_sq = numpy.multiply(distances, distances)
print "Done computing distances"





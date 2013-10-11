import whitematteranalysis as wma
import numpy
import vtk
from joblib import Parallel, delayed
import glob
import os

number_of_fibers = 10000
number_of_fibers = 3000
#number_of_fibers = 6000
#number_of_fibers = 1000
fiber_length  = 30

inputDirectory = '/Users/lauren/Data/TBI/Tracts/'
inputMask = "{0}/*.vtp".format(inputDirectory)
inputPolyDatas = glob.glob(inputMask)

# 9 to 11 are above the others. need to center?
#inputPolyDatas = inputPolyDatas[0:8] + inputPolyDatas[11:13]
inputPolyDatas = inputPolyDatas[0:8] + inputPolyDatas[11:15] + inputPolyDatas[16:17] 
#inputPolyDatas = inputPolyDatas[0:3]

print inputPolyDatas

#fname = 'input/AG_5114.vtk'
#fname = inputPolyDatas[0]

#for fname in inputPolyDatas:
for fname in [inputPolyDatas[0]]:
    pd = wma.io.read_polydata(fname)
    print 'processing'
    pd = wma.filter.preprocess(pd,fiber_length)
    print 'random downsampling'
    pd = wma.filter.downsample(pd, number_of_fibers)

    #pd2 = wma.filter.anisotropic_smooth(pd, 5, points_per_fiber=20, n_jobs=10)
    #pd2, count = wma.filter.anisotropic_smooth(pd, 2.5, points_per_fiber=20, n_jobs=3)
    #pd2, count = wma.filter.anisotropic_smooth(pd, 1, points_per_fiber=25, n_jobs=3)
    #pd2, count, pdclusters, cluster_num, cluster_count = wma.filter.anisotropic_smooth(pd, 1.5, points_per_fiber=25, n_jobs=5)
    # higher threshold for hausdorff distance
    pd2, count, pdclusters, cluster_num, cluster_count = wma.filter.anisotropic_smooth(pd, 25.0, points_per_fiber=25, n_jobs=5)
    #pd2, count = wma.filter.anisotropic_smooth(pd, 2.5, points_per_fiber=25, n_jobs=3)

    # write new one in current directory
    wma.io.write_polydata(pd2, 'smooth_' + os.path.basename(fname))
    
    # write one without junk fibers also
    pdA = wma.filter.mask(pd2, count > 3, count)
    wma.io.write_polydata(pdA, 'smooth_greater_5_' + os.path.basename(fname))

    pdB = wma.filter.mask(pdclusters, cluster_count > 3, cluster_num)
    wma.io.write_polydata(pdB, 'clusters_greater_5_' + os.path.basename(fname))

    

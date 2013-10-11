import whitematteranalysis as wma
import numpy
import vtk
from joblib import Parallel, delayed
import glob
import os

number_of_fibers = 10000
#number_of_fibers = 3000
#number_of_fibers = 6000
number_of_fibers = 1500
#number_of_fibers = 100
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

for fname in inputPolyDatas:
    #for fname in [inputPolyDatas[0]]:
    pd = wma.io.read_polydata(fname)
    print 'processing'
    pd = wma.filter.preprocess(pd,fiber_length)
    print 'random downsampling'
    pd = wma.filter.downsample(pd, number_of_fibers)

    pd2, weights = wma.filter.smooth(pd, fiber_distance_sigma = 20.0, points_per_fiber=35, n_jobs=10, upper_thresh=40)
    
    #wma.render.render(pd2)

    # write new one in current directory
    wma.io.write_polydata(pd2, 'gaussian_smooth_sig20mm_' + os.path.basename(fname))
    
    # write one without junk fibers also
    pdA = wma.filter.mask(pd2, weights >= 2, weights)
    wma.io.write_polydata(pdA, 'gaussian_smooth_sig20mm_wgr2' + os.path.basename(fname))

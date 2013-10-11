import whitematteranalysis as wma
import numpy
import vtk
from joblib import Parallel, delayed
import glob

number_of_fibers = 10000

inputDirectory = '/Users/lauren/Data/TBI/Tracts/'
inputMask = "{0}/*.vtp".format(inputDirectory)
inputPolyDatas = glob.glob(inputMask)
input_pds = list()

# 9 to 11 are above the others. need to center?
#inputPolyDatas = inputPolyDatas[0:8] + inputPolyDatas[11:13]
inputPolyDatas = inputPolyDatas[0:8] + inputPolyDatas[11:15] + inputPolyDatas[16:17] 
#inputPolyDatas = inputPolyDatas[0:3]

print inputPolyDatas

fname = 'input/AG_5114.vtk'
fname = inputPolyDatas[0]

pd = wma.io.read_polydata(fname)
pd = wma.filter.downsample(pd, number_of_fibers)

pd2 = wma.filter.anisotropic_smooth(pd, 5, points_per_fiber=20, n_jobs=10)

ren = wma.render.render(pd2)

scalars = pd2.GetCellData().GetScalars()
nfib = scalars.GetNumberOfTuples()
colors = list()
for fidx in range(0, nfib):
    colors.append(scalars.GetTuple1(fidx))

colors = numpy.array(colors)
mask = colors > 6

pd3 = wma.filter.mask(pd2, mask, colors)
ren3 = wma.render.render(pd3, scalar_bar=True)

    

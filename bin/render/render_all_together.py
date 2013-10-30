import glob
import os
import numpy
import vtk
import whitematteranalysis as wma

number_of_fibers = 1000
minimum_length = 30

#inputDirectory='/Users/lauren/Desktop/TBI-RESULTS/midsag_align'
indir = '.'

inputMask = "{0}/*.vtk".format(indir)

inputPolyDatas = glob.glob(inputMask)

input_pds = list()

def view_polydatas(polydata_list):
    appender = vtk.vtkAppendPolyData()
    idx = 0
    for pd in polydata_list:
        nf = pd.GetNumberOfLines()
        print idx
        print nf
        mask = numpy.ones(nf)
        colors = numpy.multiply(mask, idx-1)
        pd2 = wma.filter.mask(pd, mask, colors)
        appender.AddInputData(pd2)
        idx = idx + 1
    appender.Update()
    pd3 = appender.GetOutput()
    ren = wma.render.render(pd3)
    return ren

for fname in inputPolyDatas:
    print fname
    pd = wma.io.read_polydata(fname)
    pd2 = wma.filter.preprocess(pd, minimum_length)
    pd3 = wma.filter.downsample(pd2, number_of_fibers)
    input_pds.append(pd3)
    #input_pds.append(pd)

# view input data
ren = view_polydatas(input_pds)

ren.save_views()




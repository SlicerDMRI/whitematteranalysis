
import os
import glob
import matplotlib.pyplot as plt
import numpy

import vtk

import whitematteranalysis as wma

import multiprocessing
parallel_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', parallel_jobs
parallel_jobs *= 2
print 'Using two threads per CPU:', parallel_jobs



fname = 'test.vtp'
pd = wma.io.read_polydata(fname)

#pd_filtered, magnitudes = wma.filter.laplacian_of_gaussian(pd,5)
pd_filtered, magnitudes = wma.filter.laplacian_of_gaussian(pd,20)
pd_filtered, magnitudes = wma.filter.laplacian_of_gaussian(pd,10)

ren = wma.render.render(pd_filtered, scalar_bar=True, data_mode="Point")

glyph = vtk.vtkGlyph3D()
line = vtk.vtkLineSource()
glyph.SetSource(line.GetOutput())
glyph.SetInput(pd_filtered)
glyph.SetScaleModeToScaleByVector()
glyph.SetColorModeToColorByScalar()
glyph.Update()
tube = vtk.vtkTubeFilter()
tube.SetInput(glyph.GetOutput())
tube.SetNumberOfSides(12)

appender = vtk.vtkAppendPolyData()
#appender.AddInput(glyph.GetOutput())
appender.AddInput(tube.GetOutput())
appender.AddInput(pd_filtered)
appender.Update()

#ren2 = wma.render.render(glyph.GetOutput(), scalar_bar=True, tube=False, data_mode="Point")
ren2 = wma.render.render(appender.GetOutput(), scalar_bar=True, tube=False, data_mode="Point")

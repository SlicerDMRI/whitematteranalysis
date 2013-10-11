import vtk
import whitematteranalysis as wma
pd = wma.io.read_polydata('white_matter_0000.vtk')
pd2 = wma.filter.downsample(pd, 50)
pd2 = wma.filter.downsample(pd, 50, preserve_point_data=True)

pd2.GetPointData().SetActiveScalars('FA')

ren = wma.render.render(pd2, scalar_bar=True, scalar_range=[0,1], data_mode="Point")

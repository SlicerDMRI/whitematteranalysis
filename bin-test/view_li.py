import whitematteranalysis as wma
pd = wma.io.read_polydata('tractography_with_LI.vtk')
ren = wma.render.render(pd, scalar_range=[-.5,.5], scalar_bar=True)
ren.save_views()


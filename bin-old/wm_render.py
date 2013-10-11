import whitematteranalysis as wma

pd = wma.io.read_polydata('/Users/lauren/Data/TBI/Tracts/01035-dwi-filt-Ed-DTI-tract.vtp')
pd2 = wma.filter.downsample(pd,1000)
#ren=wma.render.RenderPolyData()
#ren.render(pd2)

ren = wma.render.render(pd2)

dir='/Users/lauren/Desktop/OUTPUT/'
ren.save_views(dir)


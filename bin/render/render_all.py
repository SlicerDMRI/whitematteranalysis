import glob
import os
import whitematteranalysis as wma

#inputDirectory='/Users/lauren/Desktop/TBI-RESULTS/midsag_align'
indir = '.'
outdir = 'rendered_data'
if not os.path.exists(outdir):
    os.makedirs(outdir)

inputMask = "{0}/*.vtp".format(indir)
inputMask2 = "{0}/*.vtk".format(indir)

inputPolyDatas = glob.glob(inputMask) + glob.glob(inputMask2)

idx = 1

renderers = list()

for pd_file in inputPolyDatas:
    print "S ", idx, "::", pd_file
    pd = wma.io.read_polydata(pd_file)
    ren = wma.render.render(pd,1000)
    ren.view_superior()
    # prevent immediate deletion of this renderer so we can look at them all
    renderers.append(ren)

    # view input data
    outdir_subdir = os.path.join(outdir, str(idx))
    if not os.path.exists(outdir_subdir):
        os.makedirs(outdir_subdir)
    ren.save_views(outdir_subdir)
    idx = idx + 1




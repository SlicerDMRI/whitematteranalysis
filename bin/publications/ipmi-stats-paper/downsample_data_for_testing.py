import whitematteranalysis as wma
import glob
import os

minimum_length = 40
number_of_fibers = 500

indir = '/Users/odonnell/Data/PNL-Delisi-Subset/Tractography_UKF_2t'
outdir = '/Users/odonnell/Data/PNL-Delisi-Subset/Tractography_UKF_2t_downsample'

#outdir = '.'

input_mask = "{0}/*.vtk".format(indir)
input_poly_datas = glob.glob(input_mask)


print input_poly_datas

import whitematteranalysis as wma

print 'Read and preprocess'
input_pds = list()
for fname in input_poly_datas:
    print fname
    pd = wma.io.read_polydata(fname)
    pd2 = wma.filter.preprocess(pd, minimum_length)
    pd3 = wma.filter.downsample(pd2, number_of_fibers)
    fname_out = os.path.join(outdir, os.path.basename(fname))
    print fname_out
    wma.io.write_polydata(pd3, fname_out)


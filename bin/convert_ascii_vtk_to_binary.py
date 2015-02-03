#!/usr/bin/env python
indir = '/Users/odonnell/Data/tbi_with_scalars'
outdir = '/Users/odonnell/Data/tbi_with_scalars/binary_data'

import glob
import os.path
from joblib import Parallel, delayed
import whitematteranalysis as wma

input_mask = "{0}/*.vtk".format(indir)
input_poly_datas = glob.glob(input_mask)
print input_poly_datas

def convert_pd_to_binary(fname, outdir):
    #print fname
    head, tail = os.path.split(fname)
    print tail, '---->'
    pd = wma.io.read_polydata(fname)
    fname2, ext = os.path.splitext(tail)
    fname2 += '.vtp'
    fname2 = os.path.join(outdir, fname2)
    print fname2
    wma.io.write_polydata(pd, fname2)

Parallel(n_jobs=4,
         verbose=1)(
    delayed(convert_pd_to_binary)(fname, outdir)
    for fname in input_poly_datas)


#!/usr/bin/env python


import glob
import os
import argparse
import multiprocessing
import numpy
import vtk

try:
    import whitematteranalysis as wma
except:
    print("Error importing white matter analysis package\n")
    raise

try:
    from joblib import Parallel, delayed
except:
    print("Error importing joblib package\n")
    raise


#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Make sign change of the gradient direction",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")

parser.add_argument(
    'inputnrrd',
    help='Nrrd header in .nhdr format')
parser.add_argument(
    'outputnrrd',
    help='New Nrrd header')
parser.add_argument(
    '-d', action="store", dest="dim", type=str,
    help='The dimension to change: x, y, or z.')

args = parser.parse_args()


inputnrrd = args.inputnrrd
outputnrrd = args.outputnrrd

with open(inputnrrd) as f:
    contents = f.readlines()

outstr = ''
for line in contents:
    tmpstr = line
    if line.startswith('DWMRI_gradient'):
        strs = line.split(':=')
        dirs = strs[1].split(' ')

        if args.dim == 'x':
            val = float(dirs[0])
            if val != 0:
                val = -val
                dirs[0] = str(val)
        elif args.dim == 'y':
            val = float(dirs[1])
            if val != 0:
                val = -val
                dirs[1] = str(val)
        elif args.dim == 'z':
            val = float(dirs[2])
            if val != 0:
                val = -val
                dirs[2] = str(val) + '\n'

        tmpstr = strs[0] + ':=' + dirs[0] + ' ' + dirs[1] + ' ' + dirs[2]

    outstr += tmpstr

with open(outputnrrd, 'w') as o:
    o.write(outstr)

print('Done! New nrrd header: ' + outputnrrd)


    
#!/usr/bin/env python

import argparse
import os
import numpy
import vtk
import time

try:
    import whitematteranalysis as wma
except:
    print "<wm_register.py> Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Create a MRML file that includes all vtk and vtp files in the input directory. Color tracts with different colors in the scene.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='A directory of whole-brain tractography as vtkPolyData (.vtk or .vtp). Output MRML scene file scene.mrml will be stored here and will point to all tractography files.')

args = parser.parse_args()

if not os.path.isdir(args.inputDirectory):
    print "<create_mrml> Error: Input directory", args.inputDirectory, "does not exist."
    exit()

mrml_filename = "scene.mrml"

input_polydatas = wma.io.list_vtk_files(args.inputDirectory)
number_of_files = len(input_polydatas)
print "<quality_control> Found ", number_of_files, "vtk files in input directory:", args.inputDirectory

# define R, G, B colors
# hack a colormap. 0..255 values for each
step = int(100*255.0 / (number_of_files-1))
print step, number_of_files
R = numpy.array(range(0,100*255+1, step)) / 100.0
G = numpy.abs(range(100*-127,100*128+1, step))* 2.0 / 100.0
B = numpy.array(range(100*255+1,0, -step)) / 100.0

#print len(R), len (G), len(B)
#print R
#print G
#print B

colors = list()
idx = 0
for pd in input_polydatas:
    colors.append([R[idx], G[idx],B[idx]])
    idx += 1
colors = numpy.array(colors)
print colors
wma.mrml.write(input_polydatas, colors, os.path.join(args.inputDirectory, mrml_filename), ratio=1.0)

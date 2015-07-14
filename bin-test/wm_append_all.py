#!/usr/bin/env python


import argparse
import os
import vtk

try:
    import whitematteranalysis as wma
except:
    print "<wm_register.py> Error importing white matter analysis package\n"
    raise



#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Append all polydata in a directory into one polydata",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"Unbiased Groupwise Registration of White Matter Tractography. LJ O'Donnell,  WM Wells III, Golby AJ, CF Westin. Med Image Comput Comput Assist Interv. 2012;15(Pt 3):123-30.\"",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='A file of whole-brain tractography as vtkPolyData (.vtk or .vtp).')
parser.add_argument(
    'outputFile',
    help='The output file should be a vtk or vtp format.')


args = parser.parse_args()

if not os.path.isdir(args.inputDirectory):
    print "<register> Error: Input directory", args.inputDirectory, "does not exist."
    exit()


input_pd_fnames = wma.io.list_vtk_files(args.inputDirectory)

appender = vtk.vtkAppendPolyData()

# Read in all polydata and append together into one object
for fname in input_pd_fnames:

    print "<append.py> Reading input file:", fname

    pd = wma.io.read_polydata(fname)

    print "<append.py> Input number of fibers:", pd.GetNumberOfLines()

    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        appender.AddInputData(pd)
    else:
        appender.AddInput(pd)

# ensure the data are updated
appender.Update()
 
# save the output file
wma.io.write_polydata(appender.GetOutput(), args.outputFile)





import argparse
import os
import pickle

try:
    import whitematteranalysis as wma
except:
    print "<wm_register.py> Error importing white matter analysis package\n"
    raise



#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Grabs fibers with desired indices from whole brain polydata"
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"Unbiased Groupwise Registration of White Matter Tractography. LJ O'Donnell,  WM Wells III, Golby AJ, CF Westin. Med Image Comput Comput Assist Interv. 2012;15(Pt 3):123-30.\"",
    version='1.0')

parser.add_argument(
    'inputFile',
    help='A file of whole-brain tractography as vtkPolyData (.vtk or .vtp).')
parser.add_argument(
    'inputMask',
    help='The input file (pickled) holding indices for fibers to be masked.')
parser.add_argument(
    'outputFile',
    help='The output file should be a vtk or vtp format.')

if not os.path.exists(args.inputFile):
    print "<register> Error: Input file", args.inputFile, "does not exist."
    exit()



subject_id = os.path.splitext(os.path.basename(args.inputFile))[0]

print "<mask_lines.py> Reading ", args.inputFile, "..."

pd = read_polydata(args.inputFile)

print "<mask_lines.py> Input number of fibers:", pd.GetNumberOfLines()

# read in pickle file
fiber_mask = pickle.load(args.inputMask)

# mask the fibers of interest
pd2 = wma.filter.mask(pd, fiber_mask)

# save the output
wma.io.write_polydata(pd2, args.outputFile)



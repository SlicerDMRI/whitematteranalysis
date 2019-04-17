#!/usr/bin/env python


from __future__ import print_function
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
    description="Applies preprocessing to input directory. Downsamples, removes short fibers. Preserves tensors and scalar point data along retained fibers.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='Contains whole-brain tractography as vtkPolyData file(s).')
parser.add_argument(
    'outputDirectory',
    help='The output directory should be a new empty directory. It will be created if needed.')
parser.add_argument(
    'transformFile',
    help='The AFFINE transform in MNI format (for ITK format, please use Slicer to transform).')
## parser.add_argument(
##     '-f', action="store", dest="numberOfFibers", type=int,
##     help='Number of fibers to keep from each dataset.')
## parser.add_argument(
##     '-l', action="store", dest="fiberLength", type=int,
##     help='Minimum length (in mm) of fibers to keep.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')
parser.add_argument(
    '-invert', action='store_true', dest="invert_flag",
    help='Apply the inverse of the transform.')

args = parser.parse_args()


if not os.path.isdir(args.inputDirectory):
    print("Error: Input directory", args.inputDirectory, "does not exist.")
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print("Output directory", outdir, "does not exist, creating it.")
    os.makedirs(outdir)

print("")
print("=====input directory======\n", args.inputDirectory)
print("=====output directory=====\n", args.outputDirectory)
print("==========================")

print('CPUs detected:', multiprocessing.cpu_count())
if args.numberOfJobs is not None:
    parallel_jobs = args.numberOfJobs
else:
    parallel_jobs = multiprocessing.cpu_count()
print('Using N jobs:', parallel_jobs)


print("==========================")

# =======================================================================
# Above this line is argument parsing. Below this line is the pipeline.
# =======================================================================


# Loop over input DWIs
inputPolyDatas = wma.io.list_vtk_files(args.inputDirectory)

print("<wm_preprocess.py> Input number of files: ", len(inputPolyDatas))

# Read in the transform to print its contents for the user
reader = vtk.vtkMNITransformReader()
reader.SetFileName(args.transformFile)
reader.Update()
transform = reader.GetTransform()
print(transform)
    
def pipeline(inputPolyDatas, sidx, args):
    # get subject identifier from unique input filename
    # -------------------
    #subjectID = os.path.splitext(os.path.basename(inputPolyDatas[sidx]))[0]
    fname = os.path.basename(inputPolyDatas[sidx])
    print("<wm_preprocess.py> ", sidx + 1, "/", len(inputPolyDatas))

    # read input vtk data
    # -------------------
    inpd = wma.io.read_polydata(inputPolyDatas[sidx])

    # Read in the transform
    reader = vtk.vtkMNITransformReader()
    reader.SetFileName(args.transformFile)
    reader.Update()
    transform = reader.GetTransform()

    if args.invert_flag == True:
        transform.Inverse()

    # Apply the transform
    transformer = vtk.vtkTransformPolyDataFilter()
    transformer.SetInput(inpd)
    transformer.SetTransform(transform)
    transformer.Update()
    outpd = transformer.GetOutput()
        
    # outputs
    # -------------------
    fname = os.path.join(args.outputDirectory, fname)
    try:
        print("Writing output polydata", fname, "...")
        wma.io.write_polydata(outpd, fname)
    except:
        print("Unknown exception in IO")
        raise
    del outpd
    del inpd

# loop over all inputs
Parallel(n_jobs=parallel_jobs, verbose=0)(
        delayed(pipeline)(inputPolyDatas, sidx, args)
        for sidx in range(0, len(inputPolyDatas)))

print("Launched all jobs")
exit()

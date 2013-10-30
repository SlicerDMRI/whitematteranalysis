#!/usr/bin/env python


import glob
import os
import argparse
import multiprocessing

try:
    import whitematteranalysis as wma
except:
    print "<wm_laterality.py> Error importing white matter analysis package\n"
    raise

try:
    from joblib import Parallel, delayed
except:
    print "<wm_laterality.py> Error importing joblib package\n"
    raise


#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Applies preprocessing to input directory. Downsamples, removes short fibers. Does not preserve tensors or retain scalar point data.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='Contains whole-brain tractography as vtkPolyData file(s).')
parser.add_argument(
    'outputDirectory',
    help='The output directory should be a new empty directory. It will be created if needed.')
parser.add_argument(
    'numberOfFibers', type=int,
    help='Number of fibers to keep from each dataset.')
parser.add_argument(
    "fiberLength", type=int,
    help='Minimum length (in mm) of fibers to keep.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')

args = parser.parse_args()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

print "wm_laterality. Starting white matter laterality computation."
print ""
print "=====input directory======\n", args.inputDirectory
print "=====output directory=====\n", args.outputDirectory
print "=========================="

print "fibers to retain per subject: ", args.numberOfFibers

if args.fiberLength is not None:
    print "minimum length of fibers to retain (in mm): ", args.fiberLength

print 'CPUs detected:', multiprocessing.cpu_count()
if args.numberOfJobs is not None:
    parallel_jobs = args.numberOfJobs
else:
    parallel_jobs = multiprocessing.cpu_count()
print 'Using N jobs:', parallel_jobs


print "=========================="

# =======================================================================
# Above this line is argument parsing. Below this line is the pipeline.
# =======================================================================

# Loop over input DWIs
inputMask1 = "{0}/*.vtk".format(args.inputDirectory)
inputMask2 = "{0}/*.vtp".format(args.inputDirectory)

inputPolyDatas = glob.glob(inputMask1) + glob.glob(inputMask2)

print "<wm_preprocess.py> Input number of files: ", len(inputPolyDatas)

# for testing
#inputPolyDatas = inputPolyDatas[0:2]

def pipeline(inputPolyDatas, sidx, args):
    # get subject identifier from unique input filename
    # -------------------
    subjectID = os.path.splitext(os.path.basename(inputPolyDatas[sidx]))[0]
    id_msg = "<wm_preprocess.py> ", sidx + 1, "/", len(inputPolyDatas)  
    msg = "**Starting subject:", subjectID
    print id_msg, msg

    # read input vtk data
    # -------------------
    msg = "**Reading input:", subjectID
    print id_msg, msg

    wm = wma.io.read_polydata(inputPolyDatas[sidx])

    # remove short fibers
    # -------------------
    msg = "**Preprocessing:", subjectID
    print id_msg, msg

    wm2 = wma.filter.preprocess(wm, args.fiberLength)
    del wm
    print "Number of fibers retained: ", wm2.GetNumberOfLines()
        
    # downsample 
    # -------------------
    msg = "**Downsampling input:", subjectID, " number of fibers: ", args.numberOfFibers
    print id_msg, msg

    # , preserve_point_data=True needs editing of preprocess function to use mask function
    wm3 = wma.filter.downsample(wm2, args.numberOfFibers)
    del wm2
    print "Number of fibers retained: ", wm3.GetNumberOfLines()
        
    # outputs
    # -------------------
    msg = "**Writing output data for subject:", subjectID
    print id_msg, msg

    fname = os.path.join(args.outputDirectory, subjectID+'_pp.vtp')
    try:
        print "Writing output polydata", fname, "..."
        wma.io.write_polydata(wm3, fname)
        print "Wrote output", fname, "."
    except:
        print "Unknown exception in IO"
        raise
    del wm3

# loop over all inputs
Parallel(n_jobs=parallel_jobs, verbose=0)(
        delayed(pipeline)(inputPolyDatas, sidx, args)
        for sidx in range(0, len(inputPolyDatas)))

exit()

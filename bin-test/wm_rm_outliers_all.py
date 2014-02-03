#!/usr/bin/env python


import glob
import os
import argparse
import multiprocessing

try:
    import whitematteranalysis as wma
except:
    print "<wm_outlier.py> Error importing white matter analysis package\n"
    raise

try:
    from joblib import Parallel, delayed
except:
    print "<wm_outlier.py> Error importing joblib package\n"
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
    '-d', action="store", dest="distanceThreshold", type=float,
    help='Distance threshold for outlier detection')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')

args = parser.parse_args()


if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

print "wm_outlier. Starting white matter outlier computation."
print ""
print "=====input directory======\n", args.inputDirectory
print "=====output directory=====\n", args.outputDirectory
print "=========================="

if args.distanceThreshold is not None:
    print "Distance threshold: ", args.distanceThreshold
else:
    args.distanceThreshold  = 20.0
    print "Distance threshold set to default: ", args.distanceThreshold

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
    print(id_msg + msg)

    # read input vtk data
    # -------------------
    msg = "**Reading input:", subjectID
    print(id_msg + msg)

    wm = wma.io.read_polydata(inputPolyDatas[sidx])

    num_lines = wm.GetNumberOfLines()
    print "Input number of fibers", num_lines
        
    # outlier removal 
    # -------------------
    wm2 = None
    if args.distanceThreshold is not None:
        msg = "**Downsampling input:", subjectID, " number of fibers: ", args.distanceThreshold
        print(id_msg + msg)

        # , preserve_point_data=True needs editing of preprocess function to use mask function
        keep, mask, reject = wma.filter.remove_outliers(wm, args.distanceThreshold, n_jobs=0, distance_method ='Mean')
        print "Number of fibers retained: ", keep.GetNumberOfLines(), "/", num_lines

    del wm
        
    # outputs
    # -------------------
    msg = "**Writing output data for subject:", subjectID
    print id_msg, msg

    fname = os.path.join(args.outputDirectory, subjectID+'_keep.vtp')
    try:
        print "Writing output polydata", fname, "..."
        wma.io.write_polydata(keep, fname)
        print "Wrote output", fname, "."
    except:
        print "Unknown exception in IO"
        raise
    del keep

    fname = os.path.join(args.outputDirectory, subjectID+'_reject.vtp')
    try:
        print "Writing output polydata", fname, "..."
        wma.io.write_polydata(reject, fname)
        print "Wrote output", fname, "."
    except:
        print "Unknown exception in IO"
        raise
    del reject

# loop over all inputs
Parallel(n_jobs=parallel_jobs, verbose=0)(
        delayed(pipeline)(inputPolyDatas, sidx, args)
        for sidx in range(0, len(inputPolyDatas)))

exit()

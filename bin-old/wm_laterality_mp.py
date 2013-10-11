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
    description="Applies white matter laterality pipeline to input directory.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='Contains whole-brain tractography as vtkPolyData file(s).')
parser.add_argument(
    'outputDirectory',
    help='The output directory should be a new empty directory.')
parser.add_argument(
    '-f', action="store", dest="numberOfFibers", type=int,
    help='Number of fibers to analyze from each dataset.')
parser.add_argument(
    '-l', action="store", dest="fiberLength", type=int,
    help='Minimum length (in mm) of fibers to analyze.')
parser.add_argument(
    '-s', action="store", dest="sigma", type=int,
    help='Sigma for laterality computation. Useful range 5-15 (mm).')
parser.add_argument(
    '-t', action="store", dest="threshold", type=int,
    help='Threshold lower fiber distances to 0. Useful range 0-5mm.')
parser.add_argument(
    '-no_align', action='store_false', dest="flag_midsagittalAlignment")
parser.add_argument(
    '-rm_outlier', action='store_true', dest="flag_removeOutliers")

args = parser.parse_args()

if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

if not os.path.isdir(args.outputDirectory):
    print "Error: Output directory", args.outputDirectory, "does not exist."
    exit()

print "wm_laterality. Starting white matter laterality computation."
print ""
print "=====input directory======\n", args.inputDirectory
print "=====output directory=====\n", args.outputDirectory
print "=========================="

if args.numberOfFibers is not None:
    print "fibers to analyze per subject: ", args.numberOfFibers
else:
    print "fibers to analyze per subject: ALL"

if args.fiberLength is not None:
    print "minimum length of fibers to analyze (in mm): ", args.fiberLength

if args.sigma is not None:
    print "sigma for laterality index computation: ", args.sigma

if args.threshold is not None:
    print "intra-fiber distance threshold (in mm): ", args.threshold

if args.flag_midsagittalAlignment:
    print "Automatic midsagittal alignment is ON."
else:
    print "Automatic midsagittal alignment is OFF."


if args.flag_removeOutliers:
    print "Automatic outlier removal is ON."
else:
    print "Automatic outlier removal is OFF."

print "=========================="

# =======================================================================
# Above this line is argument parsing. Below this line is the pipeline.
# =======================================================================

# Loop over input DWIs
inputMask1 = "{0}/*.vtk".format(args.inputDirectory)
inputMask2 = "{0}/*.vtp".format(args.inputDirectory)

inputPolyDatas = glob.glob(inputMask1) + glob.glob(inputMask2)

print "<wm_laterality.py> Input number of files: ", len(inputPolyDatas)
number_of_jobs = multiprocessing.cpu_count()
number_of_jobs = 10
print 'CPUs detected:', number_of_jobs


def pipeline(inputPolyDatas, sidx, args):
    # get subject identifier from unique input filename
    # -------------------
    subjectID = os.path.splitext(os.path.basename(inputPolyDatas[sidx]))[0]
    id_msg = "<wm_laterality.py> ", sidx + 1, "/", len(inputPolyDatas)  
    msg = "**Starting subject:", subjectID
    print id_msg, msg

    # read input vtk data
    # -------------------
    msg = "**Reading input:", subjectID
    print id_msg, msg

    wm = wma.io.read_polydata(inputPolyDatas[sidx])

    # remove short fibers
    # -------------------
    if args.fiberLength is not None:
        msg = "**Preprocessing:", subjectID
        print id_msg, msg

        wm = wma.filter.preprocess(wm, args.fiberLength)
        print "Number of fibers retained: ", wm.GetNumberOfLines()

    # remove outlier fibers
    # -------------------
    if args.flag_removeOutliers:
        msg = "**Removing outliers:", subjectID
        print id_msg, msg

        # if it's huge downsample to twice requested size first
        if args.numberOfFibers is not None:
            if (wm.GetNumberOfLines() > args.numberOfFibers * 2):
                wm = wma.filter.downsample(
                    wm, args.numberOfFibers * 2)
                print wm.GetNumberOfLines()

        outlierThreshold = 10
        wm = wma.filter.remove_outliers(wm, outlierThreshold)


    # do midsagittal alignment on each dataset
    # -------------------
    if args.flag_midsagittalAlignment:
        msg = "**Aligning midsag.:", subjectID
        print id_msg, msg

        align = wma.midsagalign.MidsagittalAlignment()
        align.parallel_jobs = 1
        wm, transform = align.compute(wm)
    #else:
    #    wm_align = wm

    # downsample if requested
    # -------------------
    if args.numberOfFibers is not None:
        msg = "**Downsampling input:", subjectID
        print id_msg, msg

        wm = wma.filter.downsample(wm, args.numberOfFibers)

    # compute laterality on each dataset
    # -------------------
    msg = "**Computing laterality:", subjectID
    print id_msg, msg

    laterality = wma.laterality.WhiteMatterLaterality()
    if args.sigma is not None:
        laterality.sigma = args.sigma
    if args.threshold is not None:
        laterality.threshold = args.threshold

    laterality.parallel_jobs = 1
    laterality_results = laterality.compute(wm)

    # outputs
    # -------------------
    msg = "**Writing output data for subject:", subjectID
    print id_msg, msg

    outdir = os.path.join(args.outputDirectory, subjectID)
    try:
        print "Writing output files..."
        laterality_results.write(outdir)
        print "wrote output"
    except:
        print "Unknown exception in IO"
        raise


# loop over all inputs
Parallel(n_jobs=number_of_jobs, verbose=0)(
        delayed(pipeline)(inputPolyDatas, sidx, args)
        for sidx in range(0, len(inputPolyDatas)))

exit()

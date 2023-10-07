#!/usr/bin/env python


import glob
import os
import argparse
import multiprocessing
import numpy

import whitematteranalysis as wma

from joblib import Parallel, delayed


def main():
    #-----------------
    # Parse arguments
    #-----------------
    parser = argparse.ArgumentParser(
        description="Applies preprocessing to input directory. Downsamples, removes short fibers. Preserves tensors and scalar point data along retained fibers.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu")
    
    parser.add_argument(
        'inputDirectory',
        help='Contains whole-brain tractography as vtkPolyData file(s).')
    parser.add_argument(
        'outputDirectory',
        help='The output directory should be a new empty directory. It will be created if needed.')
    ## parser.add_argument(
    ##     '-f', action="store", dest="numberOfFibers", type=int,
    ##     help='Number of fibers to keep from each dataset.')
    ## parser.add_argument(
    ##     '-l', action="store", dest="fiberLength", type=int,
    ##     help='Minimum length (in mm) of fibers to keep.')
    parser.add_argument(
        '-j', action="store", dest="numberOfJobs", type=int,
        help='Number of processors to use.')
    
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
    
    print(f"<{os.path.basename(__file__)}> Input number of files: ", len(inputPolyDatas))
    
    # for testing
    #inputPolyDatas = inputPolyDatas[0:2]
    
    def pipeline(inputPolyDatas, sidx, args):
        # get subject identifier from unique input filename
        # -------------------
        #subjectID = os.path.splitext(os.path.basename(inputPolyDatas[sidx]))[0]
        fname = os.path.basename(inputPolyDatas[sidx])
        print(f"<{os.path.basename(__file__)}> ", sidx + 1, "/", len(inputPolyDatas))
    
        # read input vtk data
        # -------------------
        inpd = wma.io.read_polydata(inputPolyDatas[sidx])
        num_lines = inpd.GetNumberOfLines()
        fiber_mask = numpy.ones(num_lines)
        outpd = wma.filter.mask(inpd, fiber_mask, preserve_point_data=False, preserve_cell_data=False, verbose=False)
        print("Number of fibers retained: ", outpd.GetNumberOfLines(), "/", num_lines)
            
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

if __name__ == '__main__':
    main()

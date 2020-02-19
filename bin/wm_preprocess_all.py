#!/usr/bin/env python


import glob
import os
import argparse
import multiprocessing

try:
    import whitematteranalysis as wma
except:
    print("<wm_laterality.py> Error importing white matter analysis package\n")
    raise

try:
    from joblib import Parallel, delayed
except:
    print("<wm_laterality.py> Error importing joblib package\n")
    raise

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
    parser.add_argument(
        '-f', action="store", dest="numberOfFibers", type=int,
        help='Number of fibers to keep from each dataset.')
    parser.add_argument(
        '-l', action="store", dest="fiberLength", type=int,
        help='Minimum length (in mm) of fibers to keep.')
    parser.add_argument(
        '-j', action="store", dest="numberOfJobs", type=int,
        help='Number of processors to use.')
    parser.add_argument(
        '-retaindata', action='store_true', dest="flag_retaindata",
        help='If given, all point and cell data stored along the tractography will be retained.')
    
    args = parser.parse_args()
    
    
    if not os.path.isdir(args.inputDirectory):
        print("Error: Input directory", args.inputDirectory, "does not exist.")
        exit()
    
    outdir = args.outputDirectory
    if not os.path.exists(outdir):
        print("Output directory", outdir, "does not exist, creating it.")
        os.makedirs(outdir)
    
    print("wm_laterality. Starting white matter laterality computation.")
    print("")
    print("=====input directory======\n", args.inputDirectory)
    print("=====output directory=====\n", args.outputDirectory)
    print("==========================")
    
    if args.numberOfFibers is not None:
        print("fibers to retain per subject: ", args.numberOfFibers)
    else:
        print("fibers to retain per subject: ALL")
    
    if args.fiberLength is not None:
        print("minimum length of fibers to retain (in mm): ", args.fiberLength)
    else:
        print("minimum length of fibers to retain (in mm): 0")
    
    print('CPUs detected:', multiprocessing.cpu_count())
    if args.numberOfJobs is not None:
        parallel_jobs = args.numberOfJobs
    else:
        parallel_jobs = multiprocessing.cpu_count()
    print('Using N jobs:', parallel_jobs)
    
    if args.flag_retaindata:
        print("Retain all data stored along the tractography.")
    else:
        print("Remove all data stored along the tractography and only keep fiber streamlines.")
    retaindata = args.flag_retaindata
    
    print("==========================")
    
    # =======================================================================
    # Above this line is argument parsing. Below this line is the pipeline.
    # =======================================================================
    
    # Loop over input DWIs
    inputPolyDatas = wma.io.list_vtk_files(args.inputDirectory)
    
    print("<wm_preprocess.py> Input number of files: ", len(inputPolyDatas))
    
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
        #print "Input number of fibers", num_lines
        
        # remove short fibers
        # -------------------
        wm2 = None
        if args.fiberLength is not None:
            msg = "**Preprocessing:", subjectID
            print(id_msg + msg)
            wm2 = wma.filter.preprocess(wm, args.fiberLength, preserve_point_data=retaindata, preserve_cell_data=retaindata, verbose=False)
            print("Number of fibers retained (length threshold", args.fiberLength, "): ", wm2.GetNumberOfLines(), "/", num_lines)
    
        if wm2 is None:
            wm2 = wm
        else:
            del wm
            
        # downsample 
        # -------------------
        wm3 = None
        if args.numberOfFibers is not None:
            msg = "**Downsampling input:", subjectID, " number of fibers: ", args.numberOfFibers
            print(id_msg + msg)
    
            # , preserve_point_data=True needs editing of preprocess function to use mask function
            wm3 = wma.filter.downsample(wm2, args.numberOfFibers, preserve_point_data=retaindata, preserve_cell_data=retaindata, verbose=False)
            print("Number of fibers retained: ", wm3.GetNumberOfLines(), "/", num_lines)
    
        if wm3 is None:
            wm3 = wm2
        else:
            del wm2
            
        # outputs
        # -------------------
        msg = "**Writing output data for subject:", subjectID
        print(id_msg, msg)
    
        fname = os.path.join(args.outputDirectory, subjectID+'_pp.vtp')
        try:
            print("Writing output polydata", fname, "...")
            wma.io.write_polydata(wm3, fname)
            print("Wrote output", fname, ".")
        except:
            print("Unknown exception in IO")
            raise
        del wm3
    
    # loop over all inputs
    Parallel(n_jobs=parallel_jobs, verbose=0)(
            delayed(pipeline)(inputPolyDatas, sidx, args)
            for sidx in range(0, len(inputPolyDatas)))
    
    exit()

if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import multiprocessing
import os
import time

import numpy as np
import scipy.optimize
import vtk

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Runs multisubject unbiased group registration of tractography.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"Unbiased Groupwise Registration of White Matter Tractography. LJ O'Donnell,  WM Wells III, Golby AJ, CF Westin. Med Image Comput Comput Assist Interv. 2012;15(Pt 3):123-30.\"")
    parser.add_argument(
        'inputDirectory',
        help='A directory of whole-brain tractography as vtkPolyData (.vtk or .vtp).')
    parser.add_argument(
        'outputDirectory',
        help='The output directory will be created if it does not exist.')
    parser.add_argument(
        '-mode', action="store", dest="mode", type=str, default="affine",
        help='The mode can be affine or nonrigid. Affine is the default. It should be run first before nonrigid.')
    parser.add_argument(
        '-f', action="store", dest="numberOfFibers", type=int, default=10000,
        help='Total number of fibers to analyze from each dataset. During registration, at each iteration fibers are randomly sampled from within this data. 10000 is the default number of total fibers.')
    parser.add_argument(
        '-l', action="store", dest="fiberLength", type=int, default=20,
        help='Minimum length (in mm) of fibers to analyze. Values of 20 to 40mm are most effective for registration.')
    parser.add_argument(
        '-lmax', action="store", dest="fiberLengthMax", type=int, default=260,
        help='Maximum length (in mm) of fibers to analyze. This parameter can be used to remove extremely long fibers that may have traversed several structures. For example, a value of 200 will avoid sampling the tail end of the fiber length distribution. The default is 260 mm, which is a safe value that should have little effect, as there are few to no fibers expected to be longer than 260mm.')
    parser.add_argument(
        '-j', action="store", dest="numberOfJobs", type=int, default=1,
        help='Number of processors to use.')
    parser.add_argument(
        '-verbose', action='store_true', dest="flag_verbose",
        help='Verbose. Run with -verbose to store more files and images of intermediate and final polydatas.')
    ## parser.add_argument(
    ##     '-pf', action="store", dest="pointsPerFiber", type=int, default=15,
    ##     help='Number of points for fiber representation during registration. The default of 15 is reasonable.')
    parser.add_argument(
        '-norender', action='store_true', dest="flag_norender",
        help='No Render. Prevents rendering of images that would require an X connection.')
    parser.add_argument(
        '-midsag_symmetric', action="store_true", dest="flag_midsag_symmetric",
        help='Register all subjects including reflected copies of input subjects, for a symmetric registration.')
    parser.add_argument(
        '-advanced_only_random_seed', action='store', dest="randomSeed", type=int,
        help='(Advanced parameter for testing only.) Set random seed for reproducible sampling in software tests.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    print("\n\n<register> =========GROUP REGISTRATION============")
    print(f"<{os.path.basename(__file__)}> Performing unbiased group registration.")
    print(f"<{os.path.basename(__file__)}> Input  directory: {args.inputDirectory}")
    print(f"<{os.path.basename(__file__)}> Output directory: {args.outputDirectory}")
    print("\n<register> ============PARAMETERS=================")
    
    mode = args.mode
    print(f"<{os.path.basename(__file__)}> Registration mode: {mode}")
    
    if not os.path.isdir(args.inputDirectory):
        print(f"<{os.path.basename(__file__)}> Error: Input directory {args.inputDirectory} does not exist.")
        exit()
    
    outdir = args.outputDirectory
    if not os.path.exists(outdir):
        print(f"<{os.path.basename(__file__)}> Output directory {outdir} does not exist, creating it.")
        os.makedirs(outdir)
    
    number_of_fibers = args.numberOfFibers
    print(f"<{os.path.basename(__file__)}> Number of fibers to analyze per subject: {number_of_fibers}")
    
    fiber_length = args.fiberLength
    print(f"<{os.path.basename(__file__)}> Minimum length of fibers to analyze (in mm): {fiber_length}")
    
    fiber_length_max = args.fiberLengthMax
    print(f"<{os.path.basename(__file__)}> Maximum  length of fibers to analyze (in mm): {fiber_length_max}")
    
    parallel_jobs = args.numberOfJobs
    print(f"<{os.path.basename(__file__)}> Number of jobs to use: {parallel_jobs}")
    
    if args.flag_verbose:
        print(f"<{os.path.basename(__file__)}> Verbose display and intermediate image saving ON.")
    else:
        print(f"<{os.path.basename(__file__)}> Verbose display and intermediate image saving OFF.")
    verbose = args.flag_verbose
    
    
    #points_per_fiber = args.pointsPerFiber
    #print(f"<{os.path.basename(__file__)}> Number of points for fiber representation: {points_per_fiber}")
    
    if args.flag_norender:
        print(f"<{os.path.basename(__file__)}> No rendering (for compute servers without X connection).")
    else:
        print(f"<{os.path.basename(__file__)}> Rendering. For intermediate image saving to check progress.")
    no_render = args.flag_norender
    
    if args.flag_midsag_symmetric:
        print(f"<{os.path.basename(__file__)}> Midsag_Symmetric registration ON.")
    else:
        print(f"<{os.path.basename(__file__)}> Midsag_Symmetric registration OFF.")
    midsag_symmetric = args.flag_midsag_symmetric
    
    if args.randomSeed is not None:
        print(f"<{os.path.basename(__file__)}> Setting random seed to: {args.randomSeed}")
    random_seed = args.randomSeed
    
    # -------------
    # SETTINGS
    # -------------
    
    if mode == "affine":
        sigma_per_scale = [20, 10, 10, 5, 5]
        iterations_per_scale=[3, 3, 3, 3, 3]
        maxfun_per_scale = [50, 80, 80, 80, 100]
        mean_brain_size_per_scale = [2000, 3000, 3000, 3500, 4000]
        subject_brain_size_per_scale = [1000, 1500, 1500, 1500, 1500]
        initial_step_per_scale = [10, 8, 5, 3, 2]
        final_step_per_scale = [2, 2, 1, 1, 1]
        # 10ppf gives approximately equivalent results to 20 ppf, so we use 10.
        points_per_fiber = 10
        nonrigid = False
        # This affine mode tested originally as cobylatest8
        # Cobyla optimizer is safer regarding scaling
        # Cobyla needs smaller steps than Powell to avoid leaving local minima on next iteration
    
    elif mode == "affine_fast_cheap":
        # For large datasets e.g. n=100
        sigma_per_scale = [20, 10, 10, 10]
        iterations_per_scale=[3, 3, 3, 3]
        maxfun_per_scale = [50, 80, 80, 80]
        mean_brain_size_per_scale = [2000, 3000, 3000, 3000]
        subject_brain_size_per_scale = [1000, 1000, 1000, 1000]
        initial_step_per_scale = [10, 8, 5, 3]
        # small final step is important for better convergence
        final_step_per_scale = [1, 1, 0.5, 0.5]
        points_per_fiber = 5
        nonrigid = False
    
    elif mode == "affine_neonate":
        # More fibers are needed for neonate registration or noisy data.
        # This mode is the same as affine, above, but with more fibers sampled.
        sigma_per_scale = [20, 10, 10, 5, 5]
        iterations_per_scale=[3, 3, 3, 3, 3]
        maxfun_per_scale = [50, 80, 80, 80, 100]
        mean_brain_size_per_scale = [4000, 4000, 4000, 4500, 6000]
        subject_brain_size_per_scale = [1500, 1750, 2000, 2000, 2500]
        initial_step_per_scale = [10, 8, 5, 3, 2]
        final_step_per_scale = [2, 2, 1, 1, 1]
        points_per_fiber = 10
        nonrigid = False
    
    elif mode == "nonrigid_neonate":
        # This mode uses the maximum data possible to register
        # smaller, lower-quality or noisier datasets.
        # The grid goes up to 6x6x6
        grid_resolution_per_scale = [4, 5, 6, 6]
        # this is in mm space.
        initial_step_per_scale = [5, 4, 3, 2]
        final_step_per_scale = [3, 3, 2, 1]
        # use only very local information (small sigma)
        sigma_per_scale = [5, 3, 2, 1.5]
        # how many times to repeat the process at each scale
        iterations_per_scale = [3, 3, 2, 2]
        # higher totals for variable neonate tractography
        mean_brain_size_per_scale = [4000, 4500, 5000, 6000]
        subject_brain_size_per_scale = [2000, 2500, 3000, 3500]
        # 3x3x3 grid, 27*3 = 81 parameter space.
        # 4x4x4 grid, 64*3 = 192 parameter space.
        # 5x5x5 grid, 125*3 = 375 parameter space.
        # 6x6x6 grid, 216*3 = 648 parameter space.
        # Inspection of output pdfs shows that objective decreases steadily for all subjects,
        # so stop the optimizer early and create a better current model.
        maxfun_per_scale = [205, 390, 670, 670]
        # fiber representation for computation.
        #points_per_fiber = 5
        points_per_fiber = 10
        nonrigid = True
    
    elif mode == "nonrigid":
        # This mode uses the minimum data possible to register
        # large, high-quality datasets.
        # The grid goes up to 8x8x8
        grid_resolution_per_scale = [4, 5, 6, 7, 8]
        # this is in mm space.
        initial_step_per_scale = [5, 4, 3, 2, 1.5]
        final_step_per_scale = [3, 3, 2, 1, 1]
        # use only very local information (small sigma)
        # sigma 1.25 is apparently not useful: stick with a minimum of 2mm
        sigma_per_scale = [5, 3, 2, 2, 2]
        # how many times to repeat the process at each scale
        iterations_per_scale = [10, 10, 8, 5, 2]
        # These lower than the totals in the affine because
        # this computation is more expensive
        mean_brain_size_per_scale = [1000, 1000, 1000, 1000, 1000]
        subject_brain_size_per_scale = [500, 500, 500, 500, 500]
        # 3x3x3 grid, 27*3 = 81 parameter space.
        # 4x4x4 grid, 64*3 = 192 parameter space.
        # 5x5x5 grid, 125*3 = 375 parameter space.
        # 6x6x6 grid, 216*3 = 648 parameter space.
        # 7x7x7 grid, 343*3 = 1029 parameter space.
        # 8x8x8 grid, 512*3 = 1536 parameter space.
        # 9x9x9 grid, 721*3 = 2187 parameter space.
        # 10x10x10, 1000*3 = 3000 parameter space.
        # Inspection of output pdfs shows that objective decreases steadily for all subjects,
        # so stop the optimizer early and create a better current model.
        maxfun_per_scale = [205, 390, 670, 1049, 1556]
        # fiber representation for computation.
        points_per_fiber = 3
        nonrigid = True
    
    elif mode == "nonrigid_lowsample":
        # This mode uses the minimum data possible to register
        # large, high-quality datasets. e.g. n=100
        # The grid goes up to 8x8x8
        grid_resolution_per_scale = [4, 5, 6, 7, 8]
        # this is in mm space.
        initial_step_per_scale = [5, 4, 3, 2, 1.5]
        final_step_per_scale = [3, 3, 2, 1, 1]
        # use only very local information (small sigma)
        # sigma 1.25 is apparently not useful: stick with a minimum of 2mm
        sigma_per_scale = [5, 3, 2, 2, 2]
        # how many times to repeat the process at each scale
        iterations_per_scale = [10, 10, 8, 5, 2]
        # These are lower than the totals in the affine because
        # this computation is more expensive
        mean_brain_size_per_scale = [800, 800, 800, 800, 800]
        subject_brain_size_per_scale = [450, 450, 450, 450, 450]
        # 3x3x3 grid, 27*3 = 81 parameter space.
        # 4x4x4 grid, 64*3 = 192 parameter space.
        # 5x5x5 grid, 125*3 = 375 parameter space.
        # 6x6x6 grid, 216*3 = 648 parameter space.
        # 7x7x7 grid, 343*3 = 1029 parameter space.
        # 8x8x8 grid, 512*3 = 1536 parameter space.
        # 9x9x9 grid, 721*3 = 2187 parameter space.
        # 10x10x10, 1000*3 = 3000 parameter space.
        # Inspection of output pdfs shows that objective decreases steadily for all subjects,
        # so stop the optimizer early and create a better current model.
        maxfun_per_scale = [205, 390, 670, 1049, 1556]
        # fiber representation for computation.
        points_per_fiber = 3
        nonrigid = True
    
    elif mode == "affineTEST":
        # very quick test if software is working
        sigma_per_scale = [30, 10, 7.5]
        iterations_per_scale=[1, 1, 1]
        maxfun_per_scale = [60, 80, 100]
        mean_brain_size_per_scale = [1500, 2000, 3000]
        subject_brain_size_per_scale = [100, 500, 1000]
        initial_step_per_scale = [5, 5, 5, 5]
        final_step_per_scale = [2, 2, 2, 2]
        nonrigid = False
        points_per_fiber = 5
        
    elif mode == "nonrigidTEST":
        # very quick test if software is working
        grid_resolution_per_scale = [3, 4, 5, 6, 8, 10]
        initial_step_per_scale = [5, 3, 1, 1, 1, 1]
        final_step_per_scale = [2, 1, 1, 1, 1, 1]
        sigma_per_scale = [3, 2, 1, 1, 1, 1]
        iterations_per_scale = [1, 1, 1, 1, 1, 1]
        mean_brain_size_per_scale = [1000, 1000, 1000, 1000, 1000, 1000]
        subject_brain_size_per_scale = [100, 100, 100, 100, 100, 100]
        # stop computation: this is just a quick test the software is working
        maxfun_per_scale = [10, 10, 10, 10, 10, 10]
        points_per_fiber = 15
        nonrigid = True
    
    else:
        print(f"\n<register> Error: Unknown registration mode: {mode}")
        exit()
    
    # -------------
    # REGISTRATION
    # -------------
    
    # Test the input files exist
    input_polydatas = wma.io.list_vtk_files(args.inputDirectory)
    number_of_subjects = len(input_polydatas)
    print(f"<{os.path.basename(__file__)}> Found {number_of_subjects} subjects in input directory {args.inputDirectory}")
    if number_of_subjects < 1:
        print("\n<register> Error: No .vtk or .vtp files were found in the input directory.\n")
        exit()
    
    # Get input data
    input_pds, subject_ids = wma.io.read_and_preprocess_polydata_directory(args.inputDirectory, fiber_length, number_of_fibers, random_seed, fiber_length_max)
    
    # If we are registering for symmetry, include reflected copy of each brain
    if midsag_symmetric:
        input_pds2 = []
        subject_ids2 = []
        for pd, id in zip(input_pds, subject_ids):
            # make the reflected copy
            trans = vtk.vtkTransform()
            trans.Scale(-1,1,1)
            transformer = vtk.vtkTransformPolyDataFilter()
            transformer.SetTransform(trans)
            if (vtk.vtkVersion().GetVTKMajorVersion() > 5.0):
                transformer.SetInputData(pd)
            else:
                transformer.SetInput(pd)
            transformer.Update()
            pd_reflect = transformer.GetOutput()
            # append input and its reflection to list of data to register
            input_pds2.append(pd)
            input_pds2.append(pd_reflect)
            # create subject id for reflected
            subject_ids2.append(id)
            subject_ids2.append(str(id) + '_reflect')
        input_pds = input_pds2
        subject_ids = subject_ids2
    
    print("\n<register> Starting registration...\n")
    register = wma.congeal_multisubject.MultiSubjectRegistration()
    register.input_directory = args.inputDirectory
    register.output_directory = args.outputDirectory
    register.verbose = verbose
    register.parallel_jobs = parallel_jobs
    register.render = not no_render
    if nonrigid:
        register.mode = "Nonrigid"
    # We have to add polydatas after setting nonrigid in the register object
    for (pd, id) in zip(input_pds, subject_ids):
        register.add_polydata(pd, id)
    print(f"<{os.path.basename(__file__)}> Number of points for fiber representation: {points_per_fiber}")
    register.points_per_fiber = points_per_fiber
    
    # output summary file to save information about what was run
    readme_fname = os.path.join(args.outputDirectory, 'README.txt')
    readme_file = open(readme_fname, 'w')
    outstr = "Groupwise Registration Summary\n"
    outstr += '----------------------\n'
    outstr += '\n'
    outstr += "Input Directory: "
    outstr += args.inputDirectory
    outstr += '\n'
    outstr += "Output Directory: "
    outstr += args.outputDirectory
    outstr += '\n'
    outstr += "Number of Subjects: "
    outstr += str(number_of_subjects)
    outstr += '\n'
    outstr += '\n'
    outstr +=  f"Current date: {time.strftime('%x')}"
    outstr += '\n'
    outstr +=  f"Current time: {time.strftime('%X')}"
    outstr += '\n'
    outstr += '\n'
    outstr += f"Path to Script: {os.path.realpath(__file__)}"
    outstr += '\n'
    outstr += f"Working Directory: {os.getcwd()}"
    outstr += '\n'
    outstr += '\n'
    outstr += "Description of Outputs\n"
    outstr += '---------------------\n'
    outstr += 'registration_atlas.vtk: This is the template created from groupwise registration.\n'
    outstr += 'registration_performance.txt: Parameters and objective values from each iteration.\n'
    outstr += 'progress.txt: Update of how far registration has progressed.\n'
    outstr += 'input_subjects.txt:  List of subject index, ID, and full path to input file.\n'
    outstr += 'README.txt:  This summary file.\n'
    outstr += 'output_tractography/\n'
    outstr += '\tThe output directory with registered tractography and corresponding transforms.\n'
    outstr += 'The files inside each iteration directory are for testing purposes:\n'
    outstr += '\titk_txform files output from that iteration.\n'
    outstr += '\tpdf files plot objective function changes for all subjects.\n'
    outstr += '\n'
    outstr += '\n'
    outstr += "Command Line Arguments\n"
    outstr += '----------------------\n'
    outstr += str(args)
    outstr += '\n'
    outstr += '\n'
    outstr += "Parameters\n"
    outstr += '----------------------\n'
    outstr += f"Registration mode: {mode}"
    outstr += '\n'
    outstr += f"Number of fibers to analyze per subject: {str(number_of_fibers)}"
    outstr += '\n'
    outstr += f"Minimum length of fibers to analyze (in mm): {str(fiber_length)}"
    outstr += '\n'
    outstr += f"Maximum  length of fibers to analyze (in mm): {str(fiber_length_max)}"
    outstr += '\n'
    outstr += f"Number of jobs to use: {str(parallel_jobs)}"
    outstr += '\n'
    outstr += f"verbose: {str(verbose)}"
    outstr += '\n'
    outstr += f"render: {str(not no_render)}"
    outstr += '\n'
    outstr += f"midsag_symmetric: {str(midsag_symmetric)}"
    outstr += '\n'
    outstr += f"random seed: {str(random_seed)}"
    outstr += '\n'
    outstr += '\n'
    outstr += "Input Fiber Files\n"
    outstr += '-----------------\n'
    for pd in input_polydatas:
        outstr += pd
        outstr += '\n'
    readme_file.write(outstr)
    readme_file.close()
    
    # output summary file to save information about all subjects
    subjects_qc_fname = os.path.join(args.outputDirectory, 'input_subjects.txt')
    subjects_qc_file = open(subjects_qc_fname, 'w')
    outstr = "Subject_idx\tSubject_ID\tInput Filename\n"
    subjects_qc_file.write(outstr)
    idx = 1
    for fname in input_polydatas:
        subject_id = os.path.splitext(os.path.basename(fname))[0]
        outstr =  str(idx) + '\t' + str(subject_id) + '\t' + str(fname) + '\n'
        subjects_qc_file.write(outstr)
        idx += 1
    subjects_qc_file.close()
    
    # -------------
    # Done SETTINGS. Below is computation
    # -------------
    total_iterations = np.sum(np.array(iterations_per_scale))
    iteration = 1
    # estimate percentage complete based on number of fibers compared,
    # because the times cobyla calls the objective function are approx
    # constant per scale (except first scale where they are cut short)
    total_comparisons = np.multiply(iterations_per_scale,np.multiply(np.array(mean_brain_size_per_scale), np.array(subject_brain_size_per_scale)))
    total_comparisons = np.sum(total_comparisons)
    comparisons_so_far = 0
    progress_filename = os.path.join(args.outputDirectory, 'progress.txt')
    progress_file = open(progress_filename, 'w')
    print("Beginning registration. Total iterations will be:", total_iterations, file=progress_file)
    print(f"Start date: {time.strftime('%x')}", file=progress_file)
    print(f"Start time: {time.strftime('%X')}\n", file=progress_file)
    progress_file.close()
    prev_time = time.time()
    do_scales = list(range(len(sigma_per_scale)))
    
    for scale in do_scales:
        register.sigma = sigma_per_scale[scale]
        register.initial_step = initial_step_per_scale[scale]
        register.final_step = final_step_per_scale[scale]
        register.maxfun = maxfun_per_scale[scale]
        register.mean_brain_size = mean_brain_size_per_scale[scale]
        register.subject_brain_size = subject_brain_size_per_scale[scale]
        if register.mode == "Nonrigid":
            register.nonrigid_grid_resolution = grid_resolution_per_scale[scale]
            register.update_nonrigid_grid()
        
        for idx in range(0,iterations_per_scale[scale]):
            register.iterate()
            comparisons_this_scale = mean_brain_size_per_scale[scale]*subject_brain_size_per_scale[scale]
            comparisons_so_far += comparisons_this_scale
            percent = 100*(float(comparisons_so_far)/total_comparisons)
            print(f"Done iteration {iteration} / {total_iterations}. Percent finished approx: {percent:.2f}")
            progress_file = open(progress_filename, 'a')
            curr_time = time.time()
            print(f"Done iteration {iteration} / {total_iterations}. Percent finished approx: {percent:.2f}. Time: {time.strftime('%X')}. Minutes Elapsed: {(curr_time - prev_time)/60}", file=progress_file)
            progress_file.close()
            prev_time = curr_time
    
            iteration += 1
            # Intermediate save. For testing only.
            if verbose:
                register.save_transformed_polydatas(intermediate_save=True, midsag_symmetric=midsag_symmetric)
    
    # Final save when we are done
    register.save_transformed_polydatas(midsag_symmetric=midsag_symmetric)
    
    print(f"\nDone registering. For more information on the output, please read: {readme_fname}\n")
    
    progress_file = open(progress_filename, 'a')
    print("\nFinished registration.", file=progress_file)
    print(f"End date: {time.strftime('%x')}", file=progress_file)
    print(f"End time: {time.strftime('%X')}", file=progress_file)
    progress_file.close()
    
    
    # Documentation of tested settings with other optimizers
    ## # This is the fast default mode using Powell's method for optimization
    ## if mode == "affine":
    ##     sigma_per_scale = [20, 10, 10, 5]
    ##     # Powell seems not to pay much attention to requested max.
    ##     # If sigma > 10 uses Cobyla, which does use the maxfun
    ##     maxfun_per_scale = [50, 20, 40, 80]
    ##     #mean_brain_size_per_scale = [1000, 3000, 3000, 3500]
    ##     #subject_brain_size_per_scale = [500, 1500, 1500, 1500]
    ##     mean_brain_size_per_scale = [2000, 3000, 3000, 3500]
    ##     subject_brain_size_per_scale = [1000, 1500, 1500, 1500]
    ##     initial_step_per_scale = [10, 10, 5, 3]
    ##     final_step_per_scale = [8, 8, 4, 2]
    ##     points_per_fiber = 10
    ##     register.nonrigid = False

if __name__ == '__main__':
    main()

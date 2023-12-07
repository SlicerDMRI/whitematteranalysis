#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import multiprocessing
import os
import time

import numpy as np
import vtk

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Registers a whole-brain vtk tractography file to another vtk tractography file (an atlas).",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"Unbiased Groupwise Registration of White Matter Tractography. LJ O'Donnell,  WM Wells III, Golby AJ, CF Westin. Med Image Comput Comput Assist Interv. 2012;15(Pt 3):123-30.\"")
    parser.add_argument(
        'inputSubject',
        help='One subject data: whole-brain tractography as vtkPolyData (.vtk or .vtp).')
    parser.add_argument(
        'inputAtlas',
        help='An atlas, one file containing whole-brain tractography as vtkPolyData (.vtk or .vtp).')
    parser.add_argument(
        'outputDirectory',
        help='The output directory will be created if it does not exist.')
    parser.add_argument(
        '-mode', action="store", dest="mode", type=str, default="affine",
        help='The mode can be affine or nonrigid. Affine is the default. It should be run first before nonrigid.')
    parser.add_argument(
        '-l', action="store", dest="fiberLength", type=int, default=40,
        help='Minimum length (in mm) of fibers to analyze. The default is 40mm.')
    parser.add_argument(
        '-lmax', action="store", dest="fiberLengthMax", type=int, default=260,
        help='Maximum length (in mm) of fibers to analyze. This parameter can be used to remove extremely long fibers that may have traversed several structures. For example, a value of 200 will avoid sampling the tail end of the fiber length distribution. The default is 260 mm.')
    parser.add_argument(
        '-verbose', action='store_true', dest="flag_verbose",
        help='Verbose. Run with -verbose to store more files and images of intermediate and final polydatas.')
    #parser.add_argument(
    #    '-pf', action="store", dest="pointsPerFiber", type=int, default=15,
    #    help='Number of points for fiber representation during registration. The default of 15 is reasonable.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    print("\n\n<register> =========GROUP REGISTRATION============")
    print(f"<{os.path.basename(__file__)}> Registering to atlas.")
    print(f"<{os.path.basename(__file__)}> Input  subject file: ", args.inputSubject)
    print(f"<{os.path.basename(__file__)}> Input  atlas file: ", args.inputAtlas)
    print(f"<{os.path.basename(__file__)}> Output directory: ", args.outputDirectory)
    print("\n<register> ============PARAMETERS=================")
    
    mode = args.mode
    print(f"<{os.path.basename(__file__)}> Registration mode:", mode)
    
    if not os.path.isfile(args.inputSubject):
        print(f"<{os.path.basename(__file__)}> Error: Input subject data", args.inputSubject, "does not exist.")
        exit()
    
    if not os.path.isfile(args.inputAtlas):
        print(f"<{os.path.basename(__file__)}> Error: Input atlas", args.inputAtlas, "does not exist.")
        exit()
    
    fname = args.inputSubject
    subject_id = os.path.splitext(os.path.basename(fname))[0]
    subject_pd = wma.io.read_polydata(fname)
    fname = args.inputAtlas
    atlas_id = os.path.splitext(os.path.basename(fname))[0]
    atlas_pd = wma.io.read_polydata(fname)
    
    outdir = args.outputDirectory
    if not os.path.exists(outdir):
        print(f"<{os.path.basename(__file__)}> Output directory {outdir} does not exist, creating it.")
        os.makedirs(outdir)
    subject_outdir = os.path.join(outdir, subject_id)
    if not os.path.exists(subject_outdir):
        print(f"<{os.path.basename(__file__)}> Output directory {outdir} does not exist, creating it.")
        os.makedirs(subject_outdir)
    
    fiber_length = args.fiberLength
    print(f"<{os.path.basename(__file__)}> Minimum length of fibers to analyze (in mm): ", fiber_length)
    
    fiber_length_max = args.fiberLengthMax
    print(f"<{os.path.basename(__file__)}> Maximum  length of fibers to analyze (in mm): ", fiber_length_max)
    
    if args.flag_verbose:
        print(f"<{os.path.basename(__file__)}> Verbose display and intermediate image saving ON.")
    else:
        print(f"<{os.path.basename(__file__)}> Verbose display and intermediate image saving OFF.")
    verbose = args.flag_verbose
    
    print("\n<register> Starting registration...\n")
    
    
    # -------------
    # SETTINGS
    # -------------
    
    rigid = False
    
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
        # Cobyla optimizer is safer regarding scaling
        # Cobyla needs smaller steps than Powell to avoid leaving local minima on next iteration
    
    elif mode == "affine_fast":
        # This mode is equivalent to previous parameters
        sigma_per_scale = [30, 10, 7.5, 5]
        iterations_per_scale=[2, 2, 2, 2]
        maxfun_per_scale = [45, 60, 75, 90]
        mean_brain_size_per_scale = [1000, 3000, 4000, 5000]
        subject_brain_size_per_scale = [250, 1500, 1750, 2000]
        initial_step_per_scale = [10, 5, 5, 5]
        final_step_per_scale = [5, 2, 2, 2]
        # 10ppf gives approximately equivalent results to 20 ppf, so we use 10.
        points_per_fiber = 10
        nonrigid = False
    
    elif mode == "rigid_affine_fast":
        # This mode is equivalent to previous parameters and includes an initial rigid step
        sigma_per_scale = [30, 30, 10, 7.5, 5]
        iterations_per_scale=[2, 2, 2, 2, 2]
        maxfun_per_scale = [45, 45, 60, 75, 90]
        mean_brain_size_per_scale = [1000, 1000, 3000, 4000, 5000]
        subject_brain_size_per_scale = [250, 250, 1500, 1750, 2000]
        initial_step_per_scale = [10, 10, 5, 5, 5]
        final_step_per_scale = [5, 5, 2, 2, 2]
        # 10ppf gives approximately equivalent results to 20 ppf, so we use 10.
        points_per_fiber = 10
        nonrigid = False
        rigid = True
        rigid_scale = [1, 1, 1, 0, 0]
    
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
        # This is the mode formerly known as "superfast"
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
    
    elif mode == "nonrigid_10":
        # This runs longer than default nonrigid, going to a 10x10x10 grid.
        # This mode uses the minimum data possible to register
        # large, high-quality datasets.
        # This is the mode formerly known as "superfast"
        # The grid goes up to 10x10x10
        grid_resolution_per_scale = [4, 5, 6, 7, 8, 9, 10]
        # this is in mm space.
        initial_step_per_scale = [5, 4, 3, 2, 1.5, 1.5, 1.5]
        final_step_per_scale = [3, 3, 2, 1, 1, 1, 1]
        # use only very local information (small sigma)
        # sigma 1.25 is apparently not useful: stick with a minimum of 2mm
        sigma_per_scale = [5, 3, 2, 2, 2, 2, 2]
        # how many times to repeat the process at each scale
        iterations_per_scale = [10, 10, 8, 5, 4, 3, 2]
        # These lower than the totals in the affine because
        # this computation is more expensive
        mean_brain_size_per_scale = [1000, 1000, 1000, 1000, 1000, 1000, 1000]
        subject_brain_size_per_scale = [500, 500, 500, 500, 500, 500, 500]
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
        maxfun_per_scale = [205, 390, 670, 1049, 1556, 2200, 3020]
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
    
    
    
    register = wma.congeal_to_atlas.SubjectToAtlasRegistration()
    register.output_directory = subject_outdir
    register.input_polydata_filename = args.inputSubject
    register.fiber_length = args.fiberLength
    register.fiber_length_max = args.fiberLengthMax

    if nonrigid:
        register.mode = "Nonrigid"
    # We have to add polydatas after setting nonrigid in the register object
    register.set_subject(subject_pd, subject_id)
    register.set_atlas(atlas_pd, atlas_id)
    
    register.points_per_fiber = points_per_fiber
    
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
    
    progress_filename = os.path.join(subject_outdir, 'progress.txt')
    progress_file = open(progress_filename, 'w')
    print(f"Beginning registration. Total iterations will be: {total_iterations}", file=progress_file)
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
        if rigid:
            if rigid_scale[scale]:
                register.mode = "Rigid"
            else:
                register.mode = "Affine"
    
        for idx in range(0,iterations_per_scale[scale]):
            register.iterate()
            comparisons_this_scale = mean_brain_size_per_scale[scale]*subject_brain_size_per_scale[scale]
            comparisons_so_far += comparisons_this_scale
            percent = 100*(float(comparisons_so_far)/total_comparisons)
            print(f"Done iteration {iteration} / {total_iterations}. Percent finished approx: {percent:.2f}")
            progress_file = open(progress_filename, 'a')
            curr_time = time.time()
            print(f"Done iteration {iteration} / {total_iterations}. Percent finished approx: {percent:.2f}. Time: {time.strftime('%X')}. Minutes Elapsed: {(curr_time - prev_time) / 60}", file=progress_file)
            progress_file.close()
            prev_time = curr_time
            iteration += 1
            # Intermediate save. For testing only.
            if verbose:
                register.save_transformed_polydata(intermediate_save=True)
    
    # Final save when we are done
    register.save_transformed_polydata()
    
    print(f"Done registering. See output in: {subject_outdir}")
    
    progress_file = open(progress_filename, 'a')
    print("\nFinished registration.", file=progress_file)
    print(f"End date: {time.strftime('%x')}", file=progress_file)
    print(f"End time: {time.strftime('%X')}", file=progress_file)
    progress_file.close()

if __name__ == '__main__':
    main()

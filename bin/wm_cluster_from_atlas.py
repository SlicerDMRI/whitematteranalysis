#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import multiprocessing
import os
import time

import numpy as np
import vtk

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Clusters tractography (propagates clusters) from a cluster atlas (a multi-subject/multi-atlas cluster representation).",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
    parser.add_argument(
        'inputFile',
        help='A file of whole-brain tractography as vtkPolyData (.vtk or .vtp).')
    parser.add_argument(
        'atlasDirectory',
        help='The directory where the atlas is stored. Must contain atlas.p and atlas.vtp')
    parser.add_argument(
        'outputDirectory',
        help='The output directory will be created if it does not exist.')
    parser.add_argument(
        '-f', action="store", dest="numberOfFibers", type=int,
        help='Number of fibers to analyze from each subject. If this parameter is not used, all fibers will be analyzed by default.')
    parser.add_argument(
        '-l', action="store", dest="fiberLength", type=int, default=40,
        help='Minimum length (in mm) of fibers to analyze. 60mm is default.')
    parser.add_argument(
        '-j', action="store", dest="numberOfJobs", type=str, default='1',
        help='Number of processors to use.')
    parser.add_argument(
        '-verbose', action='store_true', dest="flag_verbose",
        help='Verbose. Run with -verbose for more text output.')
    parser.add_argument(
        '-mrml_fibers', action="store", dest="showNFibersInSlicer", type=float,
        help='Approximate upper limit on number of fibers to show when MRML scene of clusters is loaded into slicer.Default is 10000 fibers; increase for computers with more memory. Note this can be edited later in the MRML file by searching for SubsamplingRatio and editing that number throughout the file. Be sure to use a text editor program (save as plain text format). An extra MRML file will be saved for visualizing 100%% of fibers.')
    parser.add_argument(
        '-reg', action='store_true', dest="registerAtlasToSubjectSpace",
        help='To cluster in individual subject space, register atlas polydata to subject. Otherwise, by default this code assumes the subject has already been registered to the atlas.')
    parser.add_argument(
        '-norender', action='store_true', dest="flag_norender",
        help='No Render. Prevents rendering of images that would require an X connection.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    if not os.path.exists(args.inputFile):
        print(f"<{os.path.basename(__file__)}> Error: Input file", args.inputFile, "does not exist.")
        exit()
    
    if not os.path.isdir(args.atlasDirectory):
        print(f"<{os.path.basename(__file__)}> Error: Atlas directory {args.atlasDirectory} does not exist.")
        exit()
    
    outdir = args.outputDirectory
    if not os.path.exists(outdir):
        print(f"<{os.path.basename(__file__)}> Output directory {outdir} does not exist, creating it.")
        os.makedirs(outdir)
    
    fname = args.inputFile
    subject_id = os.path.splitext(os.path.basename(fname))[0]
    outdir = os.path.join(outdir, subject_id)
    if not os.path.exists(outdir):
        print(f"<{os.path.basename(__file__)}> Output directory {outdir} does not exist, creating it.")
        os.makedirs(outdir)
    
    print("\n==========================")
    print("input file:", args.inputFile)
    print(f"atlas directory {args.atlasDirectory}")
    print("output directory:", args.outputDirectory)
    
    if args.numberOfFibers is not None:
        print("fibers to analyze per subject: ", args.numberOfFibers)
    else:
        print("fibers to analyze per subject: ALL")
    number_of_fibers = args.numberOfFibers
    
    fiber_length = args.fiberLength
    print("minimum length of fibers to analyze (in mm): ", fiber_length)
    
    number_of_jobs = int(args.numberOfJobs)
    print(f'Using N jobs: {number_of_jobs}')
    
    if args.flag_verbose:
        print("Verbose ON.")
    else:
        print("Verbose OFF.")
    verbose = args.flag_verbose
    
    if args.showNFibersInSlicer is not None:
        show_fibers = args.showNFibersInSlicer
    else:
        show_fibers = 10000.0
    print("Maximum total number of fibers to display in MRML/Slicer: ", show_fibers)
    
    if args.registerAtlasToSubjectSpace:
        print("Registration of atlas fibers to subject fibers is ON.")
        print("Warning: the registration pipeline is under improvements--the registration implemented here is not currently supported. Please register to atlas first, then call this script without the -reg option.")
        exit()
    else:
        print("Registration of atlas fibers to subject fibers is OFF. Subject must be in atlas space before calling this script.")
    
    if args.flag_norender:
        print("No rendering (for compute servers without X connection).")
    else:
        print("Rendering. After clustering, will create colorful jpg images.")
    render = not args.flag_norender
    
    print("==========================\n")
    
    # read atlas
    print(f"<{os.path.basename(__file__)}> Loading input atlas: {args.atlasDirectory}")
    atlas = wma.cluster.load_atlas(args.atlasDirectory, 'atlas')
    
    # read data
    print(f"<{os.path.basename(__file__)}> Reading input file:", args.inputFile)
    pd = wma.io.read_polydata(args.inputFile)
        
    # preprocessing step: minimum length
    print(f"<{os.path.basename(__file__)}> Preprocessing by length:", fiber_length, "mm.")
    pd2 = wma.filter.preprocess(pd, fiber_length, return_indices=False, preserve_point_data=True, preserve_cell_data=True,verbose=False)
    
    # preprocessing step: fibers to analyze
    if number_of_fibers is not None:
        print(f"<{os.path.basename(__file__)}> Downsampling to ", number_of_fibers, "fibers.")
        input_data = wma.filter.downsample(pd2, number_of_fibers, return_indices=False, preserve_point_data=True, preserve_cell_data=True,verbose=False)
    else:
        input_data = pd2
    
    #-----------------
    # Register if needed
    #-----------------
    
    if args.registerAtlasToSubjectSpace:
    
        if 0:
            # test: transform input data to test this works
            translation = vtk.vtkTransform()
            translation.Translate(100.0, 50.0, 22.0)
            transformer = vtk.vtkTransformPolyDataFilter()
            if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                transformer.SetInputData(input_data)
            else:
                transformer.SetInput(input_data)
            transformer.SetTransform(translation)
            transformer.Update()
            input_data = transformer.GetOutput()
    
        # defaults for two subject registration
        # (may be added as parameters later)
        fiber_sample_fractions = [.4, .5, .6, .8]
        #sigma_per_scale = [30, 10, 10, 5]
        sigma_per_scale = [10, 10, 10, 5]
        #steps_per_scale=[10, 3, 2, 2]
        steps_per_scale=[5, 2, 2, 2]
        fibers_rendered = 100
        # figure out maximum function evals for optimizer if not requested
        # default for two dataset registration
        maxfun_per_scale = [20, 40, 60, 80]
        number_of_fibers = 300
        fiber_length = 75
        points_per_fiber = 5
        distance_method='Hausdorff'
        # figure out numbers of fibers to sample
        fiber_sample_sizes = (number_of_fibers * np.array(fiber_sample_fractions)).astype(int)
        # create registration object and apply settings
        register = wma.congeal.CongealTractography()
        register.parallel_jobs = number_of_jobs
        register.threshold = 0
        register.points_per_fiber = points_per_fiber
        register.distance_method = distance_method
        register.add_subject(atlas.nystrom_polydata)
        register.add_subject(input_data)
        scales = ["Coarse", "Medium", "Fine", "Finest"]
        scale_idx = 0
        elapsed = list()
        for scale in scales:
            start = time.time()
            # run the basic iteration of translate, rotate, scale
            wma.registration_functions.compute_multiscale_registration(register, scale, steps_per_scale[scale_idx], fiber_sample_sizes[scale_idx], sigma_per_scale[scale_idx], maxfun_per_scale[scale_idx])
            elapsed.append(time.time() - start)
            scale_idx += 1
            print("Registration Time:", elapsed)
        #  NEED TO OUTPUT THE COMBINED TRANSFORM APPLIED TO ATLAS NOT TO PD
        # now apply the appropriate transform to the input data.
        # transform 0 times transform 1 inverse
        tx = register.convert_transforms_to_vtk()
        tx[1].Inverse()
        tx[0].Concatenate(tx[1])
        # apply this transform to the atlas polydata
        transformer = vtk.vtkTransformPolyDataFilter()
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            transformer.SetInputData(atlas.nystrom_polydata)
        else:
            transformer.SetInput(atlas.nystrom_polydata)
        transformer.SetTransform(tx[0])
        print(tx[0])
        transformer.Update()
        pd = transformer.GetOutput()
    
        # update the atlas polydata with the new polydata in subject space
        atlas.nystrom_polydata = pd
    
    #-----------------
    # Cluster the data using clusters from the atlas
    #-----------------
    output_polydata_s, cluster_numbers_s, color, embed = \
        wma.cluster.spectral_atlas_label(input_data, atlas, number_of_jobs=number_of_jobs)
    
    # Write the polydata with cluster indices saved as cell data
    # print '<wm_cluster_atlas.py> Saving output whole brain cluster file in directory:', outdir
    # fname_output = os.path.join(outdir, 'clustered_whole_brain.vtp')
    # wma.io.write_polydata(output_polydata_s, fname_output)
    
    # Figure out file name and mean color for each cluster, and write the individual polydatas
    fnames = list()
    cluster_colors = list()
    # These lines counted clusters present in the subject only
    ##number_of_clusters = np.max(cluster_numbers_s)
    ##first_cluster = np.min(cluster_numbers_s)
    # Get the number of clusters directly from the atlas to ensure we output files for all clusters.
    [number_of_clusters, number_of_eigenvectors] = atlas.centroids.shape
    print(f"<{os.path.basename(__file__)}> Cluster indices range from:", 0, "to", number_of_clusters)
    pd_c_list = wma.cluster.mask_all_clusters(output_polydata_s, cluster_numbers_s, number_of_clusters, preserve_point_data=True, preserve_cell_data=True,verbose=True)
    
    print(f'<{os.path.basename(__file__)}> Saving output cluster files in directory:', outdir)
    cluster_sizes = list()
    cluster_fnames = list()
    for c in range(number_of_clusters):
        mask = cluster_numbers_s == c
        cluster_size = np.sum(mask)
        cluster_sizes.append(cluster_size)
        #pd_c = wma.filter.mask(output_polydata_s, mask, preserve_point_data=True, preserve_cell_data=True,verbose=False)
        pd_c = pd_c_list[c]
        # The clusters are stored starting with 1, not 0, for user friendliness.
        fname_c = f'cluster_{c+1:05d}.vtp'
        # save the filename for writing into the MRML file
        fnames.append(fname_c)
        # prepend the output directory
        fname_c = os.path.join(outdir, fname_c)
        cluster_fnames.append(fname_c)
        wma.io.write_polydata(pd_c, fname_c)
        color_c = color[mask,:]
        if cluster_size:
            cluster_colors.append(np.mean(color_c,0))
        else:
            # avoid error if empty mean above
            cluster_colors.append(color[0,:])
            
    # Notify user if some clusters empty
    print(f"<{os.path.basename(__file__)}> Checking for empty clusters (can be due to anatomical variability or too few fibers analyzed).")
    for sz, fname in zip(cluster_sizes,cluster_fnames):
        if sz == 0:
            print(sz, ":", fname)
        
    cluster_sizes = np.array(cluster_sizes)
    print(f"<{os.path.basename(__file__)}> Mean number of fibers per cluster:", np.mean(cluster_sizes), "Range:", np.min(cluster_sizes), "..", np.max(cluster_sizes))
    
    # Estimate subsampling ratio to display approx. show_fibers total fibers
    number_fibers = len(cluster_numbers_s)
    if number_fibers < show_fibers:
        ratio = 1.0
    else:
        ratio = show_fibers / number_fibers
    print(f"<{os.path.basename(__file__)}> Subsampling ratio for display of", show_fibers, "total fibers estimated as:", ratio)
    
    # Write the MRML file into the directory where the polydatas were already stored
    fname = os.path.join(outdir, 'clustered_tracts.mrml')
    wma.mrml.write(fnames, np.around(np.array(cluster_colors), decimals=3), fname, ratio=ratio)
    
    # Also write one with 100%% of fibers displayed
    fname = os.path.join(outdir, 'clustered_tracts_display_100_percent.mrml')
    wma.mrml.write(fnames, np.around(np.array(cluster_colors), decimals=3), fname, ratio=1.0)
    
    # View the whole thing in png format for quality control
    if render:
    
        try:
            print(f'<{os.path.basename(__file__)}> Rendering and saving images of clustered subject.')
            ren = wma.render.render(output_polydata_s, 1000, data_mode='Cell', data_name='EmbeddingColor',verbose=False)
            ren.save_views(outdir)
            del ren
        except:
            print(f'<{os.path.basename(__file__)}> No X server available.')
    
    print("\n==========================")
    print(f'<{os.path.basename(__file__)}> Done clustering subject.  See output in directory:\n ', outdir, '\n')
    
    
if __name__ == '__main__':
    main()

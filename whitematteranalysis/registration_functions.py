try:
    import matplotlib.pyplot as plt
    USE_MATPLOTLIB = 1
except ImportError:
    USE_MATPLOTLIB = 0
    print "<registration_functions.py> Failed to import  matplotlib.pyplot, cannot plot."
    print "<registration_functions.py> Please install matplotlib for this functionality."
import numpy
import os
import time

import vtk

import whitematteranalysis as wma


def compute_multiscale_registration(register, scale_mode, n_steps, fiber_sample_size, sigma, maxfun):

    print "<register> SCALE:", scale_mode, "SIGMA:", sigma, "SAMPLES:", fiber_sample_size, "MAXFUN:", maxfun

    register.fiber_sample_size = fiber_sample_size
    register.sigma = sigma
    register.maxfun = maxfun
    
    if scale_mode == "Coarse":
        # relatively large steps
        inc_rot = (5.0 / 180.0) * numpy.pi
        inc_trans = 5.0
        inc_scale = 0.01
        inc_shear = (2.0 / 180.0) * numpy.pi
        register.set_rhobeg(inc_rot, inc_trans, inc_scale, inc_shear)    
        # relatively easy threshold to converge
        inc_rot = (4.5 / 180.0) * numpy.pi
        inc_trans = 4.5
        inc_scale = 0.01
        inc_shear = (2.0 / 180.0) * numpy.pi        
        register.set_rhoend(inc_rot, inc_trans, inc_scale, inc_shear)    
        # n = 5
        # only translation and rotation. initialization.
        for idx in range(0, n_steps):
            print "<register> SCALE:", scale_mode, idx+1, "/", n_steps
            register.translate_only()
            register.compute()
            register.rotate_only()
            register.compute()
    elif scale_mode == "Medium":
        # medium steps
        inc_rot = (4.0 / 180.0) * numpy.pi
        inc_trans = 4.0
        inc_scale = 0.01
        inc_shear = (2.0 / 180.0) * numpy.pi
        register.set_rhobeg(inc_rot, inc_trans, inc_scale, inc_shear)    
        # relatively easy threshold to converge
        inc_rot = (3.0 / 180.0) * numpy.pi
        inc_trans = 3.0
        inc_scale = 0.008
        inc_shear = (1.5 / 180.0) * numpy.pi        
        register.set_rhoend(inc_rot, inc_trans, inc_scale, inc_shear)    
        # n = 1
        for idx in range(0, n_steps):
            print "<register> SCALE:", scale_mode, idx+1, "/", n_steps
            register.translate_only()
            register.compute()
            register.rotate_only()
            register.compute()
            register.scale_only()
            register.compute()
            register.shear_only()
            register.compute()
    elif scale_mode == "Fine":
        # finer steps
        inc_rot = (3.0 / 180.0) * numpy.pi
        inc_trans = 3.0
        inc_scale = 0.008
        inc_shear = (1.5 / 180.0) * numpy.pi        
        register.set_rhobeg(inc_rot, inc_trans, inc_scale, inc_shear)    
        # smaller threshold to converge
        inc_rot = (2.0 / 180.0) * numpy.pi
        inc_trans = 2.0
        inc_scale = 0.006
        inc_shear = (1.0 / 180.0) * numpy.pi        
        register.set_rhoend(inc_rot, inc_trans, inc_scale, inc_shear)    
        # n = 1
        for idx in range(0, n_steps):
            print "<register> SCALE:", scale_mode, idx+1, "/", n_steps
            register.translate_only()
            register.compute()
            register.rotate_only()
            register.compute()
            register.scale_only()
            register.compute()
            register.shear_only()
            register.compute()
    elif scale_mode == "Finest":
        inc_rot = (1.0 / 180.0) * numpy.pi
        inc_trans = 1.0
        inc_scale = 0.005
        inc_shear = (1.0 / 180.0) * numpy.pi
        register.set_rhobeg(inc_rot, inc_trans, inc_scale, inc_shear)
        inc_rot = (0.5 / 180.0) * numpy.pi
        inc_trans = 0.5
        inc_scale = 0.001
        inc_shear = (0.75 / 180.0) * numpy.pi
        register.set_rhoend(inc_rot, inc_trans, inc_scale, inc_shear)    
        # n = 1
        for idx in range(0, n_steps):
            print "<register> SCALE:", scale_mode, idx+1, "/", n_steps
            register.translate_only()
            register.compute()
            register.rotate_only()
            register.compute()
            register.scale_only()
            register.compute()
            register.shear_only()
            register.compute()
            
    
def run_multisubject_registration(input_directory, outdir,
                                  number_of_fibers=150,
                                  fiber_sample_fractions=[.10, .20, .30, .40],
                                  parallel_jobs=2,
                                  points_per_fiber=5,
                                  sigma_per_scale=[30, 10, 10, 5],
                                  maxfun_per_scale=None,
                                  distance_method='Hausdorff', 
                                  verbose=True, 
                                  fiber_length=75,
                                  fibers_rendered=100,
                                  steps_per_scale=[10, 3, 2, 2]):

    elapsed = list()

    input_pds, subject_ids = wma.io.read_and_preprocess_polydata_directory(input_directory, fiber_length, number_of_fibers)

    number_of_datasets = len(input_pds)
    
    # figure out maximum function evals for optimizer if not requested
    if maxfun_per_scale is None:
        # figure out how many cobyla iterations are needed
        minfun = number_of_datasets 
        maxfun_per_scale = [minfun*10, minfun*10, minfun*15, minfun*30]        

    # figure out numbers of fibers to sample
    fiber_sample_sizes = (number_of_fibers * numpy.array(fiber_sample_fractions)).astype(int)

    # create registration object and apply settings
    register = wma.congeal.CongealTractography()
    register.parallel_jobs = parallel_jobs
    register.threshold = 0
    register.points_per_fiber = points_per_fiber
    register.distance_method = distance_method
    
    # add inputs to the registration
    for pd in input_pds:
        register.add_subject(pd)

    # view output data from the initialization
    outdir_current =  os.path.join(outdir, 'iteration_0')
    if not os.path.exists(outdir_current):
        os.makedirs(outdir_current)
    output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
    # save the current atlas representation to disk
    wma.registration_functions.save_atlas(output_pds, outdir_current)
    # save pictures of the current 'atlas' or registered data
    ren = wma.registration_functions.view_polydatas(output_pds, fibers_rendered)
    ren.save_views(outdir_current)
    del ren
    wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir_current)
    
    scales = ["Coarse", "Medium", "Fine", "Finest"]
    scale_idx = 0
    for scale in scales:
        start = time.time()
        # run the basic iteration of translate, rotate, scale
        compute_multiscale_registration(register, scale, steps_per_scale[scale_idx], fiber_sample_sizes[scale_idx], sigma_per_scale[scale_idx], maxfun_per_scale[scale_idx])
        elapsed.append(time.time() - start)
        scale_idx += 1
        
        # view output data from this big iteration
        if verbose | (scale == "Finest"):
            outdir_current =  os.path.join(outdir, 'iteration_'+str(scale_idx))
            if not os.path.exists(outdir_current):
                os.makedirs(outdir_current)
            output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
            # save the current atlas representation to disk
            wma.registration_functions.save_atlas(output_pds, outdir_current)
            # save pictures of the current 'atlas' or registered data
            ren = wma.registration_functions.view_polydatas(output_pds, fibers_rendered)
            ren.save_views(outdir_current)
            del ren
            if scale == "Finest":
                wma.registration_functions.transform_polydatas_from_disk(input_directory, register.convert_transforms_to_vtk(), outdir_current)
            wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir_current)
    
            plt.figure() # to avoid all results on same plot
            plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
            plt.savefig(os.path.join(outdir_current, 'objective_function.pdf'))
        
    return register, elapsed


def run_atlas_registration(input_dir,
                           input_atlas,
                           output_dir,
                           number_of_fibers=150,
                           fiber_sample_fractions=[.10, .20, .30, .40],
                           parallel_jobs=2,
                           points_per_fiber=5,
                           sigma_per_scale=[30, 10, 10, 5],
                           maxfun_per_scale=None,
                           distance_method='Hausdorff', 
                           verbose=True, 
                           fiber_length=75,
                           fibers_rendered=100,
                           steps_per_scale=[10, 3, 2, 2]):

    elapsed = list()

    input_pds, subject_ids = wma.io.read_and_preprocess_polydata_directory(input_dir, fiber_length, number_of_fibers)

    atlas = wma.io.read_polydata(input_atlas)

    transforms = list()
    
    number_of_datasets = len(input_pds)

    # figure out maximum function evals for optimizer if not requested
    if maxfun_per_scale is None:
        # default for two dataset registration
        maxfun_per_scale = [20, 40, 60, 80]

    # figure out numbers of fibers to sample
    fiber_sample_sizes = (number_of_fibers * numpy.array(fiber_sample_fractions)).astype(int)

    # register each pd to the atlas
    sidx = 0
    for pd in input_pds:
        # create registration object and apply settings
        register = wma.congeal.CongealTractography()
        register.parallel_jobs = parallel_jobs
        register.threshold = 0
        register.points_per_fiber = points_per_fiber
        register.distance_method = distance_method
    
        register.add_subject(pd)
        register.add_subject(atlas)
        
        # view output data from the initialization
        outdir_current =  os.path.join(output_dir, subject_ids[sidx], 'iteration_0')
        if not os.path.exists(outdir_current):
            os.makedirs(outdir_current)
        output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
        # save the current atlas representation to disk
        wma.registration_functions.save_atlas(output_pds, outdir_current)
        # save pictures of the current 'atlas' or registered data
        ren = wma.registration_functions.view_polydatas(output_pds, fibers_rendered)
        ren.save_views(outdir_current)
        del ren
        wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir_current)
        
        scales = ["Coarse", "Medium", "Fine", "Finest"]
        scale_idx = 0
        for scale in scales:
            start = time.time()
            # run the basic iteration of translate, rotate, scale
            compute_multiscale_registration(register, scale, steps_per_scale[scale_idx], fiber_sample_sizes[scale_idx], sigma_per_scale[scale_idx], maxfun_per_scale[scale_idx])
            elapsed.append(time.time() - start)
            scale_idx += 1
        
            # view output data from this big iteration
            if verbose | (scale == "Finest"):
                outdir_current =  os.path.join(output_dir, 'iteration_'+str(scale_idx))
                if not os.path.exists(outdir_current):
                    os.makedirs(outdir_current)
                output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
                # save the current atlas representation to disk
                wma.registration_functions.save_atlas(output_pds, outdir_current)
                # save pictures of the current 'atlas' or registered data
                ren = wma.registration_functions.view_polydatas(output_pds, fibers_rendered)
                ren.save_views(outdir_current)
                del ren
                if scale == "Finest":
                    wma.registration_functions.transform_polydatas_from_disk(input_directory, register.convert_transforms_to_vtk(), outdir_current)
                wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir_current)
    
                plt.figure() # to avoid all results on same plot
                plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
                plt.savefig(os.path.join(outdir_current, 'objective_function.pdf'))
        
        #  NEED TO OUTPUT THE COMBINED TRANSFORM APPLIED TO PD NOT TO ATLAS
        # now apply the appropriate transform to the input data.
        # transform 0 times transform 1 inverse
        tx = register.convert_transforms_to_vtk()
        tx[1].Inverse()
        tx[0].Concatenate(tx[1])

        transforms.append(tx[0])

        sidx = sidx + 1
        print elapsed
        del register

    transform_polydatas_from_disk(input_dir, transforms, output_dir)
    

def write_transforms_to_itk_format(transform_list, outdir):
    # use with the slicer module, something like this applies things
    # and handles the ras lps issues
    #./ResampleScalarVectorDWIVolume --spaceChange -b -c  -f /Users/odonnell/LinearTransform-2.tfm /Users/odonnell/Dropbox/Data/TBI_FE_PNL/controls_images/01231-dwi-filt-Ed-B0.nhdr test.nhdr
    idx = 0
    tx_fnames = list()
    for tx in transform_list:
        three_by_three = list()
        translation = list()
        for i in range(0,3):
            for j in range(0,3):
                three_by_three.append(tx.GetMatrix().GetElement(i,j))
        translation.append(tx.GetMatrix().GetElement(0,3))
        translation.append(tx.GetMatrix().GetElement(1,3))
        translation.append(tx.GetMatrix().GetElement(2,3))
        
        fname = 'txform_{0:05d}.tfm'.format(idx)
        fname = os.path.join(outdir, fname)
        tx_fnames.append(fname)
        f = open(fname, 'w')
        f.write('#Insight Transform File V1.0\n')
        f.write('# Transform 0\n')
        f.write('Transform: AffineTransform_double_3_3\n')
        f.write('Parameters: ')
        for el in three_by_three:
            f.write('{0} '.format(el))
        for el in translation:
            f.write('{0} '.format(el))
        f.write('\nFixedParameters: 0 0 0\n')

        idx +=1
    return(tx_fnames)

    
def save_atlas(polydata_list, out_dir):
    print "<registration_functions.py>: Save current atlas (sampled registered polydata)."
    appender = vtk.vtkAppendPolyData()
    idx = 0
    n_subj = len(polydata_list)
    for pd in polydata_list:
        nf = pd.GetNumberOfLines()
        print "<registration_functions.py> subject:", idx+1, "/" , n_subj, "fibers:", nf    
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            appender.AddInputData(pd)
        else:
            appender.AddInput(pd)
        idx = idx + 1
    appender.Update()
    wma.io.write_polydata(appender.GetOutput(), os.path.join(out_dir, 'atlas.vtk'))
    del appender

def view_polydatas(polydata_list, number_of_fibers=None):
    print "<registration_functions.py>: Appending downsampled polydatas for rendering with color by subject."
    appender = vtk.vtkAppendPolyData()
    idx = 0
    n_subj = len(polydata_list)
    for pd in polydata_list:
        nf0 = pd.GetNumberOfLines()
        if number_of_fibers is not None:
            # downsample if requested
            pd = wma.filter.downsample(pd, number_of_fibers)
        nf = pd.GetNumberOfLines()
        print "<registration_functions.py> subject:", idx+1, "/" , n_subj, "fibers:", nf,  "/" , nf0    
        mask = numpy.ones(nf)
        colors = numpy.multiply(mask, idx-1)
        pd2 = wma.filter.mask(pd, mask, colors)
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            appender.AddInputData(pd2)
        else:
            appender.AddInput(pd2)
        idx = idx + 1
    appender.Update()
    pd3 = appender.GetOutput()
    ren = wma.render.render(pd3)
    return ren

def transform_polydatas(input_pds, register):
    transforms = register.convert_transforms_to_vtk()
    idx = 0
    output_pds = list()
    for transform in transforms:
        transformer = vtk.vtkTransformPolyDataFilter()
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            transformer.SetInputData(input_pds[idx])
        else:
            transformer.SetInput(input_pds[idx])
        transformer.SetTransform(transform)
        transformer.Update()
        pd = transformer.GetOutput()
        output_pds.append(pd)
        idx = idx + 1
    return output_pds


def transform_polydatas_from_disk(input_dir, transforms, output_dir):

    # Find input files
    input_pd_fnames = wma.io.list_vtk_files(input_dir)
    num_pd = len(input_pd_fnames) 
    print "<registration_functions.py> ======================================="
    print "<registration_functions.py> Transforming vtk and vtp files from directory: ", input_dir
    print "<registration_functions.py> Total number of files found: ", num_pd
    print "<registration_functions.py> Writing output to directory: ", output_dir
    print "<registration_functions.py> ======================================="

    if not os.path.exists(output_dir):
        print "<registration_functions.py> ERROR: Output directory does not exist."
        return
    if not os.path.exists(input_dir):
        print "<registration_functions.py> ERROR: Output directory does not exist."        
        return
    
    #transforms = register.convert_transforms_to_vtk()
    
    for idx in range(0, len(input_pd_fnames)):

        fname = input_pd_fnames[idx]
        subject_id = os.path.splitext(os.path.basename(fname))[0]
        out_fname = os.path.join(output_dir, subject_id + '_reg.vtk')
        print "<registration_functions.py>  ", idx + 1, "/",  num_pd, subject_id, " Transforming ", fname, "->", out_fname, "..."

        pd = wma.io.read_polydata(fname)

        transformer = vtk.vtkTransformPolyDataFilter()
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            transformer.SetInputData(pd)
        else:
            transformer.SetInput(pd)
        transformer.SetTransform(transforms[idx])
        transformer.Update()
        
        pd2 = transformer.GetOutput()
        wma.io.write_polydata(pd2, out_fname)
        
        del transformer
        del pd2
        del pd
   
def smooth_polydatas_from_disk(input_pd_fnames, sigma, number_of_fibers, fiber_length, parallel_jobs):
    for fname in input_pd_fnames:
        print 'reading' , fname
        pd = wma.io.read_polydata(fname)
        print 'processing'
        pd = wma.filter.preprocess(pd,fiber_length)
        print 'random downsampling'
        pd = wma.filter.downsample(pd, number_of_fibers)

        pd2, weights = \
            wma.filter.smooth(pd, \
            fiber_distance_sigma = sigma,\
            points_per_fiber=35, \
            n_jobs=parallel_jobs, \
            upper_thresh=sigma*2)

        # write new one in current directory
        wma.io.write_polydata(pd2, 'gaussian_smooth_sigma{0}mm_'.format(sigma) + os.path.basename(fname))
    
        # write one without junk fibers also
        pdA = wma.filter.mask(pd2, weights >= 2, weights)
        wma.io.write_polydata(pdA, 'gaussian_smooth_sigma{0}mm_weight_gr_2_'.format(sigma) + os.path.basename(fname))

        # write one without junk fibers also
        pdA = wma.filter.mask(pd2, weights >= 5, weights)
        wma.io.write_polydata(pdA, 'gaussian_smooth_sigma{0}mm_weight_gr_5_'.format(sigma) + os.path.basename(fname))

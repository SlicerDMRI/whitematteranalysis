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
                                  steps_per_scale=[5, 3, 2, 1],
                                  no_render=False,
                                  midsag_symmetric=False,
                                  random_seed=None):

    # the big gain in objective from the 1st scale is early,
    # and actually after the fine (3rd scale) registration
    # the result looks good enough. So, save output after third step,
    # and also reduce the time spent in first and last scales.
    # note: this is still more computation than in the publication.
    # was:
    # steps_per_scale=[10, 3, 2, 2]
    # now using:
    # steps_per_scale=[5, 3, 2, 1]

    elapsed = list()

    input_pds, subject_ids = wma.io.read_and_preprocess_polydata_directory(input_directory, fiber_length, number_of_fibers, random_seed)

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
    register.random_seed = random_seed
    
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
    if no_render:
        print "<register> Intermediate rendering OFF"
    else:
        ren = wma.registration_functions.view_polydatas(output_pds, fibers_rendered)
        ren.save_views(outdir_current)
        del ren
    wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir_current, subject_ids)
    
    scales = ["Coarse", "Medium", "Fine", "Finest"]
    scale_idx = 0
    iter_count = 0
    for scale in scales:
        start = time.time()
        # run the basic iteration of translate, rotate, scale
        compute_multiscale_registration(register, scale, steps_per_scale[scale_idx], fiber_sample_sizes[scale_idx], sigma_per_scale[scale_idx], maxfun_per_scale[scale_idx])
        elapsed.append(time.time() - start)
        scale_idx += 1
        
        # view output data from this big iteration
        outdir_current =  os.path.join(outdir, 'iteration_'+str(scale_idx))
        if not os.path.exists(outdir_current):
            os.makedirs(outdir_current)
        # Only save and render transformed data if near the end of the registration or verbose is on
        if verbose | (scale == "Fine") | (scale == "Finest"):
            output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
            # save the current atlas representation to disk
            wma.registration_functions.save_atlas(output_pds, outdir_current)
            # save pictures of the current 'atlas' or registered data
            if no_render:
                print "<register> Intermediate rendering OFF"
            else:
                ren = wma.registration_functions.view_polydatas(output_pds, fibers_rendered)
                ren.save_views(outdir_current)
                del ren
            # if half of the brains were reflected find the corrrect transforms for the data on disk
            transform_list = register.convert_transforms_to_vtk()
            if midsag_symmetric:
                transform_list = transform_list[::2]
            # transform input polydata from disk to keep all attributes (tensors etc) in place
            # save these inside a subdirectory to avoid accidental clustering with the "registration_atlas.vtk"
            outdir_pd = os.path.join(outdir_current, 'registered_subjects')
            if not os.path.exists(outdir_pd):
                os.makedirs(outdir_pd)
            wma.registration_functions.transform_polydatas_from_disk(input_directory, transform_list, outdir_pd)
        # Always save the current transforms and the objective function plot to check progress
        wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir_current, subject_ids)
        # Number of objective function computations this big iteration
        iter_count_new = len(register.objective_function_values)
        if no_render:
            print "<register> Intermediate saving of objective function plot off (no render)."
        else:
            plt.figure() # to avoid all results on same plot
            #plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
            # plot just the most recent scale's objective. otherwise the y axis range is too large.
            plt.plot(range(iter_count_new - iter_count), register.objective_function_values[iter_count:iter_count_new])
            plt.savefig(os.path.join(outdir_current, 'objective_function.pdf'))
        iter_count = iter_count_new


    # If we are registering for symmetry, make sure group coordinate system is midsagitally aligned
    # Empirically the group registration is well aligned and the bilateral clustering 
    # was very successful even before symmetric registration, so this is not needed now.
    # But it could be implemented/tested later especially if registration is extended to a warp.
    #if midsag_symmetric:
    #    print "LAUREN ALIGN ATLAS TO SELF ACROSS MIDSAG LINE FOR COMPLETENESS AND APPLY TXFORM TO ALL SUBJECTS"

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
        reg_pds = list()
        reg_pds.append(pd)
        reg_pds.append(atlas)
        
        # view output data from the initialization
        outdir_current =  os.path.join(output_dir, subject_ids[sidx], 'iteration_0')
        if not os.path.exists(outdir_current):
            os.makedirs(outdir_current)
        output_pds = wma.registration_functions.transform_polydatas(reg_pds, register)
        # save pictures of the current data
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
                outdir_current =  os.path.join(output_dir, subject_ids[sidx], 'iteration_'+str(scale_idx))
                if not os.path.exists(outdir_current):
                    os.makedirs(outdir_current)
                output_pds = wma.registration_functions.transform_polydatas(reg_pds, register)
                # save pictures of the current registered data
                ren = wma.registration_functions.view_polydatas(output_pds, fibers_rendered)
                ren.save_views(outdir_current)
                del ren
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
    return elapsed

def write_transforms_to_itk_format(transform_list, outdir, subject_ids=None):
    # This now outputs an ITK transform that works correctly to transform the tracts (or any volume in the same space) in Slicer.
    # Also output is a vtk transform, using MNI format to have vtk reading/writing for wma code that does not depend on itk now.

    ## OLD comment. May not work now, as no space change is needed.
    ## use with the slicer module, something like this applies things
    ## and handles the ras lps issues
    ##./ResampleScalarVectorDWIVolume --spaceChange -b -c  -f /Users/odonnell/LinearTransform-2.tfm /Users/odonnell/Dropbox/Data/TBI_FE_PNL/controls_images/01231-dwi-filt-Ed-B0.nhdr test.nhdr

    idx = 0
    tx_fnames = list()
    for tx in transform_list:

        # save out the vtk transform to a text file as it is
        # The MNI transform reader/writer are available in vtk so use those:
        writer = vtk.vtkMNITransformWriter()
        writer.AddTransform(tx)
        if subject_ids is not None:
            fname = 'txform_vtk_' + str(subject_ids[idx]) + '.xfm'
        else:
            fname = 'txform_vtk_{0:05d}.xfm'.format(idx)
        writer.SetFileName(os.path.join(outdir, fname))
        writer.Write()
        
        # Save the itk transform as the inverse of this transform (resampling transform) and in LPS.
        # This will show the same transform in the slicer GUI as the vtk transform we internally computed
        # that is stored in the .xfm text file, above.
        # To apply our transform to resample a volume in LPS:
        # convert to RAS, use inverse of transform to resample, convert back to LPS
        tx_inverse = vtk.vtkTransform()
        tx_inverse.DeepCopy(tx)
        tx_inverse.Inverse()
        ras_2_lps = vtk.vtkTransform()
        ras_2_lps.Scale(-1, -1, 1)
        lps_2_ras = vtk.vtkTransform()
        lps_2_ras.Scale(-1, -1, 1)
        tx2 = vtk.vtkTransform()
        tx2.Concatenate(lps_2_ras)
        tx2.Concatenate(tx_inverse)
        tx2.Concatenate(ras_2_lps)

        three_by_three = list()
        translation = list()
        for i in range(0,3):
            for j in range(0,3):
                three_by_three.append(tx2.GetMatrix().GetElement(i,j))
        translation.append(tx2.GetMatrix().GetElement(0,3))
        translation.append(tx2.GetMatrix().GetElement(1,3))
        translation.append(tx2.GetMatrix().GetElement(2,3))
        
        if subject_ids is not None:
            fname = 'txform_itk_' + str(subject_ids[idx]) + '.tfm'
        else:
            fname = 'txform__itk_{0:05d}.tfm'.format(idx)
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
    wma.io.write_polydata(appender.GetOutput(), os.path.join(out_dir, 'registration_atlas.vtk'))
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
            pd = wma.filter.downsample(pd, number_of_fibers, verbose=False)
        nf = pd.GetNumberOfLines()
        print "<registration_functions.py> subject:", idx+1, "/" , n_subj, "fibers:", nf,  "/" , nf0    
        mask = numpy.ones(nf)
        colors = numpy.multiply(mask, idx-1)
        pd2 = wma.filter.mask(pd, mask, colors, verbose=False)
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            appender.AddInputData(pd2)
        else:
            appender.AddInput(pd2)
        idx = idx + 1
    appender.Update()
    pd3 = appender.GetOutput()
    ren = wma.render.render(pd3, verbose=False)
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


def run_midsag_align(input_poly_data, outdir, number_of_fibers=150,
    number_of_fibers_per_step=[75, 75, 75, 100],
    parallel_jobs=2,
    points_per_fiber=5,
    sigma_per_step=[30, 10, 10, 5],
    maxfun=150,
    distance_method='Hausdorff',
    fiber_length=60,
    verbose=True):

    # note: doing this twice and averaging the output transforms seems
    # to work well, so consider that when this code is improved.

    elapsed = list()
    number_of_fibers_step_one = number_of_fibers_per_step[0]
    number_of_fibers_step_two = number_of_fibers_per_step[1]
    number_of_fibers_step_three = number_of_fibers_per_step[2]
    number_of_fibers_step_four = number_of_fibers_per_step[3]

    minfun = maxfun/5.0
    #maxfun_per_step = [50, 75, 200]
    maxfun_per_step = [minfun*1.5, minfun*2, minfun*3, minfun*5]
    maxfun_per_step = [minfun*2, minfun*3, minfun*4, minfun*5]

    sigma_step_one = sigma_per_step[0]
    sigma_step_two = sigma_per_step[1]
    sigma_step_three = sigma_per_step[2]
    sigma_step_four = sigma_per_step[3]

    print 'Read and preprocess'
    input_pds = list()

    pd = wma.io.read_polydata(input_poly_data)
    pd2 = wma.filter.preprocess(pd, fiber_length)
    pd3 = wma.filter.downsample(pd2, number_of_fibers)
    pd4 = wma.filter.downsample(pd2, number_of_fibers)
    # release memory
    del pd
    del pd2
    trans = vtk.vtkTransform()
    trans.Scale(-1,1,1)
    transformer = vtk.vtkTransformPolyDataFilter()
    transformer.SetTransform(trans)
    if (vtk.vtkVersion().GetVTKMajorVersion() > 5.0):
        transformer.SetInputData(pd4)
    else:
        transformer.SetInput(pd4)
    transformer.Update()
    pd5 = transformer.GetOutput()
    # append input and its reflection to list of data to register
    input_pds.append(pd3)
    input_pds.append(pd5)

    # create registration object and apply settings
    register = wma.congeal.CongealTractography()
    register.parallel_jobs = parallel_jobs
    register.threshold = 0
    register.points_per_fiber = points_per_fiber
    register.distance_method = distance_method
    #register.maxfun = maxfun
  
    # add inputs to the registration
    for pd in input_pds:
        register.add_subject(pd)

    # view input data from the initialization
    if verbose:
        outdir_current =  os.path.join(outdir, 'iteration_0')
        if not os.path.exists(outdir_current):
            os.makedirs(outdir_current)
        output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
        ren = wma.registration_functions.view_polydatas(output_pds, 1000)
        ren.save_views(outdir_current)
        #wma.registration_functions.transform_polydatas_from_disk(input_poly_data, register, outdir_current)
        wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir)
    
    # STEP ONE
    # make sure we take very small steps, the brains are already overlapping
    # step one has larger rotate and translate to avoid local minima if brains quite rotated
    inc_rot = (5.0 / 180.0) * numpy.pi
    inc_trans = 2.0
    inc_scale = 0.001
    inc_shear = (.5 / 180.0) * numpy.pi
    register.set_rhobeg(inc_rot, inc_trans, inc_scale, inc_shear)
  
    start = time.time()
    # run the basic iteration of translate, rotate, scale
    register.fiber_sample_size = number_of_fibers_step_one
    register.sigma = sigma_step_one
    register.maxfun = maxfun_per_step[0]
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    # Don't scale for midsagittal alignment
    #register.scale_only()
    #register.compute()
    elapsed.append(time.time() - start)

    # STEP TWO
    # make sure we take very small steps, the brains are already overlapping
    # step one has larger rotate and translate to avoid local minima if brains quite rotated
    inc_rot = (0.5 / 180.0) * numpy.pi
    inc_trans = 0.5
    inc_scale = 0.001
    inc_shear = (.5 / 180.0) * numpy.pi
    register.set_rhobeg(inc_rot, inc_trans, inc_scale, inc_shear)
    start = time.time()
    # run the basic iteration of translate, rotate, scale AGAIN
    register.fiber_sample_size = number_of_fibers_step_two
    register.sigma = sigma_step_two
    register.maxfun = maxfun_per_step[1]    
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    #register.scale_only()
    register.compute()
    register.shear_only()
    register.compute()
    elapsed.append(time.time() - start)

    #if 0:
    # STEP THREE
    start = time.time()
    # run the basic iteration of translate, rotate, scale AGAIN
    register.fiber_sample_size = number_of_fibers_step_three
    register.sigma = sigma_step_three
    register.maxfun = maxfun_per_step[2]
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    #register.scale_only()
    register.compute()
    register.shear_only()
    register.compute()
    elapsed.append(time.time() - start)

    # STEP FOUR
    start = time.time()
    # run the basic iteration of translate, rotate, scale AGAIN
    register.fiber_sample_size = number_of_fibers_step_four
    register.sigma = sigma_step_four
    register.maxfun = maxfun_per_step[3]
    register.translate_only()
    register.compute()
    register.rotate_only()
    register.compute()
    #register.scale_only()
    register.compute()
    register.shear_only()
    register.compute()
    elapsed.append(time.time() - start)

    # STEP FIVE
    # run more with the final settings.
    # See if shear is improved
    start = time.time()
    for idx in range(3):
        # run the basic iteration of translate, rotate, scale AGAIN
        register.translate_only()
        register.compute()
        register.rotate_only()
        register.compute()
        #register.scale_only()
        register.compute()
        register.shear_only()
        register.compute()
    elapsed.append(time.time() - start)
    
    # view output data from this big iteration
    subjectID = os.path.splitext(os.path.basename(input_poly_data))[0]
    if verbose:
        outdir_current =  os.path.join(outdir, 'iteration_4')
        if not os.path.exists(outdir_current):
            os.makedirs(outdir_current)
        output_pds = wma.registration_functions.transform_polydatas(input_pds, register)
        ren = wma.registration_functions.view_polydatas(output_pds, 500)
        ren.save_views(outdir_current)
        #wma.registration_functions.transform_polydatas_from_disk(input_poly_data, register, outdir_current)
        wma.registration_functions.write_transforms_to_itk_format(register.convert_transforms_to_vtk(), outdir)
    
    plt.figure() # to avoid all results on same plot
    plt.plot(range(len(register.objective_function_values)), register.objective_function_values)
    plt.savefig(os.path.join(outdir, 'objective_function_' + str(subjectID) + '.pdf'))
    plt.close()

    #return register, elapsed, input_pds

    print "TIME:", elapsed

    print "Re-reading and transforming original data to aligned. Writing outputs."
    
    # now apply the appropriate transform to the input data.
    # half of transform 1 times transform 2 inverse
    
    tx = register.convert_transforms_to_vtk()
    print tx[0]
    print tx[1]
    tx[1].Inverse()
    tx[0].Concatenate(tx[1])

    print tx[0]

    txs = vtk.vtkTransform()
    m = vtk.vtkMatrix4x4()
    for i in range(0,4):
        for j in range(0,4):
            # scale is 1, and value 3,3 is 1
            if i == j:
                m.SetElement(i,j, 1.0)
            else:
                el = tx[0].GetMatrix().GetElement(i,j)
                m.SetElement(i,j, el / 2.0)

    txs.SetMatrix(m)

    print m

    # save this transform to disk
    trans_list = []
    trans_list.append(txs)
    id_list = []
    id_list.append('midsag_' + str(subjectID))
    wma.registration_functions.write_transforms_to_itk_format(trans_list, outdir, id_list)

    # now apply it to the input data
    pd = wma.io.read_polydata(input_poly_data)

    trans = vtk.vtkTransformPolyDataFilter()
    trans.SetTransform(txs)
    if (vtk.vtkVersion().GetVTKMajorVersion() > 5.0):
        trans.SetInputData(pd)
    else:
        trans.SetInput(pd)
    trans.Update()


    #fname1 = os.path.split(input_poly_data)[1]
    #fname1 = os.path.splitext(fname1)[0]
    fname1 = os.path.join(outdir, subjectID+'_sym.vtp')
    print "Writing output polydata..."
    wma.io.write_polydata(trans.GetOutput(), fname1)

    # get subject identifier from unique input filename
    # this way the views will not all overwrite each other 
    outdir_current = os.path.join(outdir, subjectID)
    if not os.path.exists(outdir_current):
        os.makedirs(outdir_current)

    # make the picture. View input, reflection, and output
    print "Rendering result..."
    all_pds = list()
    all_pds.append(input_pds[0])
    all_pds.append(input_pds[1])
    all_pds.append(trans.GetOutput())
    ren = wma.registration_functions.view_polydatas(all_pds, 250)
    ren.save_views(outdir_current)

    # also render just the output
    outdir_current2 = os.path.join(outdir_current, 'view_output')
    if not os.path.exists(outdir_current2):
        os.makedirs(outdir_current2)
    ren = wma.render.render(trans.GetOutput(), 1000)
    ren.save_views(outdir_current2)
 

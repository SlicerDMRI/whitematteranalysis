try:
    import matplotlib.pyplot as plt
    USE_MATPLOTLIB = 1
except ImportError:
    USE_MATPLOTLIB = 0
    print "<registration_functions.py> Failed to import  matplotlib.pyplot, cannot plot."
    print "<registration_functions.py> Please install matplotlib for this functionality."
import numpy
import os

import vtk

import whitematteranalysis as wma

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


def transform_polydatas_from_disk(input_dir, register, output_dir):

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
    
    transforms = register.convert_transforms_to_vtk()
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

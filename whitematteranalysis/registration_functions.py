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
    appender = vtk.vtkAppendPolyData()
    idx = 0
    for pd in polydata_list:
        if number_of_fibers is not None:
            # downsample if requested
            pd = wma.filter.downsample(pd, number_of_fibers)
        nf = pd.GetNumberOfLines()
        print idx
        print nf
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


def transform_polydatas_from_disk(input_pd_fnames, register, outdir):
    transforms = register.convert_transforms_to_vtk()
    for idx in range(0, len(input_pd_fnames)):
        transformer = vtk.vtkTransformPolyDataFilter()
        fname = input_pd_fnames[idx]
        print fname
        pd = wma.io.read_polydata(fname)
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            transformer.SetInputData(pd)
        else:
            transformer.SetInput(pd)
        
        transformer.SetTransform(transforms[idx])
        transformer.Update()
        pd2 = transformer.GetOutput()
        writer = vtk.vtkPolyDataWriter()
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            writer.SetInputData(pd2)
        else:
            writer.SetInput(pd2)
        
        # this seemed to not work correctly for re-reading by wma or slicer
        #fname = 'white_matter_{:04}.vtp'.format(idx)
        fname = 'white_matter_{:04}.vtk'.format(idx)
        writer.SetFileName(os.path.join(outdir, fname))
        #writer.SetFileName(os.path.join(outdir, str(idx) + '.vtk'))
        writer.SetFileTypeToBinary()
        writer.Write()
        del writer
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

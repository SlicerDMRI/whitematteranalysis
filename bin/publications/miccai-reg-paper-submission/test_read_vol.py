import glob
import os

import vtk

#indir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/MICCAI2012/control_image_data'

#fname = '/Users/odonnell/Data/TBI_FE_PNL/Tracts/controls/intermediate_data/01395-dwi-filt-Ed-B0.raw.gz'
#fname = '/Users/odonnell/Data/TBI_FE_PNL/Tracts/controls/intermediate_data/01395-dwi-filt-Ed-B0.raw'
#input_mask = "{0}/*.raw".format(indir)
#input_image_files = glob.glob(input_mask)

def write_command_line_to_convert_images(fname, image_fnames, txform_fnames):
    idx = 0
    vol_fnames = list()
    f = open(fname, 'w')
    for image in image_fnames:
        f.write('/Applications/Slicer.app/Contents/cli-modules/ResampleScalarVectorDWIVolume --spaceChange -b -c  -f ')
        f.write(' ./')
        f.write(txform_fnames[idx])
        f.write(' ')
        raw_fname = image
        basename, extension = os.path.splitext(raw_fname)
        f.write(basename + '.nhdr')
        f.write(' ')
        out_fname = 'tx_image_{0:05d}.nhdr'.format(idx)
        vol_fnames.append(out_fname)
        f.write(out_fname)
        f.write('\n')
        idx += 1
    return vol_fnames

def write_transforms_to_itk_format(transform_list):
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
        
        fname = 'txform_{0}.tfm'.format(idx)
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
    
    
def read_all_images(image_files):
    image_list = list()
    for fname in image_files:
        print "Reading:", fname
        reader = vtk.vtkImageReader()
        # if they gave us a nhdr get the raw file,
        # this is very basic 
        basename, extension = os.path.splitext(fname)
        if extension == '.nhdr':
            fname = basename + '.raw'
        print fname
        reader.SetFileName(fname)
        # header info (should be same for ALL TBI/control data)
        #reader.SetDataScalarTypeToShort()
        reader.SetDataScalarTypeToFloat()
        reader.SetDataSpacing(1.6667, 1.6667, 1.7)
        reader.SetDataOrigin(119.169, 119.169, 71.4)
        reader.SetDataByteOrderToLittleEndian()
        #reader.SetDataExtent(0, 144, 0, 144, 0, 85)
        reader.SetDataExtent(0, 143, 0, 143, 0, 84)
        reader.SetFileDimensionality(3)
        reader.Update()
        image_list.append(reader.GetOutput())
    del reader
    return image_list

def average_all_images(image_datas, fname='average_image.raw'):
    # how cool is this, this class was from my masters' thesis.
    # though the new vtk coding is not recognizable anymore...
    sum = vtk.vtkImageWeightedSum()
    weights = vtk.vtkDoubleArray()
    for image in image_datas:
        weights.InsertNextTuple1(1.0)
        sum.AddInput(image)
    sum.SetWeights(weights)
    sum.NormalizeByWeightOn()
    sum.Update()
    writer = vtk.vtkImageWriter()
    writer.SetInput(sum.GetOutput())
    writer.SetFileName(fname)
    writer.SetFileDimensionality(3)
    writer.Write()
    return(sum.GetOutput())

def transform_all_images_lps_to_ras(image_datas):
    image_list = list()
    for image in image_datas:
        resample = vtk.vtkImageReslice()
        resample.SetResliceAxesDirectionCosines(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0)
        resample.SetInput(image)
        #resample.SetResliceTransform(trans)
        #print trans
        #resample.SetInformationInput(info.GetOutput())
        resample.Update()
        #resample.SetOutputOrigin(origin)
        image_list.append(resample.GetOutput())
    del resample
    return image_list

       
def transform_all_images(register, image_datas, defining_image):
    image_list = list()
    idx = 0
    transform_list = register.convert_transforms_to_vtk()
    image_datas = [image_datas[0]]
    info = vtk.vtkImageChangeInformation()
    info.SetInput(defining_image)
    info.CenterImageOn()
        
    for image in image_datas:
        change = vtk.vtkImageChangeInformation()
        change.SetInput(image)
        change.CenterImageOn()
        #change.SetOutputSpacing(1,1,1)
        #spacing = image.GetSpacing()
        
        if 1:
            trans = vtk.vtkTransform()
            # this is right except scaling is 1/correct
            tx = transform_list[idx]
            trans.Concatenate(tx)
            #trans = tx
            #trans.Scale(-1,-1,-1)
            trans.Inverse()
            #trans.Scale(-1,-1,-1)
        if 0:
            transform = numpy.zeros(15)
            transform2 = numpy.zeros(15)
            tform = register._subjects[idx].transform
            # rotation and translation are "already inverted"
            # from images being in LPS
            transform[0:6] = tform[0:6]
            # invert scaling
            transform[6:9] = numpy.divide(1.0, tform[6:9])
            # invert shear
            transform2[6:9] = 1.0
            transform2[9:15] = tform[9:15]             
            reg_ops = wma.register.RegistrationInformation()
            trans = reg_ops.convert_transform_to_vtk(transform)
            trans2 = reg_ops.convert_transform_to_vtk(transform2)
            trans2.Scale(-1,-1,-1)
            trans2.Inverse()
            trans.Concatenate(trans2)
            
        resample = vtk.vtkImageReslice()
        #resample.SetResliceAxesDirectionCosines(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
        #resample.SetResliceAxesDirectionCosines(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0)
        #resample.SetInput(image)
        resample.SetInput(change.GetOutput())
        resample.SetResliceTransform(trans)
        print trans
        #resample.SetInformationInput(info.GetOutput())
        #resample.SetOutputSpacing(spacing)
        resample.Update()
        #resample.SetOutputOrigin(origin)
        image_list.append(resample.GetOutput())
        idx += 1
    del resample
    return image_list

def transform_all_images2(register, image_datas, defining_image):
    image_list = list()
    idx = 0
    transform_list = register.convert_transforms_to_vtk()
    #image_datas = [image_datas[0]]
    for image in image_datas:
        # Calculate the center of the volume
        #reader.GetOutput().UpdateInformation()
        (xMin, xMax, yMin, yMax, zMin, zMax) =image.GetWholeExtent()
        (xSpacing, ySpacing, zSpacing) = image.GetSpacing()
        (x0, y0, z0) = image.GetOrigin()
        
        center = [x0 + xSpacing * 0.5 * (xMin + xMax),
                  y0 + ySpacing * 0.5 * (yMin + yMax),
                  z0 + zSpacing * 0.5 * (zMin + zMax)]
        print center
        print image.GetCenter()
        ijk_to_ras = vtk.vtkMatrix4x4()
        ijk_to_ras.SetElement(0,0,-1)
        ijk_to_ras.SetElement(1,1,-1)
        ijk_to_ras.SetElement(2,2,-1)
        ras_to_ijk = ijk_to_ras
        
        trans = vtk.vtkTransform()

        #
        # from ras space in output image-land
        # 
        #trans.Concatenate(ijk_to_ras)
        trans.Translate(center[0], center[1], center[2])


        #trans.Translate(100, 0, 0)
        # here we need to use the register object instead of this test
        #trans = vtk.vtkTransform()
        #trans.GetMatrix().SetElement(0,0,0.8)
        # ijk to ras, and voxel size removal? 1.6667, 1.6667, 1.7

        origin = image.GetOrigin()
        print origin
        #image.SetOrigin(-origin[0], -origin[1], -origin[2])
        #origin = image.GetOrigin()
        #print origin
        #trans.Translate(origin[0]/2, origin[1]/2, origin[2]/2)
        #trans.Translate(origin[0], origin[1], origin[2])
        #trans.Scale(1.7, 1.7, 1.7)

        #trans.Scale(-1, -1, -1)
        #trans = transform_list[idx]
        tx = transform_list[idx]
        #tx.Inverse()
        #tx.Identity()

        #trans.Concatenate(tx)

        trans.Concatenate(ijk_to_ras)
        
        tform = register._subjects[idx].transform
        if 0:
            trans.Scale(1/tform[6], 1/tform[7], 1/tform[8])
            trans.Translate(-tform[3], -tform[4], -tform[5])
            trans.RotateZ(-tform[2] * 180/numpy.pi)
            trans.RotateY(-tform[1] * 180/numpy.pi)            
            trans.RotateX(-tform[0] * 180/numpy.pi)

        if 1:
            trans.RotateX(tform[0] * 180/numpy.pi)
            trans.RotateY(tform[1] * 180/numpy.pi)
            trans.RotateZ(tform[2] * 180/numpy.pi)
            trans.Translate(tform[3], tform[4], tform[5])
            trans.Scale(tform[6], tform[7], tform[8])

        #trans.RotateZ(10)

        trans.Concatenate(ras_to_ijk)
        
        trans.Translate(-center[0], -center[1], -center[2])
        #trans.Concatenate(ras_to_ijk)
        
        #trans.Translate(-origin[0]/2, -origin[1]/2, -origin[2]/2)
        #trans.Inverse()
        
        #trans.Scale(tform[6]/1.7, tform[7]/1.7, tform[8]/1.7)
        #trans.Scale(1/tform[6], 1/tform[7], 1/tform[8])
        #, tform[1], tform[2])
        #trans.Translate(-origin[0], -origin[1], -origin[2])
        #trans.Translate(-origin[0], origin[1], origin[2])        
        # ras to ijk
        #trans.Scale(-1, -1, -1)
        #trans.Scale(-1.7, -1.7, -1.7)
        #trans.Scale(-1/1.7, -1/1.7, -1/1.7)
        #trans.Translate(origin[0], origin[1], origin[2])
        # the transform to resample should be the inverse of the one
        # applied to tract points
        #trans.Inverse()
        #resample = vtk.vtkImageResample()
        resample = vtk.vtkImageReslice()
        resample.SetInput(image)
        resample.SetResliceTransform(trans)
        print trans
        resample.SetInformationInput(defining_image)
        resample.Update()
        #resample.SetOutputOrigin(origin)
        image_list.append(resample.GetOutput())
        idx += 1
    del resample
    return image_list

#image_list = read_all_images(input_image_files)
#average_image = average_all_images(image_list)

#test_list = transform_all_images(None, image_list, image_list[0])
#average_image2 = average_all_images(test_list)

# to view the data
#execfile('/Users/odonnell/Dropbox/Coding/Python/WhiteMatterAnalysis/bin/test_image_reslice_view.py')
#ret = view_image(average_image)
#execfile('/Users/odonnell/Dropbox/Coding/Python/WhiteMatterAnalysis/bin/test_image_reslice_view.py')

#writer = vtk.vtkImageWriter()
#writer.SetInput(average_image)
#writer.SetFileName('average_image.raw')
#writer.SetFileDimensionality(3)
#writer.Write()

# Calculate the center of the volume
#image.UpdateInformation()
#(xMin, xMax, yMin, yMax, zMin, zMax) = image.GetWholeExtent()
#(xSpacing, ySpacing, zSpacing) = image.GetSpacing()
#(x0, y0, z0) = image.GetOrigin()

#center = [x0 + xSpacing * 0.5 * (xMin + xMax),
#          y0 + ySpacing * 0.5 * (yMin + yMax),
#          z0 + zSpacing * 0.5 * (zMin + zMax)]



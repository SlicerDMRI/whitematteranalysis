import glob

import scipy.stats
import numpy
import matplotlib.pyplot as plt
import vtk
import scipy.signal
import time

import whitematteranalysis as wma


def image_movie(image, slice1, slice2, skip, orientation='ax'):
    limits = [numpy.min(image), numpy.max(image)]
    print "limits:", limits
    #limits = [-1, 1]
    limits = [-.6, .6]
    #plt.figure()
    for slice in range(slice1,slice2,skip):
        plt.imshow(image[:,:,slice,:])
        #plt.clim(limits)            
        #plt.imshow(image[:,slice,:])
        #plt.imshow(image[slice,:,:])
        plt.draw()
        fname = 'slice_{:04}.png'.format(slice)
        #plt.savefig('slice%d.png' % slice)
        plt.savefig(fname)
        time.sleep(0.2)

#inputDirectory = '/Users/lauren/Desktop/OUTPUT'
#inputDirectory = '/Users/lauren/Dropbox/Coding/OUTPUTS/laterality_test_mp'
#inputDirectory = '/Users/lauren/Dropbox/Coding/OUTPUTS/laterality_test_register/output2'
#inputDirectory = '/Users/odonnell/Desktop/OUTPUT/laterality_register/iteration_2'
inputDirectory = '/Users/odonnell/Desktop/OUTPUT/laterality_register/iteration_3'

inputMask1 = "{0}/*.vtk".format(inputDirectory);
input_fnames = glob.glob(inputMask1)

#input_fnames = input_fnames[0:3]

#input_fnames = input_fnames[22:-1]
#readPolyDatas = 0

data_fiber_counts = list()
data_fiber_LI = list()
# Read data, and check if the data are all present
for fname in input_fnames:
    try:
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(fname)
        reader.Update()
        vol, LI = wma.filter.pd_to_array(reader.GetOutput())
        data_fiber_counts.append(vol)
        data_fiber_LI.append(LI)
        del reader
    except:
        print "unexpected error" , fname
    else:
        print "Read success", fname
        #print resultIO.laterality_index[0]

# test summing all
sum_counts = numpy.zeros(data_fiber_counts[0].shape)
sum_LI = numpy.zeros(data_fiber_counts[0].shape)
sum_LIR = numpy.zeros(data_fiber_counts[0].shape)
sum_LIL = numpy.zeros(data_fiber_counts[0].shape)
#for sidx in range(0,len(input_fnames)):
for sidx in range(0,len(data_fiber_counts)):
    sum_counts += (data_fiber_counts[sidx] > 0)
    # LI average across fibers, so from -1 to 1
    # only do this division where there is a fiber
    idx = numpy.nonzero(data_fiber_counts[sidx])
    LI = numpy.zeros(data_fiber_counts[0].shape)
    LI[idx] = numpy.divide(data_fiber_LI[sidx][idx], data_fiber_counts[sidx][idx])
    sum_LI += LI
    sum_LIR[LI > 0] += LI[LI > 0]
    sum_LIL[LI < 0] += LI[LI < 0]


flip_sum_LIL = sum_LIL[::-1,:,:]
symm_LIL = sum_LIL + flip_sum_LIL

flip_sum_LIR = sum_LIR[::-1,:,:]
symm_LIR = sum_LIR + flip_sum_LIR

right_img = wma.filter.array_to_vtk(symm_LIR)
writer = vtk.vtkImageWriter()
writer.SetFilePrefix('right')
writer.SetInput(right_img)
writer.SetFileDimensionality(2)
writer.SetFilePattern('%s.%03d')
writer.Write()

left_img = wma.filter.array_to_vtk(symm_LIL)
writer.SetFilePrefix('left')
writer.SetInput(left_img)
writer.Write()



# Read groups file, in theory
handedness = numpy.array([75, 85, 55, -35, -90, -85, -95, 70, -75, -25, 60, -70, 75, 100, 75, -35, -45, 90, 85, -75, 5, -80, 50, -40, -95, 100])
fMRIlat  = numpy.array([0.25, 0.388, 0.291, 0.1885, -0.1605, 0.4435, 0.6385, 0.0455, 0.3805, -0.358, 0.6595, 0.3045, 0.537, 0.236, 0.202, 0.321, 0.4405, 0.038, 0.435, 0.304, 0.1385, 0.667, 0.337, 0.699, 0.31, 0.3925])

gL = numpy.nonzero(handedness < 0)[0]
gR = numpy.nonzero(handedness >= 0)[0]

# test remove r lang lat by fMRI
#gL=numpy.array([ 3,5,  6,  8, 11, 15, 16, 19, 21, 23, 24])

# 3 groups
thresh = 60
gL3 = numpy.nonzero(handedness < -thresh)[0]
gIC3 = numpy.nonzero(numpy.abs(handedness) <= thresh)[0]
gR3 = numpy.nonzero(handedness > thresh)[0]

thresh = 75
gCH2 = numpy.nonzero(numpy.abs(handedness) >= thresh)[0]
gIC2 = numpy.nonzero(numpy.abs(handedness) < thresh)[0]

ksize = 3
ksize = 4
kernel=numpy.ones([ksize, ksize, ksize])
kernel = kernel/numpy.sum(kernel)
smooth_symm_LIL = scipy.signal.convolve(symm_LIL, kernel, 'same')
smooth_symm_LIR = scipy.signal.convolve(symm_LIR, kernel, 'same')
image_movie(smooth_symm_LIR+smooth_symm_LIL, 50, 150, 5)

# test summing all
sum_counts_groups = []
sum_LIR_groups = []
sum_LIL_groups = []

#groups = [gL3, gIC3, gR3]
groups = [gL,gR]

for gidx in range(len(groups)):
    sum_counts_groups.append(numpy.zeros(data_fiber_counts[0].shape))
    sum_LIR_groups.append(numpy.zeros(data_fiber_counts[0].shape))
    sum_LIL_groups.append(numpy.zeros(data_fiber_counts[0].shape))
    g = groups[gidx]
    for sidx in g:
        sum_counts_groups[gidx] += (data_fiber_counts[sidx] > 0)
        # LI average across fibers, so from -1 to 1
        # only do this division where there is a fiber
        idx = numpy.nonzero(data_fiber_counts[sidx])
        LI = numpy.zeros(data_fiber_counts[0].shape)
        LI[idx] = numpy.divide(data_fiber_LI[sidx][idx], data_fiber_counts[sidx][idx])
        sum_LIR_groups[gidx][LI > 0] += LI[LI > 0]
        sum_LIL_groups[gidx][LI < 0] += LI[LI < 0]

smooth_LIL_groups = []
smooth_LIR_groups = []
smooth_counts_groups = []
for gidx in range(len(groups)):
    smooth_LIL_groups.append(scipy.signal.convolve(sum_LIL_groups[gidx], kernel, 'same'))
    smooth_LIR_groups.append(scipy.signal.convolve(sum_LIR_groups[gidx], kernel, 'same'))
    smooth_counts_groups.append(scipy.signal.convolve(sum_counts_groups[gidx], kernel, 'same'))

rgbL = numpy.abs(numpy.array([smooth_LIL_groups[0].T,smooth_LIL_groups[1].T, smooth_LIL_groups[2].T]).T)
rgbR = numpy.array([smooth_LIR_groups[0].T,smooth_LIR_groups[1].T, smooth_LIR_groups[2].T]).T
rgb_counts = numpy.array([smooth_counts_groups[0].T,smooth_counts_groups[1].T, smooth_counts_groups[2].T]).T

#image_movie(smooth_LIR_groups[0]+smooth_LIL_groups[0], 50, 150, 5)
#image_movie(smooth_LIR_groups[1]+smooth_LIL_groups[1], 50, 150, 5)
#image_movie(smooth_LIR_groups[2]+smooth_LIL_groups[2], 50, 150, 5)

#image_movie(sum_LIR_groups[0]+sum_LIL_groups[0], 50, 150, 5)
#image_movie(sum_LIR_groups[1]+sum_LIL_groups[1], 50, 150, 5)
#image_movie(sum_LIR_groups[2]+sum_LIL_groups[2], 50, 150, 5)

#    flip_sum_LIL = sum_LIL[::-1,:,:]
#    symm_LIL = sum_LIL + flip_sum_LIL

#        flip_sum_LIR = sum_LIR[::-1,:,:]
#        symm_LIR = sum_LIR + flip_sum_LIR

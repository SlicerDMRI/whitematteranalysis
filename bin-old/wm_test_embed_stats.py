import glob

import scipy.stats
import numpy
import matplotlib.pyplot as plt

import whitematteranalysis as wma

#inputDirectory = '/Users/lauren/Desktop/OUTPUT'
#inputDirectory = '/Users/lauren/Dropbox/Coding/OUTPUTS/laterality_test_mp'
inputDirectory = '/Users/odonnell/Desktop/OUTPUT/laterality_test_register'
inputMask = "{0}/*.vtk".format(inputDirectory);
inputs = glob.glob(inputMask)

inputs = inputs[0]
print inputs


data = []
# Read data, and check if the data are all present
for sidx in range(0,len(inputs)) :

    fname = inputs[sidx]
    try:
        pd = wma.io.read_polydata(fname)
        data.append(pd)
    except:
        print "unexpected error" , fname
    else:
        print "Read success", fname


# try to embed.
#appender = vtk.vtkAppendPolyData()
#for pd in input_pds:
#    appender.AddInput(pd)
#appender.Update()
#pd_all_registered = appender.GetOutput()

pd_c, clusters, colors, embed, distortion = wma.cluster.spectral(pd_all_registered,number_of_jobs=number_of_jobs)


# grab first 3 components of embedding

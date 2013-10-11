import glob

import scipy.stats
import numpy
import matplotlib.pyplot as plt

import whitematteranalysis as wma

#inputDirectory = '/Users/lauren/Desktop/OUTPUT'
inputDirectory = '/Users/odonnell/Desktop/OUTPUT/laterality_test_distsq'
#'/Users/lauren/Dropbox/Coding/OUTPUTS/laterality_test_mp'
inputMask1 = "{0}/*".format(inputDirectory);
inputDirs = glob.glob(inputMask1)

# Read data, and render
for sidx in range(0,len(inputDirs)) :

    dir = inputDirs[sidx]
    try:
        resultIO = wma.io.LateralityResults()
        resultIO.read(dir, readpd=True)
        ren = wma.render.render(resultIO.polydata)
        ren.save_views(dir)
        del ren
        del resultIO        
    except:
        print "unexpected error" , dir
    else:
        print "Read success", dir

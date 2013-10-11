import glob

import scipy.stats
import numpy
import matplotlib.pyplot as plt
import scipy.signal

import whitematteranalysis as wma

#inputDirectory = '/Users/lauren/Desktop/OUTPUT'
#inputDirectory = '/Users/lauren/Dropbox/Coding/OUTPUTS/laterality_test_mp'
#inputDirectory = '/Users/odonnell/Desktop/OUTPUT/laterality_test_mp_allfibers'
#inputDirectory = '/Users/odonnell/Desktop/OUTPUT/laterality_test_distsq'
#inputDirectory = '/Users/lauren/Dropbox/Coding/OUTPUTS/laterality_all_fibers/dist_sq_withhemtot'
inputDirectory = '/Users/lauren/Dropbox/Coding/OUTPUTS/laterality_all_fibers/dist_sq_constnumfibers'
inputMask1 = "{0}/*".format(inputDirectory);
inputDirs = glob.glob(inputMask1)


#readPolyDatas = 0
def robust_skewness(data):
    # q3 + q1 - 2Q2 / q3 - q1
    # from http://weber.ucsd.edu/~hwhite/pub_files/hwcv-092.pdf
    integral = []
    for r in numpy.arange(0,0.5,0.005):
        q1 = scipy.stats.scoreatpercentile(data, (1-r)*100)
        q2 = scipy.stats.scoreatpercentile(data, 50)
        q3 = scipy.stats.scoreatpercentile(data, r*100)
        integral.append(numpy.divide((q3 + q1 -2*q2), (q3 - q1)))
    return numpy.sum(numpy.array(integral))

def robust_kurtosis(data):
    # e7-e5 + e3-e1 / e6-e2
    e1 = scipy.stats.scoreatpercentile(data, 12.5)
    e2 = scipy.stats.scoreatpercentile(data, 25)
    e3 = scipy.stats.scoreatpercentile(data, 37.5)
    #e4 = scipy.stats.scoreatpercentile(data, 50)
    e5 = scipy.stats.scoreatpercentile(data, 62.5)
    e6 = scipy.stats.scoreatpercentile(data, 75)
    e7 = scipy.stats.scoreatpercentile(data, 87.5)
    #e8 = scipy.stats.scoreatpercentile(data, 1)
    return(numpy.divide((e7 - e5 + e3 - e1), (e6 - e2)))

    
data = []
# Read data, and check if the data are all present
for sidx in range(0,len(inputDirs)):
    dir = inputDirs[sidx]
    try:
        resultIO = wma.io.LateralityResults()
        resultIO.read(dir)
        data.append(resultIO)
    except:
        print "unexpected error" , dir
    else:
        print "Read success", dir
        print resultIO.laterality_index[0]

# compute summary stats
median = numpy.zeros(len(data))
iqr = numpy.zeros(len(data))
skew = numpy.zeros(len(data))
kurtosis = numpy.zeros(len(data))
num_fibers = numpy.zeros(len(data))
num_fibers_hem = numpy.zeros(len(data))
num_fibers_comm = numpy.zeros(len(data))

for sidx in range(0,len(data)):
    print sidx
    print data[sidx].directory
    LI = data[sidx].laterality_index
    confidence = numpy.divide(data[sidx].left_hem_similarity + data[sidx].right_hem_similarity, 2)
    confidence = numpy.minimum(data[sidx].left_hem_similarity, data[sidx].right_hem_similarity)

    # test define new LI that does not have a bias toward hemisphere fiber is in
    #lhem = numpy.divide(data[sidx].left_hem_similarity, numpy.max(data[sidx].left_hem_similarity))
    #rhem = numpy.divide(data[sidx].right_hem_similarity, numpy.max(data[sidx].right_hem_similarity))    
    #lhem = numpy.divide(data[sidx].left_hem_similarity, numpy.median(data[sidx].left_hem_similarity))
    #rhem = numpy.divide(data[sidx].right_hem_similarity, numpy.median(data[sidx].right_hem_similarity))
    #lhem = numpy.divide(data[sidx].left_hem_similarity, numpy.mean(data[sidx].left_hem_similarity))
    #rhem = numpy.divide(data[sidx].right_hem_similarity, numpy.mean(data[sidx].right_hem_similarity))        
    #lhem = data[sidx].left_hem_similarity - numpy.min(data[sidx].left_hem_similarity)
    #rhem = data[sidx].right_hem_similarity - numpy.min(data[sidx].right_hem_similarity)
    #lhem = data[sidx].left_hem_similarity
    #rhem = data[sidx].right_hem_similarity
    #numerator = rhem - lhem
    #denominator = rhem + lhem
    #denominator = denominator + 0.5*numpy.mean(denominator)
    #maxdiff = numpy.max(numpy.abs(numerator))
    #numerator[numpy.abs(numerator) < 0.01*maxdiff] = 0.0
    #LI = numpy.divide(numerator, denominator)
    #LI = numerator
    
    #LI =confidence
    #data[sidx].laterality_index=confidence
    #confidence = numpy.maximum(data[sidx].left_hem_similarity, data[sidx].right_hem_similarity)
    confidence_NZ = confidence[numpy.nonzero(LI)]
    #confidence_NZ = confidence_NZ / numpy.max(confidence_NZ)
    LI_NZ = LI[numpy.nonzero(LI)]
    #plt.figure()
    #plt.plot(confidence_NZ, LI_NZ, 'o')
    # throw out fibers in tiny structures, with low total similarity. these can cause an outlier
    #select = confidence_NZ > scipy.stats.scoreatpercentile(confidence_NZ, 2.5)
    select = confidence_NZ > scipy.stats.scoreatpercentile(confidence_NZ, 5)
    #select = confidence_NZ > scipy.stats.scoreatpercentile(confidence_NZ, 75)
    #select = confidence_NZ > scipy.stats.scoreatpercentile(confidence_NZ, 1)

    #select = confidence_NZ > scipy.stats.scoreatpercentile(confidence_NZ, 50)
    #select = confidence_NZ > 550
    #select = confidence_NZ > 200
    #select = numpy.multiply((confidence_NZ > 400) , (confidence_NZ <800 ))
    #select = numpy.multiply((confidence_NZ > 0.5) , (confidence_NZ < 0.7))
    #select = numpy.multiply(select, confidence_NZ < scipy.stats.scoreatpercentile(confidence_NZ, 97.5))
    # look at just left hem fibers
    #select = numpy.nonzero((data[sidx].hemisphere[numpy.nonzero(LI)] == -1))[0]
    trimmed_LI_NZ = LI_NZ[select]
    LI_NZ = trimmed_LI_NZ
    data[sidx].laterality_index = LI_NZ
    median[sidx] = numpy.median(LI_NZ)
    print sidx, ":" , median[sidx]
    iqr[sidx] = scipy.stats.scoreatpercentile(LI_NZ, 75) - scipy.stats.scoreatpercentile(LI_NZ, 25)
    #skew[sidx] = scipy.stats.skew(LI_NZ)
    skew[sidx] = robust_skewness(LI_NZ)
    #kurtosis[sidx] = scipy.stats.kurtosis(LI_NZ)
    kurtosis[sidx] = robust_kurtosis(LI_NZ)
    num_fibers[sidx] = len(LI)
    num_fibers_hem[sidx] = len(LI_NZ)

num_fibers_comm =  num_fibers - num_fibers_hem

    
# Read groups file, in theory
handedness = numpy.array([75, 85, 55, -35, -90, -85, -95, 70, -75, -25, 60, -70, 75, 100, 75, -35, -45, 90, 85, -75, 5, -80, 50, -40, -95, 100])
fMRIlat  = numpy.array([0.25, 0.388, 0.291, 0.1885, -0.1605, 0.4435, 0.6385, 0.0455, 0.3805, -0.358, 0.6595, 0.3045, 0.537, 0.236, 0.202, 0.321, 0.4405, 0.038, 0.435, 0.304, 0.1385, 0.667, 0.337, 0.699, 0.31, 0.3925])

gL = numpy.nonzero(handedness < 0)[0]
gR = numpy.nonzero(handedness >= 0)[0]

plt.figure()
for sidx in gL:
    LI = data[sidx].laterality_index
    confidence = numpy.minimum(data[sidx].left_hem_similarity, data[sidx].left_hem_similarity)
    confidence_NZ = confidence[numpy.nonzero(LI)]
    confidence_NZ = confidence_NZ / numpy.max(confidence_NZ)
    LI_NZ = LI[numpy.nonzero(LI)]
    select = confidence_NZ>200
    plt.plot(confidence_NZ, LI_NZ, 'bo')
    #plt.plot(confidence_NZ[select], LI_NZ[select], 'yo')

plt.title('LH')

plt.figure()
for sidx in gR:
    LI = data[sidx].laterality_index
    confidence = numpy.minimum(data[sidx].left_hem_similarity, data[sidx].right_hem_similarity)
    confidence_NZ = confidence[numpy.nonzero(LI)]
    confidence_NZ = confidence_NZ / numpy.max(confidence_NZ)
    LI_NZ = LI[numpy.nonzero(LI)]
    plt.plot(confidence_NZ, LI_NZ, 'ro')

plt.title('RH')

# exclude 21
#subjects = [range(0,21) + range(22,26)]
subjects = range(0,26)

varlist = [[median, "median"], [iqr, "iqr"], [skew, "skew"], [kurtosis, "kurtosis"]]

print "handedness"
for var in varlist:
    plt.figure()
    plt.plot(handedness[subjects], var[0][subjects],'o')
    (a_s, b_s, r, p, stderr) = scipy.stats.linregress(handedness[subjects], var[0][subjects])
    xr=numpy.polyval([a_s,b_s],handedness[subjects])
    plt.plot(handedness[subjects],xr)
    plt.title(var[1])
    #r,p = scipy.stats.pearsonr(handedness[subjects], var[0][subjects])
    print var[1],":", r, p
    
print "absolute fMRI lat"
for var in varlist:
    plt.figure()
    plt.plot(fMRIlat[subjects], var[0][subjects],'yo')
    plt.plot(numpy.abs(fMRIlat[subjects]), var[0][subjects],'o')
    (a_s, b_s, r, p, stderr) = scipy.stats.linregress(numpy.abs(fMRIlat[subjects]), var[0][subjects])
    xr=numpy.polyval([a_s,b_s],numpy.abs(fMRIlat[subjects]))
    plt.plot(numpy.abs(fMRIlat[subjects]),xr)
    plt.title(var[1])
    #r,p = scipy.stats.pearsonr(numpy.abs(fMRIlat[subjects]), var[0][subjects])
    print var[1],":", r, p

print "fMRI lat"
for var in varlist:
    plt.figure()
    # show left lateralized in yellow
    plt.plot(fMRIlat[subjects], var[0][subjects],'yo')
    plt.plot(numpy.abs(fMRIlat[subjects]), var[0][subjects],'o')
    (a_s, b_s, r, p, stderr) = scipy.stats.linregress(fMRIlat[subjects], var[0][subjects])
    xr=numpy.polyval([a_s,b_s],fMRIlat[subjects])
    plt.plot(fMRIlat[subjects],xr)
    plt.title(var[1])
    #r,p = scipy.stats.pearsonr(fMRIlat[subjects], var[0][subjects])
    print var[1],":", r, p


histograms = list()
for sidx in range(0,len(data)) :
    print sidx
    LI = data[sidx].laterality_index
    LI_NZ = LI[numpy.nonzero(LI)]
    #trimmed_LI_NZ = numpy.sort(LI_NZ)[50:-50]
    #LI_NZ = trimmed_LI_NZ
    #ret = plt.hist(LI_NZ,range=[-1,1],normed=True,bins=20)
    #ret = plt.hist(LI_NZ,range=[-1,1],normed=True,bins=60)
    ret = plt.hist(LI_NZ,range=[-1,1],normed=True,bins=40)
    histograms.append(ret[0])
    bins = ret[1]
    plt.close()


# 3 groups (thresh 60 divides into 3 almost even groups 9,9,8)
thresh = 60
gL3 = numpy.nonzero(handedness < -thresh)[0]
gIC3 = numpy.nonzero(numpy.abs(handedness) <= thresh)[0]
gR3 = numpy.nonzero(handedness > thresh)[0]

histograms_smooth=[]
#kernel = scipy.signal.gaussian(20,5)
kernel = scipy.signal.gaussian(10,1.5)
kernel = scipy.signal.gaussian(10,0.5)
kernel = kernel/numpy.sum(kernel)

for sidx in range(0,len(data)):
    histograms_smooth.append(scipy.signal.convolve(histograms[sidx], kernel, mode='same'))

plt.figure()
for sidx in gL3:
    plt.plot(histograms_smooth[sidx])

plt.title('L')
plt.savefig('group_L_histograms.pdf')
plt.close()

plt.figure()
for sidx in gIC3:
    plt.plot(histograms_smooth[sidx])

plt.title('IC')
plt.savefig('group_IC_histograms.pdf')
plt.close()

plt.figure()
for sidx in gR3:
    plt.plot(histograms_smooth[sidx])

plt.title('R')
plt.savefig('group_R_histograms.pdf')
plt.close()




plt.figure()

hist_sum = numpy.zeros(len(histograms_smooth[0]))
for sidx in gL3:
    #print sidx
    #plt.plot(histograms_smooth[sidx])
    hist_sum = hist_sum + histograms_smooth[sidx]

hist_sum = numpy.divide(hist_sum, len(gL3))
plt.plot(hist_sum,'b')

hist_sum = numpy.zeros(len(histograms_smooth[0]))
for sidx in gIC3:
    #print sidx
    #plt.plot(histograms_smooth[sidx])
    hist_sum = hist_sum + histograms_smooth[sidx]

hist_sum = numpy.divide(hist_sum, len(gIC3))
plt.plot(hist_sum,'m')

hist_sum = numpy.zeros(len(histograms_smooth[0]))
for sidx in gR3:
    #print sidx
    #plt.plot(histograms_smooth[sidx])
    hist_sum = hist_sum + histograms_smooth[sidx]

hist_sum = numpy.divide(hist_sum, len(gR3))

plt.plot(hist_sum,'r')
plt.title('CL-ICH-CR')
plt.savefig('groups_mean_histogram.pdf')
plt.close()


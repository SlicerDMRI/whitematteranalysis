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
inputDirectory = '/Users/lauren/Dropbox/Coding/OUTPUTS/laterality_all_fibers/dist_sq_withhemtot'
#inputDirectory = '/Users/lauren/Dropbox/Coding/OUTPUTS/laterality_all_fibers/dist_sq_constnumfibers'
inputMask1 = "{0}/*".format(inputDirectory);
inputDirs = glob.glob(inputMask1)


#readPolyDatas = 0

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
trimmed_skew = numpy.zeros(len(data))
num_fibers = numpy.zeros(len(data))
num_fibers_hem = numpy.zeros(len(data))
num_fibers_comm = numpy.zeros(len(data))

for sidx in range(0,len(data)):
    print sidx
    print data[sidx].directory
    LI = data[sidx].laterality_index
    confidence = numpy.minimum(data[sidx].left_hem_similarity, data[sidx].right_hem_similarity)
    confidence_NZ = confidence[numpy.nonzero(LI)]
    LI_NZ = LI[numpy.nonzero(LI)]
    plt.figure()
    plt.plot(confidence_NZ, LI_NZ, 'o')
    #mode = scipy.stats.mode(numpy.array(LI_NZ*1000).astype(int))[0]/1000
    #print mode
    #print numpy.median(LI_NZ)
    #LI_NZ =  LI_NZ - mode
    
    #trimmed_LI_NZ = numpy.sort(LI_NZ)[100:-100]
    #trimmed_LI_NZ = numpy.sort(LI_NZ)[50:-50]
    # throw out fibers in tiny structures, with low total similarity. these can cause an outlier
    select = confidence_NZ > scipy.stats.scoreatpercentile(confidence_NZ, 2.5)
    select = confidence_NZ > scipy.stats.scoreatpercentile(confidence_NZ, 1)
    trimmed_LI_NZ = LI_NZ[select]
    #LI_NZ = trimmed_LI_NZ 
    median[sidx] = numpy.median(LI_NZ)
    print sidx, ":" , median[sidx]
    iqr[sidx] = scipy.stats.scoreatpercentile(LI_NZ, 75) - scipy.stats.scoreatpercentile(LI_NZ, 25)
    skew[sidx] = scipy.stats.skew(LI_NZ)
    trimmed_skew[sidx] = scipy.stats.skew(trimmed_LI_NZ)
    kurtosis[sidx] = scipy.stats.kurtosis(LI_NZ)
    #trimmed_kurtosis[sidx] = scipy.stats.kurtosis(trimmed_LI_NZ)
    num_fibers[sidx] = len(LI)
    num_fibers_hem[sidx] = len(LI_NZ)

num_fibers_comm =  num_fibers - num_fibers_hem
    
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

    
# Read groups file, in theory
handedness = numpy.array([75, 85, 55, -35, -90, -85, -95, 70, -75, -25, 60, -70, 75, 100, 75, -35, -45, 90, 85, -75, 5, -80, 50, -40, -95, 100])
fMRIlat  = numpy.array([0.25, 0.388, 0.291, 0.1885, -0.1605, 0.4435, 0.6385, 0.0455, 0.3805, -0.358, 0.6595, 0.3045, 0.537, 0.236, 0.202, 0.321, 0.4405, 0.038, 0.435, 0.304, 0.1385, 0.667, 0.337, 0.699, 0.31, 0.3925])

gL = numpy.nonzero(handedness < 0)[0]
gR = numpy.nonzero(handedness >= 0)[0]

# test remove r lang lat by fMRI
#gL=numpy.array([ 3,5,  6,  8, 11, 15, 16, 19, 21, 23, 24])

# 3 groups (thresh 60 divides into 3 almost even groups 9,9,8)
thresh = 60
gL3 = numpy.nonzero(handedness < -thresh)[0]
gIC3 = numpy.nonzero(numpy.abs(handedness) <= thresh)[0]
gR3 = numpy.nonzero(handedness > thresh)[0]

# two groups (thresh 75 divides into 2 even groups)
thresh = 75
gCH2 = numpy.nonzero(numpy.abs(handedness) >= thresh)[0]
gIC2 = numpy.nonzero(numpy.abs(handedness) < thresh)[0]

# significance
# ANOVA and TTEST
p_vals_LR = list()
p_vals_LIR = list()
p_vals_IC = list()
values_tested = list()

plt.close('all')
f = median
values_tested.append('MEDIAN')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plt.boxplot((f[gL3],f[gIC3],f[gR3]));
plt.savefig('group_3_median.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plt.boxplot((f[gL],f[gR]));
plt.savefig('group_2_median.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plt.boxplot((f[gIC2],f[gCH2]));
plt.savefig('group_2mag_median.pdf')
plt.close()

f = iqr
values_tested.append('IQR')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plt.boxplot((f[gL3],f[gIC3],f[gR3]));
plt.savefig('group_3_iqr.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plt.boxplot((f[gL],f[gR]));
plt.savefig('group_2_iqr.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plt.boxplot((f[gIC2],f[gCH2]));
plt.savefig('group_2mag_iqr.pdf')
plt.close()

f = skew
values_tested.append('SKEW')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plt.boxplot((f[gL3],f[gIC3],f[gR3]));
plt.savefig('group_3_skew.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plt.boxplot((f[gL],f[gR]));
plt.savefig('group_2_skew.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plt.boxplot((f[gIC2],f[gCH2]));
plt.savefig('group_2mag_skew.pdf')
plt.close()

f = kurtosis
values_tested.append('KURTOSIS')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plt.boxplot((f[gL3],f[gIC3],f[gR3]));
plt.savefig('group_3_kurtosis.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plt.boxplot((f[gL],f[gR]));
plt.savefig('group_2_kurtosis.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plt.boxplot((f[gIC2],f[gCH2]));
plt.savefig('group_2mag_kurtosis.pdf')
plt.close()

f = num_fibers
values_tested.append('NUM_FIBERS')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plt.boxplot((f[gL3],f[gIC3],f[gR3]));
plt.savefig('group_3_num_fibers.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plt.boxplot((f[gL],f[gR]));
plt.savefig('group_2_num_fibers.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plt.boxplot((f[gIC2],f[gCH2]));
plt.savefig('group_2mag_num_fibers.pdf')
plt.close()


f = num_fibers_comm
values_tested.append('NUM_FIBERS_COMM')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plt.boxplot((f[gL3],f[gIC3],f[gR3]));
plt.savefig('group_3_num_fibers_comm.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plt.boxplot((f[gL],f[gR]));
plt.savefig('group_2_num_fibers_comm.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plt.boxplot((f[gIC2],f[gCH2]));
plt.savefig('group_2mag_num_fibers_comm.pdf')
plt.close()

f = num_fibers_hem
values_tested.append('NUM_FIBERS_HEM')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plt.boxplot((f[gL3],f[gIC3],f[gR3]));
plt.savefig('group_3_num_fibers_hem.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plt.boxplot((f[gL],f[gR]));
plt.savefig('group_2_num_fibers_hem.pdf')
plt.close()
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plt.boxplot((f[gIC2],f[gCH2]));
plt.savefig('group_2mag_num_fibers_hem.pdf')
plt.close()

print values_tested
print p_vals_LIR
print p_vals_LR
print p_vals_IC


plt.figure()

hist_sum = numpy.zeros(len(histograms[0]))
for sidx in gL:
    #print sidx
    #plt.plot(histograms[sidx])
    hist_sum = hist_sum + histograms[sidx]

hist_sum = numpy.divide(hist_sum, len(gL))

plt.plot(hist_sum,'b')
#plt.title('Left')

hist_sum = numpy.zeros(len(histograms[0]))
for sidx in gR:
    #print sidx
    #plt.plot(histograms[sidx])
    hist_sum = hist_sum + histograms[sidx]

hist_sum = numpy.divide(hist_sum, len(gR))

plt.plot(hist_sum,'r')
plt.title('Left-Right')



plt.figure()

hist_sum = numpy.zeros(len(histograms[0]))
for sidx in gIC2:
    #print sidx
    #plt.plot(histograms[sidx])
    hist_sum = hist_sum + histograms[sidx]

hist_sum = numpy.divide(hist_sum, len(gIC2))

plt.plot(hist_sum,'m')
#plt.title('Left')

hist_sum = numpy.zeros(len(histograms[0]))
for sidx in gCH2:
    #print sidx
    #plt.plot(histograms[sidx])
    hist_sum = hist_sum + histograms[sidx]

hist_sum = numpy.divide(hist_sum, len(gCH2))

plt.plot(hist_sum,'g')
plt.title('C-ICH')



plt.figure()

hist_sum = numpy.zeros(len(histograms[0]))
for sidx in gL3:
    #print sidx
    #plt.plot(histograms[sidx])
    hist_sum = hist_sum + histograms[sidx]

hist_sum = numpy.divide(hist_sum, len(gL3))
plt.plot(hist_sum,'b')

hist_sum = numpy.zeros(len(histograms[0]))
for sidx in gIC3:
    #print sidx
    #plt.plot(histograms[sidx])
    hist_sum = hist_sum + histograms[sidx]

hist_sum = numpy.divide(hist_sum, len(gIC3))
plt.plot(hist_sum,'m')

hist_sum = numpy.zeros(len(histograms[0]))
for sidx in gR3:
    #print sidx
    #plt.plot(histograms[sidx])
    hist_sum = hist_sum + histograms[sidx]

hist_sum = numpy.divide(hist_sum, len(gR3))

plt.plot(hist_sum,'r')
plt.title('CL-ICH-CR')

fsizeFA3 = [187805503, 199041963, 218436703, 223583531, 153971215,
172693903, 148755955, 188697603, 206594871, 197826243, 174418667,
179026003, 174238755, 237032887, 188905851, 156799379, 203734623,
180085291, 273057119, 149181543, 209641047, 207737779, 249086479,
149728527, 198696883, 261106647]

fsizeFA2 = [197327763, 209591923, 230822323, 233944635, 162380443, 182351623,
157801035, 199051019, 219027175, 208694231, 185016615, 189349551,
182872007, 247033963, 200255167, 164361735, 215349343, 189788671,
286057503, 159945623, 220542463, 218571735, 260058123, 157069799,
206887131, 273348987]


varlist = [[median, "median"], [iqr, "iqr"], [skew, "skew"], [kurtosis, "kurtosis"], [num_fibers, "num_fibers"], [num_fibers_comm, "num_fibers_comm"], [num_fibers_hem, "num_fibers_hem"], [fsizeFA3, "fsizeFA3"], [fsizeFA2, "fsizeFA2"], [num_fibers_comm/num_fibers, "comm percent"], [num_fibers_hem/num_fibers, "hem percent"]]

for var in varlist:
    plt.figure()
    plt.plot(handedness, var[0],'o')
    (a_s, b_s, r, p, stderr) = scipy.stats.linregress(handedness, var[0])
    xr=numpy.polyval([a_s,b_s],handedness)
    plt.plot(handedness,xr)
    plt.title(var[1])
    #r,p = scipy.stats.pearsonr(handedness, var[0])
    print var[1],":", r, p


histograms_smooth=[]
#kernel = scipy.signal.gaussian(20,5)
kernel = scipy.signal.gaussian(10,1.5)
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

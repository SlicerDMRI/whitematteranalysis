import glob

import scipy.stats
import numpy
import matplotlib.pyplot as plt
import scipy.signal

import whitematteranalysis as wma

inputDirectory = '../.'

inputMask1 = "{0}/*".format(inputDirectory);
inputDirs = glob.glob(inputMask1)


def plot_box(data, labels, fname):
    f = plt.figure()
    plt.boxplot(data);
    a = f.get_axes()
    plt.setp(a,xticklabels=labels)
    plt.savefig(fname)
    plt.close()
    
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
num_fibers_Lhem = numpy.zeros(len(data))
num_fibers_Rhem = numpy.zeros(len(data))
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

    num_fibers_Rhem[sidx] = len(numpy.nonzero(numpy.array(data[sidx].hemisphere == 1))[0])
    num_fibers_Lhem[sidx] = len(numpy.nonzero(numpy.array(data[sidx].hemisphere == -1))[0])
    
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
###handedness = numpy.array([75, 85, 55, -35, -90, -85, -95, 70, -75, -25, 60, -70, 75, 100, 75, -35, -45, 90, 85, -75, 5, -80, 50, -40, -95, 100])
###fMRIlat  = numpy.array([0.25, 0.388, 0.291, 0.1885, -0.1605, 0.4435, 0.6385, 0.0455, 0.3805, -0.358, 0.6595, 0.3045, 0.537, 0.236, 0.202, 0.321, 0.4405, 0.038, 0.435, 0.304, 0.1385, 0.667, 0.337, 0.699, 0.31, 0.3925])
# add new data
handedness = numpy.array([75, 85, 55, -35, -90, -85, -95, 70, -75, -25, 60, -70, 75, 100, 75, -35, -45, 90, 85, -75, 5, -80, 50, -40, -95, 100, -80, -75, -70])
# female=1 male=2
gender = numpy.array([1, 2, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 2,   2, 2, 2])

IDs = numpy.array(['AG_5114', 'AG_5123', 'AG_5125', 'AG_5129', 'AG_5131', 'AG_5132', 'AG_5138', 'AG_5139', 'AG_5142', 'AG_5146', 'AG_5147', 'AG_5148', 'AG_5149', 'AG_5150', 'AG_5151', 'AG_5152', 'AG_5153', 'AG_5154', 'AG_5155', 'AG_5156', 'AG_5158', 'AG_5159', 'AG_5160', 'AG_5162', 'AG_5163', 'AG_5205', 'AG_5263', 'AG_5266', 'AG_5268']) 


# still need to re-run entire pipeline with optimal parameters and max # fibers possible
# if need to exclude
#analyze = (IDs != 'AG_5159') & (IDs != 'AG_5156')
#gender = gender[numpy.nonzero(analyze)]
#handedness = handedness[numpy.nonzero(analyze)]

#gL = numpy.nonzero(handedness < 0)[0]
#gR = numpy.nonzero(handedness >= 0)[0]

# --------------
threshold = 60
#threshold = 40
# --------------
gL = numpy.nonzero(handedness <= threshold)[0]
gR = numpy.nonzero(handedness >= threshold)[0]

# test remove r lang lat by fMRI
#gL=numpy.array([ 3,5,  6,  8, 11, 15, 16, 19, 21, 23, 24])

##### 3 groups (thresh 60 divides into 3 almost even groups 9,9,8)
###thresh = 60
#thresh = 50
thresh = threshold
gL3 = numpy.nonzero(handedness < -thresh)[0]
gIC3 = numpy.nonzero(numpy.abs(handedness) <= thresh)[0]
gR3 = numpy.nonzero(handedness > thresh)[0]

# two groups (thresh 75 divides into 2 even groups)
#thresh = 75
gCH2 = numpy.nonzero(numpy.abs(handedness) >= thresh)[0]
gIC2 = numpy.nonzero(numpy.abs(handedness) < thresh)[0]

# any gender differences
gF = numpy.nonzero(gender == 1)
gM = numpy.nonzero(gender == 2)


# significance
# ANOVA and TTEST
p_vals_LR = list()
p_vals_LIR = list()
p_vals_IC = list()
p_vals_MF = list()
values_tested = list()

plt.close('all')
f = median
values_tested.append('MEDIAN')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plot_box((f[gL3],f[gIC3],f[gR3]),['L', 'IC', 'R'],'group_3_median.pdf')
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plot_box((f[gL],f[gR]), ('L', 'R'), 'group_2_median.pdf')
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plot_box((f[gIC2],f[gCH2]), ['ICH', 'CH'], 'group_2mag_median.pdf')
[s, p] = scipy.stats.ttest_ind(f[gM],f[gF])
p_vals_MF.append(p)
plot_box((f[gM],f[gF]), ['M', 'F'], 'group_2gen_median.pdf')

f = iqr
values_tested.append('IQR')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plot_box((f[gL3],f[gIC3],f[gR3]), ['L', 'IC', 'R'], 'group_3_iqr.pdf')
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plot_box((f[gL],f[gR]), ('L', 'R'), 'group_2_iqr.pdf')
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plot_box((f[gIC2],f[gCH2]), ['ICH', 'CH'], 'group_2mag_iqr.pdf')
[s, p] = scipy.stats.ttest_ind(f[gM],f[gF])
p_vals_MF.append(p)
plot_box((f[gM],f[gF]), ['M', 'F'], 'group_2gen_iqr.pdf')

f = skew
values_tested.append('SKEW')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plot_box((f[gL3],f[gIC3],f[gR3]), ['L', 'IC', 'R'], 'group_3_skew.pdf')
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plot_box((f[gL],f[gR]), ('L', 'R'), 'group_2_skew.pdf')
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plot_box((f[gIC2],f[gCH2]), ['ICH', 'CH'], 'group_2mag_skew.pdf')
[s, p] = scipy.stats.ttest_ind(f[gM],f[gF])
p_vals_MF.append(p)
plot_box((f[gM],f[gF]), ['M', 'F'], 'group_2gen_skew.pdf')

f = kurtosis
values_tested.append('KURTOSIS')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plot_box((f[gL3],f[gIC3],f[gR3]), ['L', 'IC', 'R'], 'group_3_kurtosis.pdf')
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plot_box((f[gL],f[gR]), ('L', 'R'), 'group_2_kurtosis.pdf')
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plot_box((f[gIC2],f[gCH2]), ['ICH', 'CH'], 'group_2mag_kurtosis.pdf')
[s, p] = scipy.stats.ttest_ind(f[gM],f[gF])
p_vals_MF.append(p)
plot_box((f[gM],f[gF]), ['M', 'F'], 'group_2gen_kurtosis.pdf')


f = num_fibers
values_tested.append('NUM_FIBERS')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plot_box((f[gL3],f[gIC3],f[gR3]), ['L', 'IC', 'R'], 'group_3_num_fibers.pdf')
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plot_box((f[gL],f[gR]), ('L', 'R'), 'group_2_num_fibers.pdf')
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plot_box((f[gIC2],f[gCH2]), ['ICH', 'CH'], 'group_2mag_num_fibers.pdf')
[s, p] = scipy.stats.ttest_ind(f[gM],f[gF])
p_vals_MF.append(p)
plot_box((f[gM],f[gF]), ['M', 'F'], 'group_2gen_num_fibers.pdf')

f = num_fibers_Rhem
values_tested.append('NUM_FIBERS_RHEM')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plot_box((f[gL3],f[gIC3],f[gR3]), ['L', 'IC', 'R'], 'group_3_num_fibers_Rhem.pdf')
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plot_box((f[gL],f[gR]), ('L', 'R'), 'group_2_num_fibers_Rhem.pdf')
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plot_box((f[gIC2],f[gCH2]), ['ICH', 'CH'], 'group_2mag_num_fibers_Rhem.pdf')
[s, p] = scipy.stats.ttest_ind(f[gM],f[gF])
p_vals_MF.append(p)
plot_box((f[gM],f[gF]), ['M', 'F'], 'group_2gen_num_fibers_Rhem.pdf')

f = num_fibers_Lhem
values_tested.append('NUM_FIBERS_LHEM')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plot_box((f[gL3],f[gIC3],f[gR3]), ['L', 'IC', 'R'], 'group_3_num_fibers_Lhem.pdf')
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plot_box((f[gL],f[gR]), ('L', 'R'), 'group_2_num_fibers_Lhem.pdf')
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plot_box((f[gIC2],f[gCH2]), ['ICH', 'CH'], 'group_2mag_num_fibers_Lhem.pdf')
[s, p] = scipy.stats.ttest_ind(f[gM],f[gF])
p_vals_MF.append(p)
plot_box((f[gM],f[gF]), ['M', 'F'], 'group_2gen_num_fibers_Lhem.pdf')

f = num_fibers_comm
values_tested.append('NUM_FIBERS_COMM')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plot_box((f[gL3],f[gIC3],f[gR3]), ['L', 'IC', 'R'], 'group_3_num_fibers_comm.pdf')
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plot_box((f[gL],f[gR]), ('L', 'R'), 'group_2_num_fibers_comm.pdf')
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plot_box((f[gIC2],f[gCH2]), ['ICH', 'CH'], 'group_2mag_num_fibers_comm.pdf')
[s, p] = scipy.stats.ttest_ind(f[gM],f[gF])
p_vals_MF.append(p)
plot_box((f[gM],f[gF]), ['M', 'F'], 'group_2gen_num_fibers_comm.pdf')

f = num_fibers_hem
values_tested.append('NUM_FIBERS_HEM')
[s, p] = scipy.stats.f_oneway(f[gL3],f[gIC3],f[gR3])
p_vals_LIR.append(p)
plot_box((f[gL3],f[gIC3],f[gR3]), ['L', 'IC', 'R'], 'group_3_num_fibers_hem.pdf')
[s, p] = scipy.stats.ttest_ind(f[gL],f[gR])
p_vals_LR.append(p)
plot_box((f[gL],f[gR]), ('L', 'R'), 'group_2_num_fibers_hem.pdf')
[s, p] = scipy.stats.ttest_ind(f[gIC2],f[gCH2])
p_vals_IC.append(p)
plot_box((f[gIC2],f[gCH2]), ['ICH', 'CH'], 'group_2mag_num_fibers_hem.pdf')
[s, p] = scipy.stats.ttest_ind(f[gM],f[gF])
p_vals_MF.append(p)
plot_box((f[gM],f[gF]), ['M', 'F'], 'group_2gen_num_fibers_hem.pdf')

    
print values_tested
print "L vs ICH vs R"
print p_vals_LIR
print "L vs R"
print p_vals_LR
print "CH vs ICH"
print p_vals_IC
print "M vs F"
print p_vals_MF


plt.figure()
hist_sum = numpy.zeros(len(histograms[0]))
for sidx in gL:
    #print sidx
    #plt.plot(histograms[sidx])
    hist_sum = hist_sum + histograms[sidx]
hist_sum = numpy.divide(hist_sum, len(gL))
plt.plot(hist_sum,'b')
hist_sum = numpy.zeros(len(histograms[0]))
for sidx in gR:
    #print sidx
    #plt.plot(histograms[sidx])
    hist_sum = hist_sum + histograms[sidx]
hist_sum = numpy.divide(hist_sum, len(gR))
plt.plot(hist_sum,'r')
plt.title('Left-Right')
plt.legend(['L', 'R'])
plt.savefig('leftplusright_hist.pdf')


plt.figure()
hist_sum = numpy.zeros(len(histograms[0]))
for sidx in gIC2:
    #print sidx
    #plt.plot(histograms[sidx])
    hist_sum = hist_sum + histograms[sidx]
hist_sum = numpy.divide(hist_sum, len(gIC2))
plt.plot(hist_sum,'m')
hist_sum = numpy.zeros(len(histograms[0]))
for sidx in gCH2:
    #print sidx
    #plt.plot(histograms[sidx])
    hist_sum = hist_sum + histograms[sidx]
hist_sum = numpy.divide(hist_sum, len(gCH2))
plt.plot(hist_sum,'g')
plt.legend(['ICH', 'CH'])
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
plt.legend(['CL', 'ICH', 'CR'])

varlist = [[median, "median"], [iqr, "iqr"], [skew, "skew"], [kurtosis, "kurtosis"], [num_fibers, "num_fibers"], [num_fibers_Rhem, "num_fibers_Rhem"], [num_fibers_Lhem, "num_fibers_Lhem"], [num_fibers_comm, "num_fibers_comm"], [num_fibers_hem, "num_fibers_hem"], [num_fibers_comm/num_fibers, "commissural percent"], [num_fibers_hem/num_fibers, "hemispheric percent"]]

for var in varlist:
    plt.figure()
    plt.plot(handedness, var[0],'o')
    (a_s, b_s, r, p, stderr) = scipy.stats.linregress(handedness, var[0])
    xr=numpy.polyval([a_s,b_s],handedness)
    plt.plot(handedness,xr)
    plt.title(var[1])
    plt.savefig('handedness'+var[1]+'.pdf')
    r,p = scipy.stats.pearsonr(handedness, var[0])
    print var[1],".. r:", r, " p: ", p


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

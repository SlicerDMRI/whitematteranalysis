import glob

import scipy.stats
import numpy
import matplotlib.pyplot as plt
import scipy.signal

import whitematteranalysis as wma

input_directory = '/Users/lauren/Dropbox/Coding/OUTPUTS/laterality_all_fibers/dist_sq_constnumfibers'
input_mask = "{0}/*".format(input_directory);
input_dirs = glob.glob(input_mask)


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
for sidx in range(0,len(input_dirs)):
    dir = input_dirs[sidx]
    try:
        result_io = wma.io.LateralityResults()
        result_io.read(dir)
        data.append(result_io)
        del result_io
    except:
        print "unexpected error" , dir
    else:
        print "Read success", dir


# compute summary stats
median = numpy.zeros(len(data))
iqr = numpy.zeros(len(data))
skew = numpy.zeros(len(data))
kurtosis = numpy.zeros(len(data))
LI = list()
num_fibers = numpy.zeros(len(data))
num_fibers_hem = numpy.zeros(len(data))
num_fibers_comm = numpy.zeros(len(data))

for sidx in range(0,len(data)):
    print sidx
    print data[sidx].directory
    hem_idx = numpy.nonzero(data[sidx].laterality_index)
    LI.append(data[sidx].laterality_index[hem_idx])
    num_fibers[sidx] = len(data[sidx].laterality_index)
    num_fibers_hem[sidx] = len(hem_idx)
    num_fibers_comm[sidx] =  num_fibers[sidx] - num_fibers_hem[sidx]
    # quality of each fiber is high, if there are many similar fibers
    #confidence = numpy.divide(data[sidx].left_hem_similarity + data[sidx].right_hem_similarity, 2)
    #confidence = numpy.minimum(data[sidx].left_hem_similarity, data[sidx].right_hem_similarity)
    confidence = numpy.maximum(data[sidx].left_hem_similarity[hem_idx], \
                               data[sidx].right_hem_similarity[hem_idx])
    # remove fibers likely to be outliers, make sure they don't drive results
    select = confidence > scipy.stats.scoreatpercentile(confidence, 2.5)
    #select = confidence > scipy.stats.scoreatpercentile(confidence, 50)
    LI[sidx] = LI[sidx][select]
    # summary stats                 
    median[sidx] = numpy.median(LI[sidx])
    iqr[sidx] = scipy.stats.scoreatpercentile(LI[sidx], 75) - scipy.stats.scoreatpercentile(LI[sidx], 25)
    skew[sidx] = scipy.stats.skew(LI[sidx])
    kurtosis[sidx] = scipy.stats.kurtosis(LI[sidx])
    #skew[sidx] = robust_skewness(LI[sidx])
    #kurtosis[sidx] = robust_kurtosis(LI[sidx])

    

# GROUPS
handedness = numpy.array([75, 85, 55, -35, -90, -85, -95, 70, -75, -25, 60, -70, 75, 100, 75, -35, -45, 90, 85, -75, 5, -80, 50, -40, -95, 100])
fMRIlat  = numpy.array([0.25, 0.388, 0.291, 0.1885, -0.1605, 0.4435, 0.6385, 0.0455, 0.3805, -0.358, 0.6595, 0.3045, 0.537, 0.236, 0.202, 0.321, 0.4405, 0.038, 0.435, 0.304, 0.1385, 0.667, 0.337, 0.699, 0.31, 0.3925])
# female=1 male=2
gender = numpy.array([1, 2, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1,
 1, 1, 2, 1, 2, 2])

gL = numpy.nonzero(handedness < 0)[0]
gR = numpy.nonzero(handedness >= 0)[0]

gM = numpy.nonzero(gender == 2)[0]
gF = numpy.nonzero(gender == 1)[0]

#gL = numpy.nonzero(handedness < -50)[0]
#gR = numpy.nonzero(handedness > 50)[0]


# STATS

varlist = [[median, "median"], [iqr, "iqr"], [skew, "skew"], [kurtosis, "kurtosis"]]


print "handedness"
for var in varlist:
    [s, p] = scipy.stats.ttest_ind(var[0][gL],var[0][gR])
    #[s, p] = scipy.stats.ttest_ind(var[0][gLF],var[0][gRF])
    #[s, p] = scipy.stats.ttest_ind(var[0][gLM],var[0][gRM])
    print var[1],":", s, p, "(hand)"
    if p < 0.05:
        plt.figure()
        plt.boxplot((var[0][gL],var[0][gR]));
        plt.title(var[1])
        print scipy.stats.pearsonr(handedness,var[0])
        #plt.savefig('group_2_num_fibers_comm.pdf')
        #plt.close()

print "gender"
for var in varlist:
    [s, p] = scipy.stats.ttest_ind(var[0][gM],var[0][gF])
    print var[1],":", s, p, "(gender)"
    if p < 0.05:
        plt.figure()
        plt.boxplot((var[0][gM],var[0][gF]));
        plt.title(var[1])
        print scipy.stats.pearsonr(gender,var[0])        



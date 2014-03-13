import glob
import os
import argparse

import scipy.stats
import numpy
import matplotlib.pyplot as plt
import scipy.signal

try:
    import whitematteranalysis as wma
except:
    print "<wm_laterality_stats.py> Error importing white matter analysis package\n"
    raise


fname_to_test_if_output = 'tractography_with_LI.vtk'


#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Computes summary statistics about laterality results.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputDirectory',
    help='Contains laterality results, one directory per subject.')
parser.add_argument(
    'outputDirectory',
    help='The output directory should be a new empty directory. It will be created if needed.')
parser.add_argument(
    '-verbose', action='store_true', dest="flag_verbose",
    help='Verbose. Print extra diagnostic information to the screen.')


args = parser.parse_args()

if not os.path.isdir(args.inputDirectory):
    print "<wm_laterality_stats.py> Error: Input directory", args.inputDirectory, "does not exist."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "<wm_laterality_stats.py> Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

print ""
print "--------------------------------------------------------------------"
print "<wm_laterality_stats.py> Analyzing laterality results..."
print "--------------------------------------------------------------------"
print "input directory:\t", args.inputDirectory
print "output directory:\t", args.outputDirectory

input_mask = "{0}/*".format(args.inputDirectory);
inputs = sorted(glob.glob(input_mask))
input_dirs = list()
for dir in inputs:
    if os.path.isdir(dir):
        test_file = os.path.join(dir, fname_to_test_if_output)
        # test this is the output of laterality
        if os.path.exists(test_file):
            input_dirs.append(dir)
            
verbose = args.flag_verbose
print "verbose:\t\t", verbose

# Read groups file, in theory
###handedness = numpy.array([75, 85, 55, -35, -90, -85, -95, 70, -75, -25, 60, -70, 75, 100, 75, -35, -45, 90, 85, -75, 5, -80, 50, -40, -95, 100])
###fMRIlat  = numpy.array([0.25, 0.388, 0.291, 0.1885, -0.1605, 0.4435, 0.6385, 0.0455, 0.3805, -0.358, 0.6595, 0.3045, 0.537, 0.236, 0.202, 0.321, 0.4405, 0.038, 0.435, 0.304, 0.1385, 0.667, 0.337, 0.699, 0.31, 0.3925])
# add new data
handedness = numpy.array([75, 85, 55, -35, -90, -85, -95, 70, -75, -25, 60, -70, 75, 100, 75, -35, -45, 90, 85, -75, 5, -80, 50, -40, -95, 100, -80, -75, -70])
# female=1 male=2
gender = numpy.array([1, 2, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 2,   2, 2, 2])

IDs = numpy.array(['AG_5114', 'AG_5123', 'AG_5125', 'AG_5129', 'AG_5131', 'AG_5132', 'AG_5138', 'AG_5139', 'AG_5142', 'AG_5146', 'AG_5147', 'AG_5148', 'AG_5149', 'AG_5150', 'AG_5151', 'AG_5152', 'AG_5153', 'AG_5154', 'AG_5155', 'AG_5156', 'AG_5158', 'AG_5159', 'AG_5160', 'AG_5162', 'AG_5163', 'AG_5205', 'AG_5263', 'AG_5266', 'AG_5268']) 

exclude = ['5125', '5129']
#exclude = ['5125', '5129', '5139', '5268']

#exclude = []
keep = numpy.empty(len(IDs),dtype=bool)
keep.fill(True)
for exc in exclude:
    for id_idx in range(len(IDs)):
        if exc in IDs[id_idx]:
            keep[id_idx] = False

print "Excluding:"
print IDs[~keep]
print gender[~keep]
print handedness[~keep]
# exclude
IDs = IDs[keep]
gender = gender[keep]
handedness = handedness[keep]

# test that IDs correspond to the input directories
number_IDs = len(IDs)
number_dirs = len(input_dirs)

# test for matches
match_IDs = numpy.empty(number_IDs,dtype=bool)
match_dirs = numpy.empty(number_dirs,dtype=bool)
match_IDs.fill(False)
match_dirs.fill(False)
for id_idx in range(number_IDs):
    for dir_idx in range(number_dirs):
        if IDs[id_idx] in input_dirs[dir_idx]:
            match_IDs[id_idx] = True
            match_dirs[dir_idx] = True

if 0:
    print "<wm_laterality_stats.py> Matched result:"
    print match_IDs
    print "<wm_laterality_stats.py> Matched IDs:"
    print IDs[match_IDs]
    print "<wm_laterality_stats.py> Data are missing for subject IDs:"
    print IDs[~match_IDs]
    print "<wm_laterality_stats.py> Matched directories:"
    print numpy.array(input_dirs)[match_dirs]
    print "<wm_laterality_stats.py> Extra directories that do not correspond to subject IDs:"
    print numpy.array(input_dirs)[~match_dirs]


if number_IDs == number_dirs:
    print "<wm_laterality_stats.py> Number of subject IDs matches the number of directories."
else:
    print "--------------------------------------------------------------------"
    print "<wm_laterality_stats.py> Warning: number of subject IDs:", number_IDs, "does not match the number of directories:", number_dirs
    #print "<wm_laterality_stats.py> Input subject IDs:"
    #print IDs
    print "<wm_laterality_stats.py> Data directories are missing for subject IDs:"
    print IDs[~match_IDs]
    print "<wm_laterality_stats.py> Extra directories that do not correspond to subject IDs:"
    print numpy.array(input_dirs)[~match_dirs]
    print "--------------------------------------------------------------------"
 
# there could be extra directories. Verify the matches are correct.
number_IDs_matched = len(numpy.nonzero(match_IDs)[0])
number_dirs_matched = len(numpy.nonzero(match_dirs)[0])

# grab only the matched ones
IDs = IDs[match_IDs]
input_dirs = numpy.array(input_dirs)[match_dirs]
handedness = handedness[match_IDs]
gender = gender[match_IDs]

print "--------------------------------------------------------------------"
print "<wm_laterality_stats.py> Subject IDs and matching directories for analysis:"
print "--------------------------------------------------------------------"
found_error = False
if number_IDs_matched == number_dirs_matched:
    # test they are in the same order, perfectly matched
    for id, dir in zip(IDs, input_dirs):
        if id in dir:
            print id, ":", dir
        else:
            print "<wm_laterality_stats.py> Error in subject ID and input list.", id, ":", dir
            found_error = True
else:
    print "<wm_laterality_stats.py> Error in matching subjects with input directories. Potential duplicate IDs or directories."
    found_error = True
print "--------------------------------------------------------------------"
if found_error:
    print "<wm_laterality_stats.py> Please check your subject ID list against the subject directories within the input directory."
    raise


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
print "<wm_laterality_stats.py> Reading data..."
for sidx in range(0,len(input_dirs)):
    dir = input_dirs[sidx]
    if verbose:
        print sidx, "/", len(input_dirs), ":", dir
    try:
        resultIO = wma.io.LateralityResults()
        resultIO.read(dir)
        data.append(resultIO)
    except:
        print "<wm_laterality_stats.py>  Unexpected error reading data directory:" , dir
        raise


# change to the output directory to save things there
os.chdir(outdir)

# compute summary stats
print "<wm_laterality_stats.py> Computing summary stats..."
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
    if verbose:
        print sidx, "/", len(data), ":", data[sidx].directory
    LI = data[sidx].laterality_index
    confidence = numpy.minimum(data[sidx].left_hem_similarity, data[sidx].right_hem_similarity)
    confidence_NZ = confidence[numpy.nonzero(LI)]
    # analyze only the LI from hemispheric fibers, as commissural fibers have theirs set to 0.0
    mask_hem = data[sidx].hemisphere != 0
    LI_NZ = LI[mask_hem]
    #LI_NZ = LI[numpy.nonzero(LI)]
    plt.figure()
    plt.plot(confidence_NZ, LI_NZ, 'o')
    plt.savefig('confidence_vs_LI'+str(sidx)+'.pdf')
    plt.close()

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
    #print sidx, ":" , median[sidx]
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

print "<wm_laterality_stats.py> Computing histograms..."
histograms = list()
for sidx in range(0,len(data)) :
    if verbose:
        print sidx, "/", len(data)
    LI = data[sidx].laterality_index
    # analyze only the LI from hemispheric fibers, as commissural fibers have theirs set to 0.0
    mask_hem = data[sidx].hemisphere != 0
    LI_NZ = LI[mask_hem]
    #trimmed_LI_NZ = numpy.sort(LI_NZ)[50:-50]
    #LI_NZ = trimmed_LI_NZ
    #ret = plt.hist(LI_NZ,range=[-1,1],normed=True,bins=20)
    #ret = plt.hist(LI_NZ,range=[-1,1],normed=True,bins=60)
    ret = plt.hist(LI_NZ,range=[-1,1],normed=True,bins=40)
    histograms.append(ret[0])
    #bins = ret[1]
    plt.close()


print "<wm_laterality_stats.py> Making pdf figures..."

# --------------
threshold = 0
#threshold = 60
#threshold = 40
# --------------
gL = numpy.nonzero(handedness <= -threshold)[0]
gR = numpy.nonzero(handedness >= threshold)[0]

# test remove r lang lat by fMRI
#gL=numpy.array([ 3,5,  6,  8, 11, 15, 16, 19, 21, 23, 24])

##### 3 groups (thresh 60 divides into 3 almost even groups 9,9,8)
###thresh = 60
#thresh = 50
if threshold == 0:
    threshold = 60

thresh = threshold
gL3 = numpy.nonzero(handedness < -thresh)[0]
gIC3 = numpy.nonzero(numpy.abs(handedness) <= thresh)[0]
gR3 = numpy.nonzero(handedness > thresh)[0]

# two groups (thresh 75 divides into 2 even groups)
#thresh = 75
gCH2 = numpy.nonzero(numpy.abs(handedness) >= thresh)[0]
gIC2 = numpy.nonzero(numpy.abs(handedness) < thresh)[0]

# any gender differences
gF = numpy.nonzero(gender == 1)[0]
gM = numpy.nonzero(gender == 2)[0]


# plot dataset composition
plt.figure()
plt.hist([handedness[gM], handedness[gF]], label = ['males (n='+str(len(gM))+')', 'females (n='+str(len(gF))+')'])
plt.xlabel('handedness')
plt.ylabel('number of subjects')
plt.title('handedness by gender in current dataset')
plt.legend()
plt.savefig('handedness_by_gender_dataset.pdf')

plt.close()

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
plt.close()

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

for hand_type in ['signed', 'abs']:
    for var in varlist:
        if hand_type == 'signed':
            hand_var = handedness
        else:
            hand_var = numpy.abs(handedness)
        plt.figure()
        #plt.plot(hand_var, var[0],'o')
        plt.plot(hand_var[gM], var[0][gM],'bo')
        plt.plot(hand_var[gF], var[0][gF],'ro')
        (a_s, b_s, r, p, stderr) = scipy.stats.linregress(hand_var, var[0])
        xr=numpy.polyval([a_s,b_s],hand_var)
        plt.plot(hand_var,xr,'k')
        (a_s, b_s, r, p, stderr) = scipy.stats.linregress(hand_var[gM], var[0][gM])
        xr=numpy.polyval([a_s,b_s],hand_var)
        plt.plot(hand_var,xr,'b')
        (a_s, b_s, r, p, stderr) = scipy.stats.linregress(hand_var[gF], var[0][gF])
        xr=numpy.polyval([a_s,b_s],hand_var)
        plt.plot(hand_var,xr,'r')
        plt.title(var[1])
        plt.savefig('handedness_'+hand_type+var[1]+'.pdf')
        plt.close()
        r,p = scipy.stats.pearsonr(hand_var, var[0])
        print var[1],hand_type,".. r:", r, " p: ", p


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
print gM
for sidx in gM:
    print sidx
    plt.plot(histograms_smooth[sidx])

plt.title('M')
plt.savefig('group_M_histograms.pdf')
plt.close()


plt.figure()
for sidx in gF:
    plt.plot(histograms_smooth[sidx])

plt.title('F')
plt.savefig('group_F_histograms.pdf')
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

# median and skew, how related?
plt.figure()
plt.plot(median[gM], skew[gM], 'bo')
plt.plot(median[gF], skew[gF], 'ro')
plt.savefig('groups_median_vs_skew.pdf')

# iqr and kurtosis, how related?
plt.figure()
plt.plot(iqr[gM], kurtosis[gM], 'bo')
plt.plot(iqr[gF], kurtosis[gF], 'ro')
plt.savefig('groups_iqr_vs_kurtosis.pdf')

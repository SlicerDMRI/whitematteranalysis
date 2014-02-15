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

# tensors misoriented
exclude = ['01373']
print "Excluding:", exclude

IDs = numpy.array(sorted(['01190', '01289', '01360', '01222', '01316', '01362', '01229', '01338', '01363', '01246', '01341', '01265', '01342', '01376', '01267', '01343', '01400', '01269', '01352', '01406', '01270', '01357', '01407', '01288', '01358', '01432']))

#IDs = numpy.array(sorted(['01190', '01289', '01360', '01222', '01316', '01362', '01229', '01338', '01363', '01246', '01341', '01373', '01265', '01342', '01376', '01267', '01343', '01400', '01269', '01352', '01406', '01270', '01357', '01407', '01288', '01358', '01432']))

fe_list = ['caseD01222','caseD01229','caseD01246','caseD01267','caseD01270','caseD01289','caseD01343','caseD01352','caseD01360','caseD01362','caseD01373','caseD01405','caseD01406','caseD01450']

status = numpy.empty(len(IDs),dtype=bool)
status.fill(False)
for idx in range(len(IDs)):
    id = IDs[idx]
    for fe in fe_list:
        if id in fe:
            #print "FE:", id
            status[idx] = True
    if not status[idx]:
        #print "CTRL:", id
        pass
    
# add new data
#handedness = numpy.array([75, 85, 55, -35, -90, -85, -95, 70, -75, -25, 60, -70, 75, 100, 75, -35, -45, 90, 85, -75, 5, -80, 50, -40, -95, 100, -80, -75, -70])
# female=1 male=2
#gender = numpy.array([1, 2, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 2,   2, 2, 2])

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
#print match_IDs
input_dirs = numpy.array(input_dirs)[match_dirs]
IDs = IDs[match_IDs]
status = status[match_IDs]

# groups to test
g1 = numpy.nonzero(status)[0]
g2 = numpy.nonzero(~status)[0]

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

# any gender differences
#gF = numpy.nonzero(gender == 1)[0]
#gM = numpy.nonzero(gender == 2)[0]


# significance
# ANOVA and TTEST
p_vals_group = list()
p_vals_MF = list()
values_tested = list()

plt.close('all')
f_list = [median, iqr, skew, kurtosis, num_fibers, num_fibers_Rhem, num_fibers_Lhem, num_fibers_comm, num_fibers_hem]
f_str_list = ['median', 'iqr', 'skew', 'kurtosis', 'num_fibers', 'num_fibers_Rhem', 'num_fibers_Lhem', 'num_fibers_comm', 'num_fibers_hem']
for f, f_str in zip(f_list, f_str_list):
    values_tested.append(f_str)
    [s, p] = scipy.stats.ttest_ind(f[g1],f[g2])
    p_vals_group.append(p)
    plot_box((f[g1],f[g2]), ('1', '2'), 'group_2_'+f_str+'.pdf')
    #[s, p] = scipy.stats.ttest_ind(f[gM],f[gF])
    #p_vals_MF.append(p)
    #plot_box((f[gM],f[gF]), ['M', 'F'], 'group_2gen_median.pdf')

print values_tested
print "FE vs HC"
print p_vals_group
#print "M vs F"
#print p_vals_MF


plt.figure()
hist_sum = numpy.zeros(len(histograms[0]))
for sidx in g1:
    #print sidx
    plt.plot(histograms[sidx])
    plt.savefig('hist_g1_'+str(IDs[sidx])+'.pdf')
    hist_sum = hist_sum + histograms[sidx]
hist_sum = numpy.divide(hist_sum, len(g1))
plt.plot(hist_sum,'b')
hist_sum = numpy.zeros(len(histograms[0]))
for sidx in g2:
    #print sidx
    plt.plot(histograms[sidx])
    plt.savefig('hist_g2_'+str(IDs[sidx])+'.pdf')
    hist_sum = hist_sum + histograms[sidx]
hist_sum = numpy.divide(hist_sum, len(g2))
plt.plot(hist_sum,'r')
plt.title('g1 and g2')
plt.legend(['1', '2'])
plt.savefig('groups_hist.pdf')
plt.close()

histograms_smooth=[]
#kernel = scipy.signal.gaussian(20,5)
kernel = scipy.signal.gaussian(10,1.5)
kernel = kernel/numpy.sum(kernel)

for sidx in range(0,len(data)):
    histograms_smooth.append(scipy.signal.convolve(histograms[sidx], kernel, mode='same'))

plt.figure()
for sidx in g1:
    plt.plot(histograms_smooth[sidx])

plt.title('1')
plt.savefig('group_1_histograms.pdf')
plt.close()

plt.figure()
for sidx in g2:
    plt.plot(histograms_smooth[sidx])

plt.title('2')
plt.savefig('group_2_histograms.pdf')
plt.close()


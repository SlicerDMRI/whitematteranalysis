#!/usr/bin/env python
# If the cluster was created in bilateral clustering and is not commissural, separate into separate directories
# output should be a right_hem and a left_hem directory
# copy MRML files from toplevel directory into this one
# perhaps make a commissural directory also. why not.
# note that it could be cleaner to have a cell data array in the polydata and the option to 
# measure left and right hem parts of a cluster separately.
# Note: if this is performed at the atlas clustering stage, it can be used to separate clusters into groups,
# and this can be learned. At that point all data are midsagitally aligned, which this requires.
# For running this per-subject, the alignment should be performed to handle the tracts near the midline better.
# That should be added as an option.
import numpy
import argparse
import os
import shutil
import vtk
import glob

try:
    import whitematteranalysis as wma
except:
    print "<wm_separate_hemispheres.py> Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Separate each cluster into left/right/commissural tracts according to the percentage of each fiber. "
                "The output is three directories of fiber bundles according to left hemisphere, right hemisphere, and commissural tracts. "
                "Also copies any Slicer MRML scene file that is given as input into the three separate directories.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputDirectory',
    help='A directory of clustered whole-brain tractography as vtkPolyData (.vtk or .vtp).')
parser.add_argument(
    'outputDirectory',
    help='The output directory will be created if it does not exist.')
parser.add_argument(
    '-pthresh', action="store", dest="hemispherePercentThreshold", type=float,
    help='The percent of a fiber that has to be in one hemisphere to consider the fiber as part of that hemisphere (rather than a commissural fiber). '
         'Default number is 0.6, where a higher number will tend to label fewer fibers as hemispheric and more fibers as commissural (not strictly in one hemisphere or the other), '
         'while a lower number will be stricter about what is classified as commissural. '
         'This parameter is not applicable if -clusterLocationFile was provided.')
parser.add_argument(
    '-clusterLocationFile', action="store", dest="clusterLocationFile",
    help='A csv file defining the location of each cluster, i.e., hemispheric or commissural. '
         'The hemispherePercentThreshold will be varied across the clusters if their locations are already known.')
parser.add_argument(
    '-atlasMRML', action="store", dest="atlasMRML", 
    help='A MRML file defining the atlas clusters, to be copied into all directories.')

args = parser.parse_args()


if not os.path.isdir(args.inputDirectory):
    print "<wm_separate_hemispheres.py> Error: Input directory", args.inputDirectory, "does not exist or is not a directory."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "<wm_separate_hemispheres.py> Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

atlasMRML = args.atlasMRML
if not os.path.exists(atlasMRML):
    print "<wm_separate_hemispheres.py> Atlas MRML file", atlasMRML, "does not exist."
    exit()

clusterLocationFile = args.clusterLocationFile
location_data = None
if clusterLocationFile is not None:
    print clusterLocationFile
    if not os.path.exists(clusterLocationFile):
        print "<wm_separate_hemispheres.py> Cluster location file is assigned but the file", clusterLocationFile, "does not exist."
        exit()
    else:
        location_data = numpy.loadtxt(open(clusterLocationFile, "rb"),
                                      dtype={'names': ('Cluster Index', 'Location Label'), 'formats': ('S17', 'S1')},
                                      delimiter="\t", skiprows=1)

# default to be changed if user input is there
hemisphere_percent_threshold = 0.6
if args.hemispherePercentThreshold is not None:
    if (args.hemispherePercentThreshold > 0.5) & (args.hemispherePercentThreshold <= 1.0):
        hemisphere_percent_threshold = args.hemispherePercentThreshold
    else:
        print "<wm_separate_hemispheres.py> Hemisphere fiber percent threshold", args.hemispherePercentThreshold, "must be between 0.5 and 1. (0.6 is recommended)."
        exit()

print "<wm_separate_hemispheres.py> Starting computation."
print ""
print "=====input directory ======\n", args.inputDirectory
print "=====output directory =====\n", args.outputDirectory
print "=========================="
print ""

if location_data is None:
    print "<wm_separate_hemispheres.py> Hemisphere fiber percent threshold", hemisphere_percent_threshold
else:
    print "<wm_separate_hemispheres.py> Separating clusters using location file:", clusterLocationFile
print ""


# relatively high number of points for accuracy
points_per_fiber = 40

def list_cluster_files(input_dir):
    # Find input files
    input_mask = "{0}/cluster_*.vtk".format(input_dir)
    input_mask2 = "{0}/cluster_*.vtp".format(input_dir)
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)

input_polydatas = list_cluster_files(args.inputDirectory)

number_of_subjects = len(input_polydatas)

if location_data is not None and number_of_subjects != len(location_data):
    print "<wm_separate_hemispheres.py> Number of clusters (%d) is not the same as in the location file (%d)." % (number_of_subjects, len(location_data))
    exit()

print "<wm_separate_hemispheres.py> Input number of vtk/vtp files: ", number_of_subjects

outdir_right = os.path.join(outdir, 'tracts_right_hemisphere')
if not os.path.exists(outdir_right):
    os.makedirs(outdir_right)
outdir_left = os.path.join(outdir, 'tracts_left_hemisphere')
if not os.path.exists(outdir_left):
    os.makedirs(outdir_left)
outdir_commissure = os.path.join(outdir, 'tracts_commissural')
if not os.path.exists(outdir_commissure):
    os.makedirs(outdir_commissure)

# copy requested MRML file into these directories
fname_base = os.path.basename(atlasMRML)
fname_output = os.path.join(outdir_right, fname_base)
shutil.copyfile(atlasMRML, fname_output)
fname_output = os.path.join(outdir_left, fname_base)
shutil.copyfile(atlasMRML, fname_output)
fname_output = os.path.join(outdir_commissure, fname_base)
shutil.copyfile(atlasMRML, fname_output)

# midsagittal alignment step to get transform to apply to each polydata
# alternatively could use the transform into atlas space that is found
# during the labeling. This should be an option at the labeling step.
# this transform is the one to use.
#fname1 = os.path.join(outdir, subjectID+'_sym.vtp')

# read in data
hemi_list = []
comm_list = []
for fname, c_idx in zip(input_polydatas, range(len(input_polydatas))):

    # figure out filename and extension
    fname_base = os.path.basename(fname)

    # read data
    print "<wm_separate_hemispheres.py> Reading input file:", fname
    pd = wma.io.read_polydata(fname)

    if pd.GetNumberOfLines() == 0:
        pd_right = pd
        pd_left = pd
        pd_commissure = pd

    else:
        if clusterLocationFile is None:
            hemisphere_percent_threshold = 0.6 # default

            # internal representation for fast similarity computation
            # this also detects which hemisphere fibers are in
            fibers = wma.fibers.FiberArray()
            fibers.points_per_fiber = points_per_fiber
            fibers.hemisphere_percent_threshold = hemisphere_percent_threshold
            # must request hemisphere computation from object
            fibers.hemispheres = True
            # Now convert to array with points and hemispheres as above
            fibers.convert_from_polydata(pd)

            # separate into right and left hemispheres
            # note: this assumes RAS coordinates as far as right and left labels
            # LPS would be switched
            # -------------------------
            mask_right = numpy.zeros(pd.GetNumberOfLines())
            mask_right[fibers.index_right_hem] = 1
            pd_right = wma.filter.mask(pd, mask_right, preserve_point_data=True, preserve_cell_data=True, verbose=False)

            mask_left = numpy.zeros(pd.GetNumberOfLines())
            mask_left[fibers.index_left_hem] = 1
            pd_left = wma.filter.mask(pd, mask_left, preserve_point_data=True, preserve_cell_data=True, verbose=False)

            mask_commissure = numpy.zeros(pd.GetNumberOfLines())
            mask_commissure[fibers.index_commissure] = 1
            pd_commissure = wma.filter.mask(pd, mask_commissure, preserve_point_data=True, preserve_cell_data=True, verbose=False)

        else:
            if location_data[c_idx][1] == 'c' or location_data[c_idx][1] == 'ng':
                comm_list.append(fname_base)

                pd_right = vtk.vtkPolyData()
                pd_left = vtk.vtkPolyData()
                pd_commissure = pd

            elif location_data[c_idx][1] == 'h':
                hemisphere_percent_threshold = 0.5001
                hemi_list.append(fname_base)

                # internal representation for fast similarity computation
                # this also detects which hemisphere fibers are in
                fibers = wma.fibers.FiberArray()
                fibers.points_per_fiber = points_per_fiber
                fibers.hemisphere_percent_threshold = hemisphere_percent_threshold
                # must request hemisphere computation from object
                fibers.hemispheres = True
                # Now convert to array with points and hemispheres as above
                fibers.convert_from_polydata(pd)

                # separate into right and left hemispheres
                # note: this assumes RAS coordinates as far as right and left labels
                # LPS would be switched
                # -------------------------
                mask_right = numpy.zeros(pd.GetNumberOfLines())
                mask_right[fibers.index_right_hem] = 1

                mask_left = numpy.zeros(pd.GetNumberOfLines())
                mask_left[fibers.index_left_hem] = 1
              
                if len(fibers.index_commissure) > 0:
                    if len(fibers.index_left_hem) <= len(fibers.index_right_hem):
                        mask_left[fibers.index_commissure] = 1
                    else:
                        mask_right[fibers.index_commissure] = 1

                pd_right = wma.filter.mask(pd, mask_right, preserve_point_data=True, preserve_cell_data=True, verbose=False)
                pd_left = wma.filter.mask(pd, mask_left, preserve_point_data=True, preserve_cell_data=True, verbose=False)
                pd_commissure = vtk.vtkPolyData()

                #mask_commissure = numpy.zeros(pd.GetNumberOfLines())
                #mask_commissure[fibers.index_commissure] = 1
                #pd_commissure = wma.filter.mask(pd, mask_commissure, preserve_point_data=True, preserve_cell_data=True, verbose=False)

    print ' - fiber number: left %5d - right %5d - comm %5d' % (pd_left.GetNumberOfLines(), pd_right.GetNumberOfLines(), pd_commissure.GetNumberOfLines())

    # output polydatas into the appropriate subdirectories
    fname_output = os.path.join(outdir_right, fname_base)
    wma.io.write_polydata(pd_right, fname_output)
    fname_output = os.path.join(outdir_left, fname_base)
    wma.io.write_polydata(pd_left, fname_output)
    fname_output = os.path.join(outdir_commissure, fname_base)
    wma.io.write_polydata(pd_commissure, fname_output)

def output_mrml(cluster_list, mrml_filename):
    number_of_files = len(cluster_list)
    if number_of_files > 0:
        if number_of_files > 1:
            step = int(100 * 255.0 / (number_of_files - 1))
        elif number_of_files == 1:
            step = int(100 * 255.0)

        R = numpy.array(range(0, 100 * 255 + 1, step)) / 100.0
        G = numpy.abs(range(100 * -127, 100 * 128 + 1, step)) * 2.0 / 100.0
        B = numpy.array(range(100 * 255 + 1, 0, -step)) / 100.0

        colors = list()
        for idx in range(len(cluster_list)):
            colors.append([R[idx], G[idx], B[idx]])

        colors = numpy.array(colors)
        wma.mrml.write(cluster_list, colors, mrml_filename, ratio=1.0)

hemiMRML = os.path.join(outdir, "clustered_tracts_sep_hemispheric_n{0:04d}".format(len(hemi_list))+".mrml")
commMRML = os.path.join(outdir, "clustered_tracts_sep_commissural_n{0:04d}".format(len(comm_list))+".mrml")

print ""
output_mrml(hemi_list, hemiMRML)
output_mrml(comm_list, commMRML)

shutil.copyfile(hemiMRML, os.path.join(outdir_right, os.path.basename(hemiMRML)))
shutil.copyfile(commMRML, os.path.join(outdir_right, os.path.basename(commMRML)))
shutil.copyfile(hemiMRML, os.path.join(outdir_left, os.path.basename(hemiMRML)))
shutil.copyfile(commMRML, os.path.join(outdir_left, os.path.basename(commMRML)))
shutil.copyfile(hemiMRML, os.path.join(outdir_commissure, os.path.basename(hemiMRML)))
shutil.copyfile(commMRML, os.path.join(outdir_commissure, os.path.basename(commMRML)))

print ""
print "<wm_separate_hemispheres.py> Done!!!"

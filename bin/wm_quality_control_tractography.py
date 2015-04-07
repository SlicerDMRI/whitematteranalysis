#!/usr/bin/env python
#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/ipython

# Run registration on the test dataset.

import argparse
import os
import numpy
import vtk
import time

try:
    import whitematteranalysis as wma
except:
    print "<wm_register.py> Error importing white matter analysis package\n"
    raise

HAVE_PLT = 1
try:
    import matplotlib.pyplot as plt
except:
    print "<wm_quality_control.py> Error importing matplotlib.pyplot package, can't plot quality control data.\n"
    HAVE_PLT = 0    

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Perform quality control steps (rendering images and fiber length testing) for all vtk and vtp files in the input directory.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='A directory of whole-brain tractography as vtkPolyData (.vtk or .vtp).')
parser.add_argument(
    'outputDirectory',
    help='Quality control information will be stored in the output directory, which will be created if it does not exist.')
# for now this is not parallelized. that would complicate summary info about group.
#parser.add_argument(
#    '-j', action="store", dest="numberOfJobs", type=int,
#    help='Number of processors to use.')
 
args = parser.parse_args()

print "<quality_control> Starting..."

if not os.path.isdir(args.inputDirectory):
    print "<register> Error: Input directory", args.inputDirectory, "does not exist."
    exit()

output_dir = args.outputDirectory
if not os.path.exists(output_dir):
    print "<register> Output directory", output_dir, "does not exist, creating it."
    os.makedirs(output_dir)

print "<quality_control> Testing all vtk files for quality control (fiber length measurements and rendering to make sure header and gradient orientations are ok)."
input_polydatas = wma.io.list_vtk_files(args.inputDirectory)
number_of_subjects = len(input_polydatas)
print "<quality_control> Found ", number_of_subjects, "subjects in input directory:", args.inputDirectory

# output summary files to save information about what was run
readme_fname = os.path.join(output_dir, 'README.txt')
readme_file = open(readme_fname, 'w')
outstr = "Quality Control Summary\n"
outstr += '----------------------\n'
outstr += '\n'
outstr += "Input Directory: "
outstr += args.inputDirectory
outstr += '\n'
outstr += "Output Directory: "
outstr += args.outputDirectory
outstr += '\n'
outstr += "Number of Subjects: "
outstr += str(number_of_subjects)
outstr += '\n'
outstr += '\n'
outstr +=  "Current date: "  + time.strftime("%x")
outstr += '\n'
outstr +=  "Current time: " + time.strftime("%X")
outstr += '\n'
outstr += '\n'
outstr += "Path to Script: " + os.path.realpath(__file__)
outstr += '\n'
outstr += "Working Directory: " + os.getcwd()
outstr += '\n'
outstr += '\n'
outstr += "Description of Outputs\n"
outstr += '---------------------\n'
outstr += 'fiber_length_histograms.pdf: Distribution of fiber lengths for all subjects.\n'
outstr += 'quality_control_data.txt: Data (e.g. FA, tensors) should match for all subjects.\n'
outstr += 'quality_control_fibers.txt:  Fibers per subject for several length thresholds.\n'
outstr += 'view_*.png: All subjects (colors) should overlap if coordinate system origins are ok.\n'
outstr += 'subject directories: Subject-specific histograms and views for visual inspection.'
outstr += '\n'
outstr += '\n'
outstr += "Command Line Arguments\n"
outstr += '----------------------\n'
outstr += str(args)
outstr += '\n'
outstr += '\n'
outstr += "Input Fiber Files\n"
outstr += '-----------------\n'
for pd in input_polydatas:
    outstr += pd
    outstr += '\n'
readme_file.write(outstr)
readme_file.close()

# output summary files to save information about all subjects
fibers_qc_fname = os.path.join(output_dir, 'quality_control_fibers.txt')
data_qc_fname = os.path.join(output_dir, 'quality_control_data.txt')

fibers_qc_file = open(fibers_qc_fname, 'w')
fiber_test_lengths = [0, 1, 2, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
outstr = "SUBJECT_ID \t"
for test_length in fiber_test_lengths:
    outstr = outstr + "LEN_" + str(test_length) + '\t'
outstr = outstr + '\n'
fibers_qc_file.write(outstr)
fibers_qc_file.close()

data_qc_file = open(data_qc_fname, 'w')
outstr = "SUBJECT_ID \t DATA_INFORMATION (field name, number of components, point or cell data)"
outstr = outstr + '\n'
data_qc_file.write(outstr)
data_qc_file.close()

if HAVE_PLT:
    plt.figure(1)

# Loop over subjects and check each
subject_idx = 1
appender = vtk.vtkAppendPolyData()
for fname in input_polydatas:
    subject_id = os.path.splitext(os.path.basename(fname))[0]
    print "Subject ", subject_idx, "/", number_of_subjects, "ID:", subject_id

    # Read data
    pd = wma.io.read_polydata(fname)

    # Render individual subject
    ren = wma.render.render(pd,1000)
    output_dir_subdir = os.path.join(output_dir, 'tract_QC_' + subject_id)
    if not os.path.exists(output_dir_subdir):
        os.makedirs(output_dir_subdir)
    ren.save_views(output_dir_subdir, subject_id)

    # Compute and save stats about this subject's fiber histogram
    fibers_qc_file = open(fibers_qc_fname, 'a')
    # total number of fibers
    outstr = str(subject_id) +  '\t'
    # numbers of fibers at different possible threshold lengths
    pd2, lengths = wma.filter.preprocess(pd, 100, return_lengths=True)
    lengths = numpy.array(lengths)
    for test_length in fiber_test_lengths:
        number_fibers = numpy.count_nonzero(lengths > test_length)
        outstr = outstr + str(number_fibers) + '\t'
    outstr = outstr + '\n'
    fibers_qc_file.write(outstr)
    fibers_qc_file.close()

    # Save the subject's fiber lengths  
    if HAVE_PLT:
        plt.figure(1)
        plt.hist(lengths, bins=100, histtype='step', label=subject_id)
        plt.figure(2)
        plt.title('Histogram of fiber lengths')
        plt.hist(lengths, bins=100, histtype='step', label=subject_id)
        # Place the legend below the plot so it does not overlap it when there are many subjects
        lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=False, shadow=False, ncol=1)
        # save everything even if the legend is long and goes off the plot
        plt.savefig(os.path.join(output_dir_subdir, 'fiber_length_histogram.pdf'), bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close()

    # Append selected fibers for rendering of all subjects together to check overlap
    number_rendered_fibers = 500
    pd3 = wma.filter.downsample(pd2, number_rendered_fibers)
    mask = numpy.ones(number_rendered_fibers)
    colors = numpy.multiply(mask, subject_idx)
    pd3 = wma.filter.mask(pd3, mask, colors)
    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        appender.AddInputData(pd3)
    else:
        appender.AddInput(pd3)

    # Record what scalar/tensor data is present in this subject's file
    data_qc_file = open(data_qc_fname, 'a')
    outstr = str(subject_id) +  '\t'
    inpointdata = pd.GetPointData()
    incelldata = pd.GetCellData()
    if inpointdata.GetNumberOfArrays() > 0:
        point_data_array_indices = range(inpointdata.GetNumberOfArrays())            
        for idx in point_data_array_indices:
            array = inpointdata.GetArray(idx)
            outstr = outstr + str(array.GetName()) + '\t' + str(array.GetNumberOfComponents()) + '\t' + 'point' + '\t'
    if incelldata.GetNumberOfArrays() > 0:
        cell_data_array_indices = range(incelldata.GetNumberOfArrays())            
        for idx in cell_data_array_indices:
            array = incelldata.GetArray(idx)
            outstr = outstr + str(array.GetName()) + '\t' + str(array.GetNumberOfComponents()) +'\t'  + 'cell' + '\t'
    outstr = outstr + '\n'
    data_qc_file.write(outstr)
    data_qc_file.close()

    # release data
    del pd
    del pd2
    del pd3
    subject_idx += 1

if HAVE_PLT:
    plt.figure(1)
    plt.title('Histogram of fiber lengths')
    # Place the legend below the plot so it does not overlap it when there are many subjects
    lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=False, shadow=False, ncol=1)
    # save everything even if the legend is long and goes off the plot
    plt.savefig((os.path.join(output_dir, 'fiber_length_histograms.pdf')), bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()

print "<quality_control> Final step: rendering all vtk files together to see if any have a different origin (if far apart registration may fail)."
appender.Update()
pd_all = appender.GetOutput()
ren = wma.render.render(pd_all)
ren.save_views(output_dir, "all_subjects")
del ren
del appender

data_qc_file.close()
fibers_qc_file.close()

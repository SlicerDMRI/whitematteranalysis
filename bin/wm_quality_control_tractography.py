#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import time
import warnings

import numpy as np
import vtk

import whitematteranalysis as wma
from whitematteranalysis.utils.opt_pckg import optional_package

matplotlib, have_mpl, _ = optional_package("matplotlib")
plt, _, _ = optional_package("matplotlib.pyplot")

if not have_mpl:
    warnings.warn(matplotlib._msg)
    warnings.warn("Cannot plot quality control data.")


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Perform quality control steps (rendering images and fiber length testing) for all vtk and vtp files in the input directory.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")

    parser.add_argument(
        'inputDirectory',
        help='A directory of whole-brain tractography as vtkPolyData (.vtk or .vtp).')
    parser.add_argument(
        'outputDirectory',
        help='Quality control information will be stored in the output directory, which will be created if it does not exist.')
    # for now this is not parallelized. that would complicate summary info about group.
    # parser.add_argument(
    #    '-j', action="store", dest="numberOfJobs", type=int,
    #    help='Number of processors to use.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)
    
    print(f"<{os.path.basename(__file__)}> Starting...")
    
    if not os.path.isdir(args.inputDirectory):
        print(f"<{os.path.basename(__file__)}> Error: Input directory {args.inputDirectory} does not exist.")
        exit()
    
    output_dir = args.outputDirectory
    if not os.path.exists(output_dir):
        print(f"<{os.path.basename(__file__)}> Output directory {output_dir} does not exist, creating it.")
        os.makedirs(output_dir)
    
    input_polydatas = wma.io.list_vtk_files(args.inputDirectory)
    number_of_subjects = len(input_polydatas)
    print(f"<{os.path.basename(__file__)}> Found {number_of_subjects} subjects in input directory {args.inputDirectory}")
    
    if number_of_subjects < 1:
        print("\n<quality_control> Error: No .vtk or .vtp files were found in the input directory.\n")
        exit()
    
    print(f"<{os.path.basename(__file__)}> Testing all files for quality control (computing fiber length measurements and rendering to make sure header and gradient orientations are ok).")
    print(f"<{os.path.basename(__file__)}> See the README.txt file in the output directory for more overview information.")
    
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
    outstr +=  f"Current date: {time.strftime('%x')}"
    outstr += '\n'
    outstr +=  f"Current time: {time.strftime('%X')}"
    outstr += '\n'
    outstr += '\n'
    outstr += f"Path to Script: {os.path.realpath(__file__)}"
    outstr += '\n'
    outstr += f"Working Directory: {os.getcwd()}"
    outstr += '\n'
    outstr += '\n'
    outstr += "Description of Outputs\n"
    outstr += '---------------------\n'
    outstr += 'fiber_length_histograms.pdf: Distribution of fiber lengths for all subjects.\n'
    outstr += 'quality_control_data.txt: Data (e.g. FA, tensors) should match for all subjects.\n'
    outstr += 'quality_control_fibers.txt:  Step size and numbers of fibers over various lengths (in mm).\n'
    outstr += 'quality_control_spatial_locations.txt:  Bounding box (min/max coordinates of fibers).\n'
    outstr += 'view_*.png: All subjects (colors) should overlap if coordinate system origins are ok.\n'
    outstr += 'subject directories: Subject-specific histograms and views for visual inspection.\n'
    outstr += '\n'
    outstr += 'Open the output text files as tab-delimited files in Excel or other spreadsheet program. \n'
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
    spatial_qc_fname = os.path.join(output_dir, 'quality_control_spatial_locations.txt')
    
    # html files to view results easily
    html_toplevel_fname = os.path.join(output_dir, 'index.html')
    html_left_fname = os.path.join(output_dir, 'view_left.html')
    html_right_fname = os.path.join(output_dir, 'view_right.html')
    html_inf_fname = os.path.join(output_dir, 'view_inf.html')
    html_sup_fname = os.path.join(output_dir, 'view_sup.html')
    html_ant_fname = os.path.join(output_dir, 'view_ant.html')
    html_post_fname = os.path.join(output_dir, 'view_post.html')
    
    html_view_fnames = [html_left_fname, html_right_fname, html_inf_fname, html_sup_fname, html_ant_fname, html_post_fname]
    html_views = ["view_left_", "view_right_", "view_inf_", "view_sup_", "view_ant_", "view_post_"]
    html_views_descrip = ["Left", "Right", "Inferior", "Superior", "Anterior", "Posterior"]
    
    for (fname, descrip) in zip(html_view_fnames, html_views_descrip):
        f = open(fname, 'w')
        outstr = "<!DOCTYPE html>\n<html>\n"
        f.write(outstr)
        outstr = "<style>\n.floated_img\n{\nfloat: left;\n}\n"
        outstr += "body {\nbackground-color: black; color: LightGray;\n}\n"
        outstr += "img {\nborder:1px solid #021a40; margin: 5px; width=\"200\"\n}\n"
        outstr += "p {margin: 5px; width=\"200\"; text-align: center;\n}\n"
        outstr += "h1 {text-align: center;\n}\n"
        outstr += "</style>\n"
        f.write(outstr)
        outstr = "<body>\n<h1>All " + descrip + " Views</h1>\n"
        f.write(outstr)
        f.close()
    
    fibers_qc_file = open(fibers_qc_fname, 'w')
    fiber_test_lengths = [0, 1, 2, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
    outstr = "SUBJECT_ID\tFIBER_STEP_SIZE\tTOTAL_POINTS\tMEAN_FIBER_LENGTH\tTOTAL_FIBERS\t"
    for test_length in fiber_test_lengths[1:]:
        outstr = outstr + "LEN_" + str(test_length) + '\t'
    outstr = f'{outstr}\n'
    fibers_qc_file.write(outstr)
    fibers_qc_file.close()
    
    data_qc_file = open(data_qc_fname, 'w')
    outstr = "SUBJECT_ID\tDATA_INFORMATION (field name, number of components, point or cell data)"
    outstr = f'{outstr}\n'
    data_qc_file.write(outstr)
    data_qc_file.close()
    
    spatial_qc_file = open(spatial_qc_fname, 'w')
    outstr = "SUBJECT_ID\tXmin\tXmax\tYmin\tYmax\tZmin\tZmax"    
    outstr = f'{outstr}\n'
    spatial_qc_file.write(outstr)
    spatial_qc_file.close()
    
    if have_mpl:
        plt.figure(1)
    
    # Loop over subjects and check each
    subject_idx = 1
    appender = vtk.vtkAppendPolyData()
    for fname in input_polydatas:
        subject_id = os.path.splitext(os.path.basename(fname))[0]
        print(f"Subject {subject_idx} / {number_of_subjects} ID: {subject_id}")
    
        # Read data
        pd = wma.io.read_polydata(fname)
    
        # Preprocess for rendering without short fibers and compute fiber lengths
        pd2, lengths, step_size = wma.filter.preprocess(pd, 5, return_lengths=True, verbose=False)
        lengths = np.array(lengths)
        
        # Render individual subject, only including fibers above 5mm.
        ren = wma.render.render(pd2, 1000, verbose=False)
        output_dir_subdir = os.path.join(output_dir, f'tract_QC_{subject_id}')
        if not os.path.exists(output_dir_subdir):
            os.makedirs(output_dir_subdir)
        ren.save_views(output_dir_subdir, subject_id)
        del ren
    
        print('Multiple views for individual subject')
        html_individual_multiviews = os.path.join(output_dir_subdir, f'view_multiple_{subject_id}.html')
        f = open(html_individual_multiviews, 'w')
        outstr = "<!DOCTYPE html>\n<html>\n"
        f.write(outstr)
        outstr = "<style>\n.floated_img\n{\nfloat: left;\n}\n"
        outstr += "body {\nbackground-color: black; color: LightGray;\n}\n"
        outstr += "img {\nborder:1px solid #021a40; margin: 5px; width=\"200\"\n}\n"
        outstr += "p {margin: 5px; width=\"200\"; text-align: center;\n}\n"
        outstr += "h1 {text-align: center;\n}\n"
        outstr += "</style>\n"
        f.write(outstr)
        outstr = f"<body>\n<h1>All {subject_id} Views</h1>\n"
        f.write(outstr)
        f.close()
    
        for (view, descrip) in zip(html_views, html_views_descrip):
            f = open(html_individual_multiviews, 'a')
            img_fname = os.path.join(f'{view}{subject_id}.jpg')
            outstr = "<div class=\"floated_img\">\n"
            outstr+= f"<a href=\"{img_fname}\" ><img src=\"{img_fname}\" alt=\"{subject_id}\" width=\"450\"></a>\n"
            outstr+= f"<p>{descrip}</p>\n</div>\n"
            f.write(outstr)
            f.close()
    
        # Save view information in html file
        for (fname, view) in zip(html_view_fnames, html_views):
            f = open(fname, 'a')
            img_fname = os.path.join(f'tract_QC_{subject_id}', f'{view}{subject_id}.jpg')
            html_fname = os.path.join(f'tract_QC_{subject_id}', f'view_multiple_{subject_id}.html')
            print(output_dir_subdir)
            print(img_fname)
            print(html_individual_multiviews)
            #outstr = "<a href=\"" + img_fname + "\" >" + subject_id + "<img src=\""+ img_fname + "\" alt=\"" + subject_id + "></a>\n"
            #outstr = "<figure>\n"
            #outstr+= "<img src=\"" + img_fname + "\" alt=\"" + subject_id + "\" width=\"300\">\n"
            #outstr+= "<figcaption>" + subject_id + "</figcaption>\n"
            #outstr+= "</figure>\n"
            outstr = "<div class=\"floated_img\">\n"
            outstr+= f"<a href=\"{html_fname}\" ><img src=\"{img_fname}\" alt=\"{subject_id}\"  width=\"300\"></a>\n"
            outstr+= f"<p>{subject_id}</p>\n</div>\n"
            f.write(outstr)
            f.close()
    
        # Compute and save stats about this subject's fiber histogram
        # numbers of fibers at different possible threshold lengths
        pd2, lengths, step_size = wma.filter.preprocess(pd, 20, return_lengths=True, verbose=False)
        lengths = np.array(lengths)
        fibers_qc_file = open(fibers_qc_fname, 'a')
        outstr = f'{str(subject_id)}\t'
        outstr = f'{outstr}{step_size:.4f}\t'
        # total points in the dataset
        outstr = f'{outstr}{str(pd.GetNumberOfPoints())}\t'
        # mean fiber length
        outstr = f'{outstr}{str(np.mean(lengths))}\t'
        # total numbers of fibers
        for test_length in fiber_test_lengths:
            number_fibers = np.count_nonzero(lengths > test_length)
            outstr = f'{outstr}{str(number_fibers)}\t'
        outstr = f'{outstr}\n'
        fibers_qc_file.write(outstr)
        fibers_qc_file.close()
    
        # Save information about the spatial location of the fiber tracts
        spatial_qc_file = open(spatial_qc_fname, 'a')
        outstr = f'{str(subject_id)}\t'
        for bound in pd.GetBounds():
            outstr = f'{outstr}{str(bound)}\t'
        outstr = f'{outstr}\n'
        spatial_qc_file.write(outstr)
        
        # Save the subject's fiber lengths  
        if have_mpl:
            plt.figure(1)
            if lengths.size > 1:
                plt.hist(lengths, bins=100, histtype='step', label=subject_id)
            plt.figure(2)
            plt.title('Histogram of fiber lengths')
            plt.xlabel('fiber length (mm)')
            plt.ylabel('number of fibers')
            if lengths.size > 1:
                plt.hist(lengths, bins=100, histtype='step', label=subject_id)
            else:
                # Plot something so that saving does not fail
                plt.hist([0,0,0,0], bins=1, histtype='step', label=subject_id)
                plt.title('No fibers in dataset.')
            # Place the legend below the plot so it does not overlap it when there are many subjects
            lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=False, shadow=False, ncol=1)
            # save everything even if the legend is long and goes off the plot
            plt.savefig(os.path.join(output_dir_subdir, 'fiber_length_histogram.pdf'), bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.close()
    
        # Append selected fibers for rendering of all subjects together to check overlap
        # number_rendered_fibers = 500
        number_rendered_fibers = 100
        pd3 = wma.filter.downsample(pd2, number_rendered_fibers, verbose=False)
        mask = np.ones(number_rendered_fibers)
        colors = np.multiply(mask, subject_idx)
        pd3 = wma.filter.mask(pd3, mask, colors, verbose=False)
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            appender.AddInputData(pd3)
        else:
            appender.AddInput(pd3)
    
        # Record what scalar/tensor data is present in this subject's file
        data_qc_file = open(data_qc_fname, 'a')
        outstr = f'{str(subject_id)} \t'
        inpointdata = pd.GetPointData()
        incelldata = pd.GetCellData()
        if inpointdata.GetNumberOfArrays() > 0:
            point_data_array_indices = list(range(inpointdata.GetNumberOfArrays()))            
            for idx in point_data_array_indices:
                array = inpointdata.GetArray(idx)
                outstr = f'{outstr}{str(array.GetName())}\t{str(array.GetNumberOfComponents())}\tpoint\t'
        if incelldata.GetNumberOfArrays() > 0:
            cell_data_array_indices = list(range(incelldata.GetNumberOfArrays()))            
            for idx in cell_data_array_indices:
                array = incelldata.GetArray(idx)
                outstr = f'{outstr}{str(array.GetName())}\t{str(array.GetNumberOfComponents())}\tcell\t'
        outstr = f'{outstr}\n'
        data_qc_file.write(outstr)
        data_qc_file.close()
    
        # release data
        del pd
        del pd2
        del pd3
        subject_idx += 1
    
    if have_mpl:
        plt.figure(1)
        plt.title('Histogram of fiber lengths for all subjects')
        plt.xlabel('fiber length (mm)')
        plt.ylabel('number of fibers')
        # Place the legend below the plot so it does not overlap it when there are many subjects
        lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=False, shadow=False, ncol=1)
        # save everything even if the legend is long and goes off the plot
        try:
            plt.savefig((os.path.join(output_dir, 'fiber_length_histograms.pdf')), bbox_extra_artists=(lgd,), bbox_inches='tight')
        except:
            print("Groupwise tract length histogram save failed. Check if the input datasets have any fibers--all may be empty.")
        plt.close()
    
    print(f"<{os.path.basename(__file__)}> Final step: rendering all vtk files together.")
    appender.Update()
    pd_all = appender.GetOutput()
    ren = wma.render.render(pd_all, verbose=False)
    ren.save_views(output_dir, "all_subjects")
    del ren
    del appender
    
    data_qc_file.close()
    fibers_qc_file.close()
    
    
    # Finish html files
    for (fname) in html_view_fnames:
        f = open(fname, 'a')
        img_fname = os.path.join(output_dir_subdir, f'{view}{subject_id}.jpg')
        outstr = "\n</body>\n</html>\n"
        f.write(outstr)
        f.close()

if __name__ == '__main__':
    main()

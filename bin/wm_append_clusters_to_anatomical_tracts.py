#!/usr/bin/env python
import os
import argparse
import vtk
import glob

try:
    import whitematteranalysis as wma
except:
    print("<wm_append_clusters_to_anatomical_tracts.py> Error importing white matter analysis package\n")
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Append multiple fiber clusters into anatomical tracts based on the ORG atlas.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")

parser.add_argument(
    'inputDirectory',
    help='Contains fiber clusters as vtkPolyData file(s).')
parser.add_argument(
    'atlasDirectory',
    help='The ORG atlas folder that contains the anatomical tract MRML files.')
parser.add_argument(
    'outputDirectory',
    help='The output directory should be a new empty directory. It will be created if needed.')

args = parser.parse_args()

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print("<wm_append_clusters_to_anatomical_tracts> Error: Input directory", args.inputDirectory, "does not exist.")
    exit()

inputdir_left = os.path.join(inputdir, 'tracts_left_hemisphere')
inputdir_right = os.path.join(inputdir, 'tracts_right_hemisphere')
inputdir_comm = os.path.join(inputdir, 'tracts_commissural')

if not os.path.isdir(inputdir_left):
    print("<wm_append_clusters_to_anatomical_tracts> Error: Input directory", inputdir_left, "does not exist.")
    exit()
if not os.path.isdir(inputdir_right):
    print("<wm_append_clusters_to_anatomical_tracts> Error: Input directory", inputdir_right, "does not exist.")
    exit()
if not os.path.isdir(inputdir_comm):
    print("<wm_append_clusters_to_anatomical_tracts> Error: Input directory", inputdir_comm, "does not exist.")
    exit()

atlasdir = os.path.abspath(args.atlasDirectory)
if not os.path.isdir(args.atlasDirectory):
    print("<wm_append_clusters_to_anatomical_tracts> Error: Atlas directory", args.atlasDirectory, "does not exist.")
    exit()

def list_mrml_files(input_dir):
    # Find input files
    input_mask = f"{input_dir}/T*.mrml"
    input_mrml_fnames = glob.glob(input_mask)
    return (input_mrml_fnames)

mrml_files = list_mrml_files(atlasdir)

if len(mrml_files) == 0:
    print("<wm_append_clusters_to_anatomical_tracts> Error: There is no mrml files in the input atlas folder.")
else:
    print("<wm_append_clusters_to_anatomical_tracts>", len(mrml_files)-1, "mrml files are detected.")

outdir = os.path.abspath(args.outputDirectory)
if not os.path.exists(args.outputDirectory):
    print("<wm_append_clusters_to_anatomical_tracts> Output directory", args.outputDirectory, "does not exist, creating it.")
    os.makedirs(outdir)

def output_appended_tract(cluster_vtp_list, outputfile):
    appender = vtk.vtkAppendPolyData()
    for c_idx in range(len(cluster_vtp_list)):
        cluster_vtp = cluster_vtp_list[c_idx]
        pd_cluster = wma.io.read_polydata(cluster_vtp)

        vtk_array = vtk.vtkIntArray()
        vtk_array.SetName('cluster_idx')
        for p_idx in range(0, pd_cluster.GetNumberOfPoints()):
            vtk_array.InsertNextTuple1(int(c_idx))

        pd_cluster.GetPointData().AddArray(vtk_array)

        #print '<wm_append_clusters_to_anatomical_tracts>', cluster_vtp, ', number of fibers', pd_cluster.GetNumberOfLines()

        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            appender.AddInputData(pd_cluster)
        else:
            appender.AddInput(pd_cluster)

    appender.Update()
    pd_appended_cluster = appender.GetOutput()

    wma.io.write_polydata(pd_appended_cluster, outputfile)

hemispheric_tracts = ["T_AF", "T_CB", "T_CPC", "T_MdLF", "T_PLIC", "T_SLF-I", "T_SLF-II", "T_SLF-III", "T_EC", "T_EmC", "T_ICP", "T_ILF", "T_IOFF", "T_UF", 
                     "T_Intra-CBLM-I&P", "T_Intra-CBLM-PaT", "T_CR-F", "T_CR-P", "T_CST", "T_SF", "T_SO", "T_SP", "T_TF", "T_TO", "T_TP", 
                     "T_Sup-F", "T_Sup-FP", "T_Sup-O", "T_Sup-OT", "T_Sup-P", "T_Sup-PO", "T_Sup-PT", "T_Sup-T"]

commissural_tracts = ["T_CC1", "T_CC2", "T_CC3", "T_CC4", "T_CC5", "T_CC6", "T_CC7", "T_MCP"]


print("<wm_append_clusters_to_anatomical_tracts> hemispheric tracts (left and right): ")
tract_idx = 1
for tract in hemispheric_tracts:

    print(" *", tract_idx, "-", tract)
    tract_idx = tract_idx + 1
    mrml = os.path.join(atlasdir, tract+".mrml")
    
    if not os.path.exists(mrml):
        print("<wm_append_clusters_to_anatomical_tracts> Error: Cannot locate", mrml)
        exit()

    cluster_vtp_list_left = list()
    cluster_vtp_list_right = list()
    f = open(mrml) 
    for line in f:
        idx = line.find('.vtp')
        if idx > 0:
            cluster_vtp_filename = line[idx-13:idx+4]

            cluster_vtp_filename_left = os.path.join(inputdir_left, cluster_vtp_filename);
            if not os.path.exists(cluster_vtp_filename_left):
                print("<wm_append_clusters_to_anatomical_tracts> Error:", cluster_vtp_filename_left, "does not exist.")
                exit()
            cluster_vtp_list_left.append(cluster_vtp_filename_left)

            cluster_vtp_filename_right = os.path.join(inputdir_right, cluster_vtp_filename);
            if not os.path.exists(cluster_vtp_filename_right):
                print("<wm_append_clusters_to_anatomical_tracts> Error:", cluster_vtp_filename_right, "does not exist.")
                exit()
            cluster_vtp_list_right.append(cluster_vtp_filename_right)

    output_tract_left = os.path.join(outdir, tract + '_left.vtp')
    output_tract_right = os.path.join(outdir, tract + '_right.vtp')

    output_appended_tract(cluster_vtp_list_left, output_tract_left)
    output_appended_tract(cluster_vtp_list_right, output_tract_right)

print("<wm_append_clusters_to_anatomical_tracts> commissural tracts: ")
for tract in commissural_tracts:

    print(" *", tract_idx, "-", tract)
    tract_idx = tract_idx + 1
    
    mrml = os.path.join(atlasdir, tract+".mrml")
    
    if not os.path.exists(mrml):
        print("<wm_append_clusters_to_anatomical_tracts> Error: Cannot locate", mrml)
        exit()

    cluster_vtp_list_comm = list()
    f = open(mrml) 
    for line in f:
        idx = line.find('.vtp')
        if idx > 0:
            cluster_vtp_filename = line[idx-13:idx+4]

            cluster_vtp_filename_comm = os.path.join(inputdir_comm, cluster_vtp_filename);
            if not os.path.exists(cluster_vtp_filename_comm):
                print("<wm_append_clusters_to_anatomical_tracts> Error:", cluster_vtp_filename_comm, "does not exist.")
                exit()
            cluster_vtp_list_comm.append(cluster_vtp_filename_comm)

    output_tract_comm = os.path.join(outdir, tract + '.vtp')

    output_appended_tract(cluster_vtp_list_comm, output_tract_comm)

def list_cluster_files(input_dir):
    # Find input files
    input_mask = f"{input_dir}/T_*.vtk"
    input_mask2 = f"{input_dir}/T_*.vtp"
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)

list_tracts= list_cluster_files(outdir)

print('')
print('<wm_append_clusters_to_anatomical_tracts> Appended tracts can be found at', outdir, '\n')
print('<wm_append_clusters_to_anatomical_tracts> A total of', len(list_tracts), 'tracts\n')


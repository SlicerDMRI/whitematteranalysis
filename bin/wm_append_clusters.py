#!/usr/bin/env python
import os
import argparse
import vtk

try:
    import whitematteranalysis as wma
except:
    print("<wm_append_clusters.py> Error importing white matter analysis package\n")
    raise

def main():
    #-----------------
    # Parse arguments
    #-----------------
    parser = argparse.ArgumentParser(
        description="Append multiple fiber clusters into one fiber tract.",
        epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu",
        version='1.0')
    
    parser.add_argument(
        'inputDirectory',
        help='Contains fiber clusters as vtkPolyData file(s).')
    parser.add_argument(
        'outputDirectory',
        help='The output directory should be a new empty directory. It will be created if needed.')
    parser.add_argument(
        '-appendedTractName', action="store", type=str,
        help="The filename of the output appended fiber tract. No need to give file suffix (.vtp or .vtk). This is an optional parameter; if provided, it will override default constructed filename. ")
    parser.add_argument(
        '-clusterList', action="store", type=int, nargs='+',
        help='A list of clusters to be appended, e.g., 1 2 3. If neither -clusterList nor -tractMRML are provided, all vtk/vtp files in the input folder will be appended.')
    parser.add_argument(
        '-tractMRML', action="store", type=str,
        help='A MRML file that contains the fiber clusters to be appended. If neither -clusterList nor -tractMRML are provided, all vtk/vtp files in the input folder will be appended.')
    
    args = parser.parse_args()
    
    inputdir = os.path.abspath(args.inputDirectory)
    if not os.path.isdir(args.inputDirectory):
        print("<wm_append_clusters> Error: Input directory", args.inputDirectory, "does not exist.")
        exit()
    
    outdir = os.path.abspath(args.outputDirectory)
    if not os.path.exists(args.outputDirectory):
        print("<wm_append_clusters> Error: Output directory", args.outputDirectory, "does not exist, creating it.")
        os.makedirs(outdir)
    
    if (args.clusterList is not None) and (args.tractMRML is not None):
        print("<wm_append_clusters> Error: Only one of -clusterList and -tractMRML can be provided.")
        exit()
    elif (args.clusterList is None) and (args.tractMRML is None):
        print("<wm_append_clusters> All vtk/vtp files in the input directory will be appended.")
    
    outputfilename = 'appended_tract.vtp'
    if (args.appendedTractName is not None):
        outputfilename = args.appendedTractName + '.vtp'
    
    outputfile = os.path.join(outdir, outputfilename)
    
    if args.clusterList is not None:
        cluster_vtp_list = list()
        for c_idx in args.clusterList:
            cluster_vtp_filename = 'cluster_' + str(c_idx).zfill(5) + '.vtp'
            if not os.path.exists(os.path.join(inputdir, cluster_vtp_filename)):
                print("<wm_append_clusters> Error:", cluster_vtp_filename, "does not exist.")
                exit()
            cluster_vtp_list.append(cluster_vtp_filename)
    
    elif args.tractMRML is not None:
        cluster_vtp_list = list()
        f = open(args.tractMRML, 'r') 
        for line in f:
            idx = line.find('.vtp')
            if idx > 0:
                cluster_vtp_filename = line[idx-13:idx+4]
                if not os.path.exists(os.path.join(inputdir, cluster_vtp_filename)):
                    print("<wm_append_clusters> Error:", cluster_vtp_filename, "does not exist.")
                    exit()
                cluster_vtp_list.append(cluster_vtp_filename)
    else:
        cluster_vtp_list = wma.io.list_vtk_files(inputdir)
    
    print("")
    print("<wm_append_clusters> Start to append clusters.")
    print("")
    print("=====input directory======\n", inputdir)
    print("=====output directory=====\n", outdir)
    print("=====clusters to be appended====\n", cluster_vtp_list)
    print("")
    
    appender = vtk.vtkAppendPolyData()
    for c_idx in range(len(cluster_vtp_list)):
        cluster_vtp = cluster_vtp_list[c_idx]
        pd_cluster = wma.io.read_polydata(os.path.join(inputdir, cluster_vtp))
    
        vtk_array = vtk.vtkIntArray()
        vtk_array.SetName('cluster_idx')
        for p_idx in range(0, pd_cluster.GetNumberOfPoints()):
            vtk_array.InsertNextTuple1(int(c_idx))
    
        pd_cluster.GetPointData().AddArray(vtk_array)
    
        print('<wm_append_clusters>', cluster_vtp, ', number of fibers', pd_cluster.GetNumberOfLines())
    
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            appender.AddInputData(pd_cluster)
        else:
            appender.AddInput(pd_cluster)
    
    appender.Update()
    pd_appended_cluster = appender.GetOutput()
    
    wma.io.write_polydata(pd_appended_cluster, outputfile)
    
    print('')
    print('<wm_append_clusters> Appended fiber tract: number of fibers', pd_appended_cluster.GetNumberOfLines())
    print('')
    print('<wm_append_clusters> Save result to', outputfile, '\n')

if __name__ == '__main__':
    main()

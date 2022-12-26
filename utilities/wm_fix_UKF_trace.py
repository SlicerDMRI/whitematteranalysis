#!/usr/bin/env python
import os
import argparse
import vtk
import glob
import numpy

import whitematteranalysis as wma

def list_files(input_dir,regstr):
    # Find input files
    input_mask = ("{0}/"+regstr).format(input_dir)
    input_pd_fnames = glob.glob(input_mask)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)

def fix_trace(inpd):

    inpoints = inpd.GetPoints()
    inpointdata = inpd.GetPointData()

    point_data_array_indices = list(range(inpointdata.GetNumberOfArrays()))
    for idx in point_data_array_indices:
        array = inpointdata.GetArray(idx)

        if array.GetName() == 'trace1' or array.GetName() == 'trace2':
            array_trace_fixed = vtk.vtkFloatArray()
            array_trace_fixed.SetName("correct_"+array.GetName())
            
            inpd.GetLines().InitTraversal()
            for lidx in range(0, inpd.GetNumberOfLines()):
                ptids = vtk.vtkIdList()
                inpd.GetLines().GetNextCell(ptids)
                for pidx in range(0, ptids.GetNumberOfIds()):
                    trace_val = array.GetTuple(ptids.GetId(pidx))[0]
                    # inverse funcion of https://github.com/pnlbwh/ukftractography/blob/fcf83e290de9feb38e6592c2fcacc107cdc029fe/ukf/tractography.cc#L1902
                    trace_correct = 1 / ( numpy.tan(trace_val / 2 * 3.14) ) / 1.0e6 

                    if array.GetName() == 'trace1':
                        array_trace_fixed.InsertNextTuple1(trace_correct)
                    elif array.GetName() == 'trace2':
                        array_trace_fixed.InsertNextTuple1(trace_correct)

            inpd.GetPointData().AddArray(array_trace_fixed)

    inpd.GetPointData().RemoveArray('trace1')
    inpd.GetPointData().RemoveArray('trace2')
    inpd.GetPointData().Update()

    return inpd


def main():
    #-----------------
    # Parse arguments
    #-----------------
    parser = argparse.ArgumentParser(
        description="Fix the incorrect trace value stored by UKF.",
        epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")
    
    parser.add_argument(
        'inputDirectory',
        help='Contains fiber clusters as vtkPolyData file(s).')
    
    args = parser.parse_args()

    fixed_log = os.path.join(args.inputDirectory, 'fixed.log')
    if os.path.exists(fixed_log):
        print(" ** Aleady fixed!")
        return

    vtk_list = list_files(args.inputDirectory, "*.vt*")
    for vtk_file in vtk_list:
        print("Correcting:", vtk_file)
        pd = wma.io.read_polydata(vtk_file)
        pd = fix_trace(pd)
        wma.io.write_polydata(pd, vtk_file)
    
    f = open(fixed_log, 'a')
    f.write('Fixed all')
    f.close()

if __name__ == '__main__':
    main()
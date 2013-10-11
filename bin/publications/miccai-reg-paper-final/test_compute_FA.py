import numpy
import numpy.linalg
import vtk

import whitematteranalysis as wma

epsilon = numpy.finfo(numpy.float).eps

def fractional_anisotropy(w):
    w = numpy.maximum(w,0)
    
    norm = numpy.sqrt(w[0]*w[0] + w[1]*w[1] +  w[2]*w[2]) 
    
    if (norm < epsilon):
        norm = norm + epsilon

    return ((0.70710678)*
            (numpy.sqrt((w[0]-w[1])*(w[0]-w[1]) + 
                        (w[2]-w[1])*(w[2]-w[1]) +
                        (w[2]-w[0])*(w[2]-w[0])))/norm);

def compute_scalar_measures(pd):
    tensors = pd.GetPointData().GetTensors()
    lines = pd.GetLines()

    fa = vtk.vtkFloatArray()
    lambda_parallel = vtk.vtkFloatArray()
    lambda_perp = vtk.vtkFloatArray()
    #fa_array = list()
    
    for tidx in range(tensors.GetNumberOfTuples()):
        if (tidx % 5000) == 0:
            print tidx, '/', tensors.GetNumberOfTuples()
        D = tensors.GetTuple9(tidx)
        D = numpy.array(D).reshape(3,3)
        w, v = numpy.linalg.eig(D)
        fa_tidx = fractional_anisotropy(w)
        lambda_parallel_tidx = w[2]
        lambda_perp_tidx = (w[0] + w[1])/2.0
        fa.InsertNextTuple1(fa_tidx)
        lambda_parallel.InsertNextTuple1(lambda_parallel_tidx)
        lambda_perp.InsertNextTuple1(lambda_perp_tidx)

    fa_avg = vtk.vtkFloatArray()
    fa_lines_list = list()
    fa_avg_list = list()
    lines.InitTraversal()
    for lidx in range(0, pd.GetNumberOfLines()):
        if (lidx % 100) == 0:
            print lidx, '/', pd.GetNumberOfLines()
        pts = vtk.vtkIdList()       
        #lines.GetCell(lidx, pts)
        lines.GetNextCell(pts)
        # compute average FA for this line
        if pts.GetNumberOfIds():
            fa_list = list()
            for pidx in range(0, pts.GetNumberOfIds()):
                fa_list.append(fa.GetTuple1(pts.GetId(pidx)))
            #fa_avg.InsertNextTuple1(numpy.mean(numpy.array(fa_list)))
            fa_avg.InsertNextTuple1(numpy.median(numpy.array(fa_list)))
            fa_lines_list.append(fa_list)
            fa_avg_list.append(numpy.median(numpy.array(fa_list)))
        else:
            fa_avg.InsertNextTuple1(0.0)
            fa_lines_list.append([0.0])
            
    fa_avg.SetName('mean_FA')
    fa.SetName('FA')
    lambda_parallel.SetName('parallel_diffusivity')
    lambda_perp.SetName('perpendicular_diffusivity')

    outpd = pd
    outpd.GetCellData().AddArray(fa_avg)
    outpd.GetCellData().SetActiveScalars('mean_FA')
    outpd.GetPointData().AddArray(fa)
    outpd.GetPointData().AddArray(lambda_parallel)
    outpd.GetPointData().AddArray(lambda_perp)
    outpd.GetPointData().SetActiveScalars('FA')

    return outpd, fa_lines_list, fa_avg_list


def compute_mean_measures(pd):
    lines = pd.GetLines()
    fa = pd.GetPointData().GetArray('FA')
    
    fa_avg = vtk.vtkFloatArray()
    fa_lines_list = list()
    fa_avg_list = list()
    lines.InitTraversal()
    for lidx in range(0, pd.GetNumberOfLines()):
        if (lidx % 100) == 0:
            print lidx, '/', pd.GetNumberOfLines()
        pts = vtk.vtkIdList()       
        lines.GetNextCell(pts)
        # compute average FA for this line
        if pts.GetNumberOfIds():
            fa_list = list()
            for pidx in range(0, pts.GetNumberOfIds()):
                fa_list.append(fa.GetTuple1(pts.GetId(pidx)))
            fa_avg.InsertNextTuple1(numpy.mean(numpy.array(fa_list)))
            #fa_avg.InsertNextTuple1(numpy.median(numpy.array(fa_list)))
        else:
            fa_avg.InsertNextTuple1(0.0)
            
    fa_avg.SetName('mean_FA')

    outpd = pd
    outpd.GetCellData().AddArray(fa_avg)
    outpd.GetCellData().SetActiveScalars('mean_FA')

    return outpd

def compute_mean_array_along_lines(pd, array_name, output_array_name):
    lines = pd.GetLines()
    point_array = pd.GetPointData().GetArray(array_name)
    
    point_array_avg = vtk.vtkFloatArray()
    point_array_lines_list = list()
    point_array_avg_list = list()
    lines.InitTraversal()
    for lidx in range(0, pd.GetNumberOfLines()):
        if (lidx % 100) == 0:
            print lidx, '/', pd.GetNumberOfLines()
        pts = vtk.vtkIdList()       
        lines.GetNextCell(pts)
        # compute average POINT_ARRAY for this line
        if pts.GetNumberOfIds():
            point_array_list = list()
            for pidx in range(0, pts.GetNumberOfIds()):
                point_array_list.append(point_array.GetTuple1(pts.GetId(pidx)))
            point_array_avg.InsertNextTuple1(numpy.mean(numpy.array(point_array_list)))
            #fa_avg.InsertNextTuple1(numpy.median(numpy.array(fa_list)))
        else:
            point_array_avg.InsertNextTuple1(0.0)
            
    point_array_avg.SetName(output_array_name)

    outpd = pd
    outpd.GetCellData().AddArray(point_array_avg)
    outpd.GetCellData().SetActiveScalars(output_array_name)

    return outpd


def compute_max_array_along_lines(pd, array_name, output_array_name):
    lines = pd.GetLines()
    point_array = pd.GetPointData().GetArray(array_name)
    
    point_array_max = vtk.vtkFloatArray()
    point_array_lines_list = list()
    point_array_max_list = list()
    lines.InitTraversal()
    for lidx in range(0, pd.GetNumberOfLines()):
        if (lidx % 100) == 0:
            print lidx, '/', pd.GetNumberOfLines()
        pts = vtk.vtkIdList()       
        lines.GetNextCell(pts)
        # compute mean POINT_ARRAY for this line
        if pts.GetNumberOfIds():
            point_array_list = list()
            for pidx in range(0, pts.GetNumberOfIds()):
                point_array_list.append(point_array.GetTuple1(pts.GetId(pidx)))
            point_array_max.InsertNextTuple1(numpy.max(numpy.array(point_array_list)))
            #fa_max.InsertNextTuple1(numpy.median(numpy.array(fa_list)))
        else:
            point_array_max.InsertNextTuple1(0.0)
            
    point_array_max.SetName(output_array_name)

    outpd = pd
    outpd.GetCellData().AddArray(point_array_max)
    outpd.GetCellData().SetActiveScalars(output_array_name)

    return outpd

if 0:
    pd = wma.io.read_polydata('iteration_4/white_matter_0000.vtk')

    pd, fa_lines_list, fa_avg_list = compute_scalar_measures(pd)

    ren = wma.render.render(pd, scalar_bar=True, data_mode="Point", tube=False, scalar_range=[0,1])

    ren2 = wma.render.render(pd, scalar_bar=True, data_mode="Cell", tube=False, scalar_range=[0,1])



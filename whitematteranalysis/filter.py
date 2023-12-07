# -*- coding: utf-8 -*-

""" filter.py

Swiss army knife of fiber polydata processing.

This module provides filtering functions that operate on vtkPolyData
objects containing whole-brain tractography.

Preprocess by removing fibers based on length, or endpoint
distance. Downsample. Mask.  Create symmetrized test data.

preprocess
downsample
mask
symmetrize
remove_hemisphere
remove_outliers

"""

import os
import warnings

import numpy as np
import vtk

from whitematteranalysis.utils.opt_pckg import optional_package

from . import fibers, similarity

joblib, have_joblib, _ = optional_package("joblib")
Parallel, _, _ = optional_package("joblib.Parallel")
delayed, _, _ = optional_package("joblib.delayed")

if not have_joblib:
    warnings.warn(joblib._msg)
    warnings.warn("Cannot multiprocess.")


verbose = 0


def add_point_data_array(inpd, data, array_name):
    """Make a vtk point data array from input data array and add to inpd.

    Input data must have dimensions of the number of lines in the
    input polydata. The output array will be added to the polydata and
    will be point data, so the per-line values will be duplicated to
    become per-point values (each point on the line) for visualization.
    """

    ptids = vtk.vtkIdList()
    inpd.GetLines().InitTraversal()
    outarray = vtk.vtkFloatArray()
    outarray.SetName(array_name)
    for lidx in range(0, inpd.GetNumberOfLines()):
        inpd.GetLines().GetNextCell(ptids)
        for pidx in range(0, ptids.GetNumberOfIds()):
            outarray.InsertNextTuple1(data[lidx])
    inpd.GetPointData().AddArray(outarray)

def flatten_length_distribution(inpd, min_length_mm=None, max_length_mm=None, num_bins=10, fibers_per_bin=1000, verbose=True):
    """In order to cluster all structures without removing more
    prevalent shorter ones, sample equal numbers across length ranges.

    This should enable clustering of smaller u-fibers than previously
    possible, while incorporating information about longer fibers they are
    near for stability.
    """

    # try to get N fibers in each length range ("bin")
    # First calculate the ranges using input or measured max/min lengths
    if (max_length_mm is None) or (min_length_mm is None):
        # call preprocess just to get the lengths measured
        inpd2, lengths, step_size = preprocess(inpd, 0.0, max_length_mm=max_length_mm, verbose=False, return_lengths=True)
        if max_length_mm is None:
            max_length_mm = np.max(lengths)
        if min_length_mm is None:
            min_length_mm = np.min(lengths)

    increment = (max_length_mm - min_length_mm) / (num_bins - 1)
    bin_ends = list()
    max_l = min_length_mm
    while max_l <= max_length_mm:
        bin_ends.append(max_l)
        max_l += increment
    if verbose:
        print(f"Bins/length ranges: {bin_ends}")

    print(bin_ends[0:-1], bin_ends[1:])

    # append the sampled fibers together into a new polydata
    appender = vtk.vtkAppendPolyData()
    for (bin_low, bin_hi) in zip (bin_ends[0:-1], bin_ends[1:]):
        pd = preprocess(inpd, bin_low, max_length_mm=bin_hi, verbose=False)
        pd2 = downsample(pd, fibers_per_bin,verbose=False)
        if verbose:
            print(f"{pd2.GetNumberOfLines()} fibers in length range [{bin_low}, {bin_hi}]")
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            appender.AddInputData(pd2)
        else:
            appender.AddInput(pd2)
    appender.Update()
    return(appender.GetOutput())

def compute_lengths(inpd):
    """Compute length of each fiber in polydata. Returns lengths and step size.

    Step size is estimated using points in the middle of a fiber with over 15 points.

    """

    # Make sure we have lines and points.
    if (inpd.GetNumberOfLines() == 0) or (inpd.GetNumberOfPoints() == 0):
        print(f"<{os.path.basename(__file__)}> No fibers found in input polydata.")
        return 0, 0
    
    # measure step size (using first line that has >=5 points)
    cell_idx = 0
    ptids = vtk.vtkIdList()
    inpoints = inpd.GetPoints()
    inpd.GetLines().InitTraversal()
    while (ptids.GetNumberOfIds() < 5) & (cell_idx < inpd.GetNumberOfLines()):
        inpd.GetLines().GetNextCell(ptids)
        ##    inpd.GetLines().GetCell(cell_idx, ptids)
        ## the GetCell function is not wrapped in Canopy python-vtk
        cell_idx += 1
    # make sure we have some points along this fiber
    # In case all fibers in the brain are really short, treat it the same as no fibers.
    if ptids.GetNumberOfIds() < 5:
        return 0, 0
    
    # Use points from the middle of the fiber to estimate step length.
    # This is because the step size may vary near endpoints (in order to include
    # endpoints when downsampling the fiber to reduce file size).
    step_size = 0.0
    count = 0.0
    for ptidx in range(1, ptids.GetNumberOfIds()-1):
        point0 = inpoints.GetPoint(ptids.GetId(ptidx))
        point1 = inpoints.GetPoint(ptids.GetId(ptidx + 1))
        step_size += np.sqrt(np.sum(np.power(np.subtract(point0, point1), 2)))
        count += 1
    step_size = step_size / count

    fiber_lengths = list()
    # loop over lines
    inpd.GetLines().InitTraversal()
    num_lines = inpd.GetNumberOfLines()
    for lidx in range(0, num_lines):
        inpd.GetLines().GetNextCell(ptids)
        # save length
        fiber_lengths.append(ptids.GetNumberOfIds() * step_size)

    return np.array(fiber_lengths), step_size

def preprocess(inpd, min_length_mm,
               remove_u=False,
               remove_u_endpoint_dist=40,
               remove_brainstem=False,
               return_indices=False,
               return_lengths=False,
               preserve_point_data=False,
               preserve_cell_data=False,
               verbose=True, max_length_mm=None):
    """Remove fibers below a length threshold and using other criteria (optional).

    Based on fiber length, and optionally on distance between
    endpoints (u-shape has low distance), and inferior location
    (likely in brainstem).

    """

    # Make sure we have lines and points.
    if (inpd.GetNumberOfLines() == 0) or (inpd.GetNumberOfPoints() == 0):
        print(f"<{os.path.basename(__file__)}> No fibers found in input polydata.")
        if return_indices:
            if return_lengths:
                return inpd, 0, 0, 0
            else:
                return inpd, 0
        else:
            if return_lengths:
                return inpd, 0, 0
            else:
                return inpd

    
    fiber_lengths, step_size = compute_lengths(inpd)

    #print "LINES:", inpd.GetNumberOfLines(), "STEP", step_size

    min_length_pts = round(min_length_mm / float(step_size))
    if verbose:
        print(f"<{os.path.basename(__file__)}> Minimum length {min_length_mm} mm. Tractography step size * minimum number of points = {step_size} * {min_length_pts}")

    # set up processing and output objects
    ptids = vtk.vtkIdList()
    inpoints = inpd.GetPoints()

    # loop over lines
    inpd.GetLines().InitTraversal()
    num_lines = inpd.GetNumberOfLines()
    # keep track of the lines we will keep
    line_indices = list()
 
    for lidx in range(0, num_lines):
        inpd.GetLines().GetNextCell(ptids)

        # first figure out whether to keep this line
        keep_curr_fiber = False

        # test for line being long enough
        if ptids.GetNumberOfIds() > min_length_pts:
            keep_curr_fiber = True

            if max_length_mm is not None:
                if ptids.GetNumberOfIds() * step_size > max_length_mm:
                    keep_curr_fiber = False
                
            if remove_u | remove_brainstem:
                # find first and last points on the fiber
                ptid = ptids.GetId(0)
                point0 = inpoints.GetPoint(ptid)
                ptid = ptids.GetId(ptids.GetNumberOfIds() - 1)
                point1 = inpoints.GetPoint(ptid)

            if remove_u:
                # compute distance between endpoints
                endpoint_dist = np.sqrt(np.sum(np.power(
                            np.subtract(point0, point1), 2)))
                if endpoint_dist < remove_u_endpoint_dist:
                    keep_curr_fiber = False

            if remove_brainstem:
                # compute average SI (third coordinate) < -40
                mean_sup_inf = (point0[2] + point1[2]) / 2
                if mean_sup_inf < -40:
                    keep_curr_fiber = False

        # if we are keeping the line record its index
        if keep_curr_fiber:
            line_indices.append(lidx)

    fiber_mask = np.zeros(num_lines)
    fiber_mask[line_indices] = 1
    outpd = mask(inpd, fiber_mask, preserve_point_data=preserve_point_data, preserve_cell_data=preserve_cell_data, verbose=verbose)

    if return_indices:
        if return_lengths:
            return outpd, np.array(line_indices), fiber_lengths, step_size
        else:
            return outpd, np.array(line_indices)
    else:
        if return_lengths:
            return outpd, fiber_lengths, step_size
        else:
            return outpd

def downsample(inpd, output_number_of_lines, return_indices=False, preserve_point_data=False, preserve_cell_data=True, initial_indices=None, verbose=True, random_seed=1234):
    """ Random (down)sampling of fibers without replacement. """

    if initial_indices is None:
        num_lines = inpd.GetNumberOfLines()
    else:
        num_lines = len(initial_indices)

    if num_lines < output_number_of_lines:
        return inpd

    # use the input random seed every time for code testing experiments
    if random_seed is not None:
        if verbose:
            print(f"<{os.path.basename(__file__)}> Setting random seed to {random_seed}")
        np.random.seed(seed=random_seed)

    # randomly pick the lines that we will keep
    line_indices = np.random.permutation(num_lines)
    if initial_indices is None:
        line_indices = line_indices[0:output_number_of_lines]
        fiber_mask = np.zeros(num_lines)
        fiber_mask[line_indices] = 1
    else:
        line_indices = initial_indices[line_indices[0:output_number_of_lines]]
        fiber_mask = np.zeros(inpd.GetNumberOfLines())
        fiber_mask[line_indices] = 1

    # don't color by line index by default, preserve whatever was there
    #outpd = mask(inpd, fiber_mask, fiber_mask)
    outpd = mask(inpd, fiber_mask, preserve_point_data=preserve_point_data, preserve_cell_data=preserve_cell_data, verbose=verbose)

    # final line count
    #print(f"<{os.path.basename(__file__)}> Number of lines selected: {outpd.GetNumberOfLines()}")
    if return_indices:
        # return sorted indices, this is the line ordering of output
        # polydata (because we mask rather than changing input line order)
        return outpd, np.sort(line_indices)
    else:
        return outpd


def mask(inpd, fiber_mask, color=None, preserve_point_data=False, preserve_cell_data=True, verbose=True):
    """ Keep lines and their points where fiber_mask == 1.

     Unlike vtkMaskPolyData that samples every nth cell, this function
     uses an actual mask, and also gets rid of points that
     are not used by any cell, reducing the size of the polydata file.

     This code also sets scalar cell data to input color data.  This
     input data is expected to be 1 or 3 components.

     If there is no input cell scalar color data, existing cell
     scalars that we have created (EmbeddingColor, ClusterNumber,
     EmbeddingCoordinate) are looked for and masked as well.

     """

    inpoints = inpd.GetPoints()
    inpointdata = inpd.GetPointData()
    incelldata = inpd.GetCellData()
    
    # output and temporary objects
    ptids = vtk.vtkIdList()
    outpd = vtk.vtkPolyData()
    outlines = vtk.vtkCellArray()
    outpoints = vtk.vtkPoints()
    outcolors = None
    outpointdata = outpd.GetPointData()
    outcelldata = outpd.GetCellData()
    tensor_names = []

    if color is not None:
        # if input is RGB
        if len(color.shape) == 2:
            if color.shape[1] == 3:
                outcolors = vtk.vtkUnsignedCharArray()
                outcolors.SetNumberOfComponents(3)

        # otherwise output floats as colors
        if outcolors == None:
            outcolors = vtk.vtkFloatArray()

    # check for cell data arrays to keep
    if preserve_cell_data:
        if incelldata.GetNumberOfArrays() > 0:
            cell_data_array_indices = list(range(incelldata.GetNumberOfArrays()))            
            for idx in cell_data_array_indices:
                array = incelldata.GetArray(idx)
                dtype = array.GetDataType()
                if dtype == 10:
                    out_array = vtk.vtkFloatArray()
                elif dtype == 6:
                    out_array = vtk.vtkIntArray()
                elif dtype == 3:
                    out_array = vtk.vtkUnsignedCharArray()
                else:
                    out_array = vtk.vtkFloatArray()
                out_array.SetNumberOfComponents(array.GetNumberOfComponents())
                out_array.SetName(array.GetName())
                if verbose:
                    print(f"Cell data array found: {array.GetName()} {array.GetNumberOfComponents()}")
                outcelldata.AddArray(out_array)
                # make sure some scalars are active so rendering works
                #outpd.GetCellData().SetActiveScalars(array.GetName())
                
            #if inpd.GetCellData().GetArray('ClusterNumber'):
            #    # this will be active unless we have embedding colors
            #    outpd.GetCellData().SetActiveScalars('ClusterNumber')
            #if inpd.GetCellData().GetArray('EmbeddingColor'):
            #    # Note Slicer can't display these cell scalars (all is red)
            #    outpd.GetCellData().SetActiveScalars('EmbeddingColor')

        else:
            preserve_cell_data = False

    # check for point data arrays to keep
    if preserve_point_data:
        if inpointdata.GetNumberOfArrays() > 0:
            point_data_array_indices = list(range(inpointdata.GetNumberOfArrays()))            
            for idx in point_data_array_indices:
                array = inpointdata.GetArray(idx)
                out_array = vtk.vtkFloatArray()
                out_array.SetNumberOfComponents(array.GetNumberOfComponents())
                out_array.SetName(array.GetName())
                if verbose:
                    print(f"Point data array found: {array.GetName()} {array.GetNumberOfComponents()}")
                outpointdata.AddArray(out_array)
                # make sure some scalars are active so rendering works
                #outpd.GetPointData().SetActiveScalars(array.GetName())
                # keep track of tensors to choose which is active
                if array.GetNumberOfComponents() == 9:
                    tensor_names.append(array.GetName())
        else:
            preserve_point_data = False

    # Set up scalars and tensors attributes for correct visualization in Slicer.
    # Slicer works with point data and does not handle cell data well.
    # This set of changes breaks old internal wma default visualization of cell scalars. 
    # Changes must be propagated through wma so that render is called with the name of the field to visualize.
    # the new way in wma is like this line, below.
    #ren = wma.render.render(output_polydata_s, 1000, data_mode="Cell", data_name='EmbeddingColor')

    # For Slicer: First set one of the expected tensor arrays as default for vis
    tensors_labeled = False
    for name in tensor_names:
        if name == "tensors":
            outpd.GetPointData().SetTensors(outpd.GetPointData().GetArray("tensors"))
            tensors_labeled = True
        if name == "Tensors":
            outpd.GetPointData().SetTensors(outpd.GetPointData().GetArray("Tensors"))
            tensors_labeled = True
        if name == "tensor1":
            outpd.GetPointData().SetTensors(outpd.GetPointData().GetArray("tensor1"))
            tensors_labeled = True
        if name == "Tensor1":
            outpd.GetPointData().SetTensors(outpd.GetPointData().GetArray("Tensor1"))
            tensors_labeled = True
    if not tensors_labeled:
        if len(tensor_names) > 0:
            print(f"Data has unexpected tensor name(s). Unable to set active for visualization: {tensor_names}")
    # now set cell data visualization inactive.
    outpd.GetCellData().SetActiveScalars(None)
                
    # loop over lines
    inpd.GetLines().InitTraversal()
    outlines.InitTraversal()

    for lidx in range(0, inpd.GetNumberOfLines()):
        inpd.GetLines().GetNextCell(ptids)

        if fiber_mask[lidx]:

            if verbose:
                if lidx % 100 == 0:
                    print(f"<{os.path.basename(__file__)}> Line: {lidx} / {inpd.GetNumberOfLines()}")

            # get points for each ptid and add to output polydata
            cellptids = vtk.vtkIdList()

            for pidx in range(0, ptids.GetNumberOfIds()):
                point = inpoints.GetPoint(ptids.GetId(pidx))
                idx = outpoints.InsertNextPoint(point)
                cellptids.InsertNextId(idx)
                if preserve_point_data:
                    for idx in point_data_array_indices:
                        array = inpointdata.GetArray(idx)
                        out_array = outpointdata.GetArray(idx)
                        out_array.InsertNextTuple(array.GetTuple(ptids.GetId(pidx)))

            outlines.InsertNextCell(cellptids)                    

            if color is not None:
                # this code works with either 3 or 1 component only
                if outcolors.GetNumberOfComponents() == 3:
                    outcolors.InsertNextTuple3(color[lidx,0], color[lidx,1], color[lidx,2])
                else:
                    outcolors.InsertNextTuple1(color[lidx])

            if preserve_cell_data:
                for idx in cell_data_array_indices:
                    array = incelldata.GetArray(idx)
                    out_array = outcelldata.GetArray(idx)
                    out_array.InsertNextTuple(array.GetTuple(lidx))

    # put data into output polydata
    outpd.SetLines(outlines)
    outpd.SetPoints(outpoints)

    # if we had an input color requested during masking, set that to be the default scalar for vis
    if color is not None:
        outpd.GetCellData().SetScalars(outcolors)

    if verbose:
        print(f"<{os.path.basename(__file__)}> Fibers sampled: {outpd.GetNumberOfLines()} / {inpd.GetNumberOfLines()}")

    return outpd


def symmetrize(inpd):
    """Generate symmetric polydata by reflecting.

    Output polydata has twice as many lines as input.

    """

    # output and temporary objects
    ptids = vtk.vtkIdList()
    points = inpd.GetPoints()
    outpd = vtk.vtkPolyData()
    outlines = vtk.vtkCellArray()
    outpoints = vtk.vtkPoints()
    outpoints.DeepCopy(points)

    # loop over lines
    inpd.GetLines().InitTraversal()
    outlines.InitTraversal()

    # set scalar cell data to 1 for orig, -1 for reflect, for vis
    outcolors = vtk.vtkFloatArray()

    # index into end of point array
    lastidx = outpoints.GetNumberOfPoints()
    print(f"<{os.path.basename(__file__)}> Input number of points: {lastidx}")

    # loop over all lines, insert line and reflected copy into output pd
    for lidx in range(0, inpd.GetNumberOfLines()):
        # progress
        if verbose:
            if lidx % 100 == 0:
                print(f"<{os.path.basename(__file__)}> Line: {lidx} / {inpd.GetNumberOfLines()}")

        inpd.GetLines().GetNextCell(ptids)

        num_points = ptids.GetNumberOfIds()

        # insert fiber (ptids are same since new points go at the end)
        outlines.InsertNextCell(ptids)
        outcolors.InsertNextTuple1(1)

        # insert reflection into END of point array and into line array
        refptids = vtk.vtkIdList()
        for pidx in range(0, num_points):

            point = points.GetPoint(ptids.GetId(pidx))

            # reflect (RAS -> reflect first value)
            refpoint = (-point[0], point[1], point[2])
            idx = outpoints.InsertNextPoint(refpoint)
            refptids.InsertNextId(idx)

        outlines.InsertNextCell(refptids)
        outcolors.InsertNextTuple1(-1)

    # put data into output polydata
    outpd.SetLines(outlines)
    outpd.SetPoints(outpoints)
    outpd.GetCellData().SetScalars(outcolors)

    return outpd


def remove_hemisphere(inpd, hemisphere=-1):
    """ Remove left (-1) or right (+1) hemisphere points. """

    # output and temporary objects
    ptids = vtk.vtkIdList()
    outpd = vtk.vtkPolyData()
    outlines = vtk.vtkCellArray()
    outlines.InitTraversal()
    outpoints = vtk.vtkPoints()

    # loop over lines
    inpd.GetLines().InitTraversal()

    # loop over all lines, inserting into output only
    # the part in the correct hemisphere
    for lidx in range(0, inpd.GetNumberOfLines()):
        # progress
        if verbose:
            if lidx % 100 == 0:
                print(f"<{os.path.basename(__file__)}> Line: {lidx} / {inpd.GetNumberOfLines()}")

        inpd.GetLines().GetNextCell(ptids)

        num_points = ptids.GetNumberOfIds()

        # insert kept points into point array and into line array
        keptptids = vtk.vtkIdList()
        for pidx in range(0, num_points):

            point = inpd.GetPoints().GetPoint(ptids.GetId(pidx))

            # if we keep this point (if not in removed hemisphere)
            if ((point[0] < 0) &
                (hemisphere == 1)) | ((point[0] > 0) &
                                      (hemisphere == -1)):

                idx = outpoints.InsertNextPoint(point)
                keptptids.InsertNextIdx(idx)

        outlines.InsertNextCell(keptptids)

    # put data into output polydata
    outpd.SetLines(outlines)
    outpd.SetPoints(outpoints)

    return outpd


def remove_outliers(inpd, min_fiber_distance, n_jobs=0, distance_method ='Mean'):
    """ Remove fibers that have no other nearby fibers, i.e. outliers.

    The pairwise fiber distance matrix is computed, then fibers
    are rejected if their average neighbor distance (using closest 3
    neighbors) is higher than min_fiber_distance.

    """

    fiber_array = fibers.FiberArray()
    #fiber_array.points_per_fiber = 5
    fiber_array.points_per_fiber = 10
    fiber_array.convert_from_polydata(inpd)

    fiber_indices = list(range(0, fiber_array.number_of_fibers))

    # squared distances are computed
    min_fiber_distance = min_fiber_distance * min_fiber_distance
    
    # pairwise distance matrix
    if have_joblib and n_jobs > 0:
        distances = Parallel(n_jobs=n_jobs, verbose=1)(
            delayed(similarity.fiber_distance)(
                fiber_array.get_fiber(lidx),
                fiber_array,
                threshold = 0,
                distance_method = distance_method)
            for lidx in fiber_indices)

        distances = np.array(distances)

        # now we check where there are no nearby fibers in d
        mindist = np.zeros(fiber_array.number_of_fibers)
        for lidx in fiber_indices:
            dist = np.sort(distances[lidx, :])
            # robust minimum distance
            mindist[lidx] = (dist[1] + dist[2] + dist[3]) / 3.0
            #mindist[lidx] = (dist[1] + dist[2]) / 2.0

    else:
        # do this in a loop to use less memory. then parallelization can 
        # happen over the number of subjects.
        mindist = np.zeros(fiber_array.number_of_fibers)
        for lidx in fiber_indices:
            distances = similarity.fiber_distance(fiber_array.get_fiber(lidx), fiber_array, 0,  distance_method = distance_method)
            dist = np.sort(distances)
            # robust minimum distance
            mindist[lidx] = (dist[1] + dist[2] + dist[3]) / 3.0
            
    # keep only fibers who have nearby similar fibers
    fiber_mask = mindist < min_fiber_distance

    if True:
        num_fibers = f"{len(np.nonzero(fiber_mask)[0])} / {len(fiber_mask)}"
        print(f"<{os.path.basename(__file__)}> Number retained after outlier removal: {num_fibers}")

    outpd = mask(inpd, fiber_mask, mindist)
    outpd_reject = mask(inpd, ~fiber_mask, mindist)

    return outpd, fiber_mask, outpd_reject

def smooth(inpd, fiber_distance_sigma = 25, points_per_fiber=30, n_jobs=2, upper_thresh=30):
    """ Average nearby fibers.
    
    The pairwise fiber distance matrix is computed, then fibers
    are averaged with their neighbors using Gaussian weighting.

    The "local density" or soft neighbor count is also output.
    """

    sigmasq = fiber_distance_sigma * fiber_distance_sigma
    
    # polydata to array conversion, fixed-length fiber representation
    current_fiber_array = fibers.FiberArray()
    current_fiber_array.points_per_fiber = points_per_fiber
    current_fiber_array.convert_from_polydata(inpd)

    # fiber list data structure initialization for easy fiber averaging
    curr_count = list()
    curr_fibers = list()
    next_fibers = list()
    next_weights = list()
    for lidx in range(0, current_fiber_array.number_of_fibers):
        curr_fibers.append(current_fiber_array.get_fiber(lidx))
        curr_count.append(1)

    fiber_indices = list(range(0, current_fiber_array.number_of_fibers))

    # compare squared distances to squared distance threshold
    upper_thresh = upper_thresh*upper_thresh
    
    print(f"<{os.path.basename(__file__)}> Computing pairwise distances...")
    
    # pairwise distance matrix
    if have_joblib:
        distances = Parallel(n_jobs=n_jobs, verbose=1)(
            delayed(similarity.fiber_distance)(
            current_fiber_array.get_fiber(lidx),
            current_fiber_array,
            0, 'Hausdorff')
            for lidx in fiber_indices)
        distances = np.array(distances)
    else:
        distances = \
            np.zeros(
            (current_fiber_array.number_of_fibers,
             current_fiber_array.number_of_fibers))
        for lidx in fiber_indices:
            distances[lidx, :] = \
                similarity.fiber_distance(
                    current_fiber_array.get_fiber(lidx),
                    current_fiber_array, 0)
   
    # gaussian smooth all fibers using local neighborhood
    for fidx in fiber_indices:
        if (fidx % 100) == 0:
            print(f'{fidx} / {current_fiber_array.number_of_fibers}')

        # find indices of all nearby fibers
        indices = np.nonzero(distances[fidx] < upper_thresh)[0]
        local_fibers = list()
        local_weights = list()

        for idx in indices:
            dist = distances[fidx][idx]
            # these are now squared distances
            weight = np.exp(-dist/sigmasq)
            #weight = np.exp(-(dist*dist)/sigmasq)
            local_fibers.append(curr_fibers[idx] * weight)
            local_weights.append(weight)
        # actually perform the weighted average
        # start with the one under the center of the kernel
        #out_fiber = curr_fibers[fidx]
        #out_weights = 1.0
        out_fiber = local_fibers[0]
        out_weights = local_weights[0]
        for fiber in local_fibers[1:]:
            out_fiber += fiber
        for weight in local_weights[1:]:
            out_weights += weight
        out_fiber = out_fiber / out_weights
        next_fibers.append(out_fiber)
        next_weights.append(out_weights)

    # set up array for output
    output_fiber_array = fibers.FiberArray()    
    output_fiber_array.number_of_fibers = len(curr_fibers)
    output_fiber_array.points_per_fiber = points_per_fiber
    dims = [output_fiber_array.number_of_fibers, output_fiber_array.points_per_fiber]
    # fiber data
    output_fiber_array.fiber_array_r = np.zeros(dims)
    output_fiber_array.fiber_array_a = np.zeros(dims)
    output_fiber_array.fiber_array_s = np.zeros(dims)
    next_fidx = 0
    for next_fib in next_fibers:
        output_fiber_array.fiber_array_r[next_fidx] = next_fib.r
        output_fiber_array.fiber_array_a[next_fidx] = next_fib.a
        output_fiber_array.fiber_array_s[next_fidx] = next_fib.s
        next_fidx += 1
            
    # convert output to polydata
    outpd = output_fiber_array.convert_to_polydata()
    
    # color by the weights or "local density"
    # color output by the number of fibers that each output fiber corresponds to
    outcolors = vtk.vtkFloatArray()
    outcolors.SetName('KernelDensity')
    for weight in next_weights:
        outcolors.InsertNextTuple1(weight)
    #outpd.GetCellData().SetScalars(outcolors)
    outpd.GetCellData().AddArray(outcolors)
    outpd.GetCellData().SetActiveScalars('KernelDensity')

    return outpd, np.array(next_weights)
    
def anisotropic_smooth(inpd, fiber_distance_threshold, points_per_fiber=30, n_jobs=2, cluster_max = 10):
    """ Average nearby fibers.
    
    The pairwise fiber distance matrix is computed, then fibers
    are averaged with their neighbors until an edge (>max_fiber_distance) is encountered.

    """

    # polydata to array conversion, fixed-length fiber representation
    current_fiber_array = fibers.FiberArray()
    current_fiber_array.points_per_fiber = points_per_fiber
    current_fiber_array.convert_from_polydata(inpd)
    original_number_of_fibers = current_fiber_array.number_of_fibers
    
    # fiber list data structure initialization for easy fiber averaging
    curr_count = list()
    curr_fibers = list()
    curr_indices = list()
    for lidx in range(0, current_fiber_array.number_of_fibers):
        curr_fibers.append(current_fiber_array.get_fiber(lidx))
        curr_count.append(1)
        curr_indices.append(list([lidx]))
        
    converged = False
    iteration_count = 0
    
    while not converged:
        print(f"<{os.path.basename(__file__)}> ITERATION: {iteration_count} SUM FIBER COUNTS: {np.sum(np.array(curr_count))}")
        print(f"<{os.path.basename(__file__)}> number indices {len(curr_indices)}")
        
        # fiber data structures for output of this iteration
        next_fibers = list()
        next_count = list()
        next_indices = list()
        
        # information for this iteration
        done = np.zeros(current_fiber_array.number_of_fibers)
        fiber_indices = list(range(0, current_fiber_array.number_of_fibers))

        # if the maximum number of fibers have been combined, stop averaging this fiber
        done[np.nonzero(np.array(curr_count) >= cluster_max)] = 1
        
        # pairwise distance matrix
        if have_joblib:
            distances = Parallel(n_jobs=n_jobs, verbose=1)(
                delayed(similarity.fiber_distance)(
                current_fiber_array.get_fiber(lidx),
                current_fiber_array,
                0, 'Hausdorff')
                for lidx in fiber_indices)
            distances = np.array(distances)
        else:
            distances = \
                np.zeros(
                (current_fiber_array.number_of_fibers,
                 current_fiber_array.number_of_fibers))
            for lidx in fiber_indices:
                distances[lidx, :] = \
                    similarity.fiber_distance(
                        current_fiber_array.get_fiber(lidx),
                        current_fiber_array, 0, 'Hausdorff')

        # distances to self are not of interest
        for lidx in fiber_indices:
            distances[lidx,lidx] = np.inf
        
        # sort the pairwise distances. 
        distances_flat = distances.flatten()
        pair_order = np.argsort(distances_flat)

        print(f"<{os.path.basename(__file__)}> DISTANCE MIN: {distances_flat[pair_order[0]]} DISTANCE COUNT: {distances.shape}")

        # if the smallest distance is greater or equal to the
        # threshold, we have converged
        if distances_flat[pair_order[0]] >= fiber_distance_threshold:
            converged = True
            print(f"<{os.path.basename(__file__)}> CONVERGED")
            break
        else:
            print(f"<{os.path.basename(__file__)}> NOT CONVERGED")
            
        # loop variables
        idx = 0
        pair_idx = pair_order[idx]
        number_of_fibers = distances.shape[0]
        number_averages = 0
        
        # combine nearest neighbors unless done, until hit threshold
        while distances_flat[pair_idx] < fiber_distance_threshold:
            # find the fiber indices corresponding to this pairwise distance
            # use div and mod
            f_row = pair_idx / number_of_fibers
            f_col = pair_idx % number_of_fibers

            # check if this neighbor pair can be combined
            combine = (not done[f_row]) and (not done[f_col])
            if combine :
                done[f_row] += 1
                done[f_col] += 1
                # weighted average of the fibers (depending on how many each one represents)
                next_fibers.append(
                    (curr_fibers[f_row] * curr_count[f_row] + \
                     curr_fibers[f_col] *curr_count[f_col]) / \
                    (curr_count[f_row] + curr_count[f_col]))
                # this was the regular average
                #next_fibers.append((curr_fibers[f_row] + curr_fibers[f_col])/2)
                next_count.append(curr_count[f_row] + curr_count[f_col])
                number_averages += 1
                #next_indices.append(list([curr_indices[f_row], curr_indices[f_col]]))
                next_indices.append(list(curr_indices[f_row] + curr_indices[f_col]))
                
            # increment for the loop
            idx += 1
            pair_idx = pair_order[idx]

        # copy through any unvisited (already converged) fibers
        unvisited = np.nonzero(done==0)[0]
        for fidx in unvisited:
            next_fibers.append(curr_fibers[fidx])
            next_count.append(curr_count[fidx])
            next_indices.append(curr_indices[fidx])
            
        # set up for next iteration
        curr_fibers = next_fibers
        curr_count = next_count
        curr_indices = next_indices
        iteration_count += 1

        # set up array for next iteration distance computation
        current_fiber_array = fibers.FiberArray()    
        current_fiber_array.number_of_fibers = len(curr_fibers)
        current_fiber_array.points_per_fiber = points_per_fiber
        dims = [current_fiber_array.number_of_fibers, current_fiber_array.points_per_fiber]
        # fiber data
        current_fiber_array.fiber_array_r = np.zeros(dims)
        current_fiber_array.fiber_array_a = np.zeros(dims)
        current_fiber_array.fiber_array_s = np.zeros(dims)
        curr_fidx = 0
        for curr_fib in curr_fibers:
            current_fiber_array.fiber_array_r[curr_fidx] = curr_fib.r
            current_fiber_array.fiber_array_a[curr_fidx] = curr_fib.a
            current_fiber_array.fiber_array_s[curr_fidx] = curr_fib.s
            curr_fidx += 1

        print(f"<{os.path.basename(__file__)}> SUM FIBER COUNTS: {np.sum(np.array(curr_count))} SUM DONE FIBERS: {np.sum(done)}")
        print(f"<{os.path.basename(__file__)}> MAX COUNT: {np.max(np.array(curr_count))} AVGS THIS ITER:  {number_averages}")

    # when converged, convert output to polydata
    outpd = current_fiber_array.convert_to_polydata()

    # color output by the number of fibers that each output fiber corresponds to
    outcolors = vtk.vtkFloatArray()
    outcolors.SetName('FiberTotal')
    for count in curr_count:
        outcolors.InsertNextTuple1(count)
    outpd.GetCellData().SetScalars(outcolors)

    # also color the input pd by output cluster number
    cluster_numbers = np.zeros(original_number_of_fibers)
    cluster_count = np.zeros(original_number_of_fibers)
    cluster_idx = 0
    for index_list in curr_indices:
        indices = np.array(index_list).astype(int)
        cluster_numbers[indices] = cluster_idx
        cluster_count[indices] = curr_count[cluster_idx]
        cluster_idx += 1
    outclusters =  vtk.vtkFloatArray()
    outclusters.SetName('ClusterNumber')
    for cluster in cluster_numbers:
        outclusters.InsertNextTuple1(cluster)
    inpd.GetCellData().AddArray(outclusters)
    inpd.GetCellData().SetActiveScalars('ClusterNumber')

    return outpd, np.array(curr_count), inpd, cluster_numbers, cluster_count
    

    
def laplacian_of_gaussian(inpd, fiber_distance_sigma = 25, points_per_fiber=30, n_jobs=2, upper_thresh=30):
    """ Filter nearby fibers, using LoG weights.
    
    The pairwise fiber distance matrix is computed, then fibers are
    averaged with their neighbors using LoG weighting.  This is
    essentially a fiber subtraction operation, giving vectors pointing
    from the center fiber under the kernel, to all nearby fibers. Thus
    the output of this operation is not a fiber, but we compute
    properties of the output that might be interesting and related to
    fibers. We summarize the result using the average vector at each
    fiber point (output is its magnitude, similar to edge
    strength). The covariance of the vectors is also
    investigated. This matrix would be spherical in an isotropic
    region such as a tract center (tube/line detector), or planar in a
    sheetlike tract (sheet detector).

    The equation is: (1-d^2/sigma^2) exp(-d^2/(2*sigma^2)), and
    weights are normalized in the neighborhood (weighted averaging).
    """

    sigmasq = fiber_distance_sigma * fiber_distance_sigma
    
    # polydata to array conversion, fixed-length fiber representation
    fiber_array = fibers.FiberArray()
    fiber_array.points_per_fiber = points_per_fiber
    fiber_array.convert_from_polydata(inpd)

    fiber_indices = list(range(0, fiber_array.number_of_fibers))

    # pairwise distance matrix
    if have_joblib:
        distances = Parallel(n_jobs=n_jobs, verbose=1)(
            delayed(similarity.fiber_distance)(
            fiber_array.get_fiber(lidx),
            fiber_array,
            0, 'Hausdorff')
            for lidx in fiber_indices)
        distances = np.array(distances)
    else:
        distances = \
            np.zeros(
            (fiber_array.number_of_fibers,
             fiber_array.number_of_fibers))
        for lidx in fiber_indices:
            distances[lidx, :] = \
                similarity.fiber_distance(
                    fiber_array.get_fiber(lidx),
                    fiber_array, 0)

    # fiber list data structure initialization for easy fiber averaging
    fiber_list = list()
    for lidx in range(0, fiber_array.number_of_fibers):
        fiber_list.append(fiber_array.get_fiber(lidx))

    filter_vectors = list()
    filter_vector_magnitudes = list()
    filter_confidences = list()
    
    # gaussian smooth all fibers using local neighborhood
    for fidx in fiber_indices:
        if (fidx % 100) == 0:
            print(f'{fidx} / {fiber_array.number_of_fibers}')

        current_fiber = fiber_list[fidx]

        # find indices of all nearby fibers
        # this includes the center fiber under the kernel
        indices = np.nonzero(distances[fidx] < upper_thresh)[0]
        local_fibers = list()
        local_weights = list()

        for idx in indices:
            dist = distances[fidx][idx]
            # compute filter kernel weights
            weight = np.exp(-(dist*dist)/sigmasq)
            #weight = (1 - (dist*dist)/sigmasq) * np.exp(-(dist*dist)/(2*sigmasq))
            local_fibers.append(fiber_list[idx])
            local_weights.append(weight)

        # actually perform the weighted average
        #mean_weight = np.mean(np.array(local_weights))
        #out_weights = local_weights[0]
        #for weight in local_weights[1:]:
        #    out_weights += weight
        # the weights must sum to 0 for LoG
        # (response in constant region is 0)
        #mean_weight = out_weights / len(local_weights)
        #local_normed_weights = list()
        #for weight in local_weights:
        #    local_normed_weights.append(weight - mean_weight)

        #match_fiber = local_fibers[0]
        #out_vector = local_fibers[0] * local_normed_weights[0]
        idx = 0
        for fiber in local_fibers:
            #out_vector += fiber
            # ensure fiber ordering by matching to current fiber only
            # otherwise the order is undefined after fiber subtraction
            matched_fiber = current_fiber.match_order(fiber)
            #filtered_fiber = matched_version * local_normed_weights[idx]
            #filtered_fiber = matched_version * local_weights[idx]
            if idx == 0:
                out_vector = fibers.Fiber()
                out_vector.points_per_fiber = points_per_fiber
                out_vector.r = np.zeros(points_per_fiber)
                out_vector.a = np.zeros(points_per_fiber)
                out_vector.s = np.zeros(points_per_fiber)
            #filtered_fiber = match_fiber.match_order(fiber)
            #out_vector.r = (out_vector.r + matched_fiber.r) * local_weights[idx]
            #out_vector.a = (out_vector.a + matched_fiber.a) * local_weights[idx]
            #out_vector.s = (out_vector.s + matched_fiber.s) * local_weights[idx]
            out_vector.r += (current_fiber.r - matched_fiber.r) * local_weights[idx]
            out_vector.a += (current_fiber.a - matched_fiber.a) * local_weights[idx]
            out_vector.s += (current_fiber.s - matched_fiber.s) * local_weights[idx]
            idx += 1

        total_weights = np.sum(np.array(local_weights))
        out_vector = out_vector / total_weights       

        filter_vectors.append(out_vector)
        filter_confidences.append(total_weights)

        filter_vector_magnitudes.append(np.sqrt(\
                np.multiply(out_vector.r, out_vector.r) + \
                    np.multiply(out_vector.a, out_vector.a) + \
                    np.multiply(out_vector.s, out_vector.s)))
        #filter_vector_magnitudes.append(np.sum(out_vector.r))


    # output a new pd!!!!
    # with fixed length fibers. and the new vector field.
    # output the vectors from the filtering
    outpd = fiber_array.convert_to_polydata()
    vectors = vtk.vtkFloatArray()
    vectors.SetName('FiberDifferenceVectors')
    vectors.SetNumberOfComponents(3)
    for vec in filter_vectors:
        for idx in range(points_per_fiber):
            vectors.InsertNextTuple3(vec.r[idx],vec.a[idx],vec.s[idx])
    magnitudes = vtk.vtkFloatArray()
    magnitudes.SetName('FiberDifferenceMagnitudes')
    magnitudes.SetNumberOfComponents(1)
    for mag in filter_vector_magnitudes:
        for idx in range(points_per_fiber):
            magnitudes.InsertNextTuple1(mag[idx])
    confidences = vtk.vtkFloatArray()
    confidences.SetName('FiberDifferenceConfidences')
    confidences.SetNumberOfComponents(1)
    for mag in filter_confidences:
        for idx in range(points_per_fiber):
            confidences.InsertNextTuple1(mag)
         
    outpd.GetPointData().AddArray(vectors)
    outpd.GetPointData().SetActiveVectors('FiberDifferenceVectors')

    outpd.GetPointData().AddArray(confidences)
    outpd.GetPointData().SetActiveScalars('FiberDifferenceConfidences')

    outpd.GetPointData().AddArray(magnitudes)
    outpd.GetPointData().SetActiveScalars('FiberDifferenceMagnitudes')

    # color by the weights or "local density"
    # color output by the number of fibers that each output fiber corresponds to
    #outcolors = vtk.vtkFloatArray()
    #outcolors.SetName('KernelDensity')
    #for weight in next_weights:
    #    outcolors.InsertNextTuple1(weight)
    #inpd.GetCellData().AddArray(outcolors)
    #inpd.GetCellData().SetActiveScalars('KernelDensity')
    #outcolors = vtk.vtkFloatArray()
    #outcolors.SetName('EdgeMagnitude')
    #for magnitude in filter_vector_magnitudes:
    #    outcolors.InsertNextTuple1(magnitude)
    #inpd.GetCellData().AddArray(outcolors)
    #inpd.GetCellData().SetActiveScalars('EdgeMagnitude')

    return outpd, np.array(filter_vector_magnitudes)

def pd_to_array(inpd, dims=225):
    count_vol = np.ndarray([dims,dims,dims])
    ptids = vtk.vtkIdList()
    points = inpd.GetPoints()
    data_vol = []    
    # check for cell data
    cell_data = inpd.GetCellData().GetScalars()
    if cell_data:
        data_vol = np.ndarray([dims,dims,dims])
    # loop over lines
    inpd.GetLines().InitTraversal()
    print(f"<{os.path.basename(__file__)}> Input number of points: {points.GetNumberOfPoints()} lines: {inpd.GetNumberOfLines()}")
    # loop over all lines
    for lidx in range(0, inpd.GetNumberOfLines()):
        # progress
        #if verbose:
        #    if lidx % 1 == 0:
        #        print(f"<{os.path.basename(__file__)}> Line: {lidx} / {inpd.GetNumberOfLines()}")
        inpd.GetLines().GetNextCell(ptids)
        num_points = ptids.GetNumberOfIds()
        for pidx in range(0, num_points):
            point = points.GetPoint(ptids.GetId(pidx))
            # center so that 0,0,0 moves to 100,100,100
            point = np.round(np.array(point) + 110)
            count_vol[point[0], point[1], point[2]] += 1
            if cell_data:
                data_vol[point[0], point[1], point[2]] += cell_data.GetTuple(lidx)[0]
    return count_vol, data_vol
    
def array_to_vtk(inarray, name='from_numpy'):
    vol = vtk.vtkImageData()
    dims = inarray.shape
    vol.SetDimensions(dims[0], dims[1], dims[2])
    vol.SetOrigin(0,0,0)
    #vol.SetSpacing(gridSpacing,gridSpacing,gridSpacing)
    sc = vtk.vtkShortArray()
    sc.SetNumberOfValues(dims[0] * dims[1] * dims[2])
    sc.SetNumberOfComponents(1)
    sc.SetName(name)
    for ii,tmp in enumerate(inarray.flatten()):
        sc.SetValue(ii,round((np.abs(tmp))*100))
    vol.GetPointData().SetScalars(sc)
    return vol
    

def measure_line_lengths(inpd):
    ptids = vtk.vtkIdList()
    points = inpd.GetPoints()
    output_lengths = np.zeros(inpd.GetNumberOfLines())
    # loop over lines
    inpd.GetLines().InitTraversal()
    print(f"<{os.path.basename(__file__)}> Input number of points: {points.GetNumberOfPoints()} lines: {inpd.GetNumberOfLines()}")
    # loop over all lines
    for lidx in range(0, inpd.GetNumberOfLines()):
        # progress
        #if verbose:
        #    if lidx % 1 == 0:
        #        print(f"<{os.path.basename(__file__)}> Line: {lidx} / {inpd.GetNumberOfLines()}")
        inpd.GetLines().GetNextCell(ptids)
        output_lengths[lidx] = ptids.GetNumberOfIds()
    return(output_lengths)

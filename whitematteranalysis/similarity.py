# -*- coding: utf-8 -*-

import math
import os
import sys

import numpy as np
import vtk

from . import fibers

sys.setrecursionlimit(1000000)

def distance_to_similarity(distance, sigmasq=100):

    # compute the similarities using Gaussian kernel
    similarities = np.exp(-distance / (sigmasq))

    return similarities

def fiber_distance(fiber, fiber_array, threshold=0, distance_method='MeanSquared', fiber_landmarks=None, landmarks=None, sigmasq=6400, bilateral=False):
    """
    Find pairwise fiber distance from fiber to all fibers in fiber_array.
    The Mean and MeanSquared distances are the average distance per
    fiber point, to remove scaling effects (dependence on number of
    points chosen for fiber parameterization). The Hausdorff distance
    is the maximum distance between corresponding points.
    input fiber should be class Fiber. fibers should be class FiberArray
    """

    # get fiber in reverse point order, equivalent representation
    fiber_equiv = fiber.get_equivalent_fiber()

    # compute pairwise fiber distances along fibers
    distance_1 = _fiber_distance_internal_use(fiber, fiber_array, threshold, distance_method, fiber_landmarks, landmarks, sigmasq)
    distance_2 = _fiber_distance_internal_use(fiber_equiv, fiber_array, threshold, distance_method, fiber_landmarks, landmarks, sigmasq)
        
    # choose the lowest distance, corresponding to the optimal fiber
    # representation (either forward or reverse order)
    if distance_method == 'StrictSimilarity':
        # for use in laterality
        # this is the product of all similarity values along the fiber
        distance = np.maximum(distance_1, distance_2)
    else:
        distance = np.minimum(distance_1, distance_2)

    if bilateral:
        fiber_reflect = fiber.get_reflected_fiber()
        # call this function again with the reflected fiber. Do NOT reflect again (bilateral=False) to avoid infinite loop.
        distance_reflect = fiber_distance(fiber_reflect, fiber_array, threshold, distance_method, fiber_landmarks, landmarks, sigmasq, bilateral=False)
        # choose the best distance, corresponding to the optimal fiber
        # representation (either reflected or not)
        if distance_method == 'StrictSimilarity':
            # this is the product of all similarity values along the fiber
            distance = np.maximum(distance, distance_reflect)
        else:
            distance = np.minimum(distance, distance_reflect)
        
    return distance

def fiber_distance_oriented(fiber, fiber_array, threshold=0, distance_method='MeanSquared', fiber_landmarks=None, landmarks=None, sigmasq=6400, bilateral=False):
    """Find pairwise fiber distance from fiber to all fibers in fiber_array, without testing both fiber orders.
    This is to be used if the fibers are ordered purposely, such as

    according to the endpoint nearest a cortical point.  The Mean and
    MeanSquared distances are the average distance per fiber point, to
    remove scaling effects (dependence on number of points chosen for
    fiber parameterization). The Hausdorff distance is the maximum
    distance between corresponding points.  input fiber should be class
    Fiber. fibers should be class FiberArray
    """

    # compute pairwise fiber distances along fibers
    distance = _fiber_distance_internal_use(fiber, fiber_array, threshold, distance_method, fiber_landmarks, landmarks, sigmasq)

    if bilateral:
        fiber_reflect = fiber.get_reflected_fiber()
        # call this function again with the reflected fiber. Do NOT reflect again (bilateral=False) to avoid infinite loop.
        distance_reflect = fiber_distance_oriented(fiber_reflect, fiber_array, threshold, distance_method, fiber_landmarks, landmarks, sigmasq, bilateral=False)
        # choose the best distance, corresponding to the optimal fiber
        # representation (either reflected or not)
        if distance_method == 'StrictSimilarity':
            # this is the product of all similarity values along the fiber
            distance = np.maximum(distance, distance_reflect)
        else:
            distance = np.minimum(distance, distance_reflect)
        
    return distance

def _fiber_distance_internal_use(fiber, fiber_array, threshold=0, distance_method='MeanSquared', fiber_landmarks=None, landmarks=None, sigmasq=None):
    """ Compute the total fiber distance from one fiber to an array of
    many fibers.
    This function does not handle equivalent fiber representations,
    for that use fiber_distance, above.
    """

    if distance_method == 'Landmarks':
        return _fiber_distance_internal_landmarks(fiber, fiber_array, fiber_landmarks, landmarks)
    if landmarks is not None:
        print("ERROR: Please use distance method Landmarks to compute landmark distances")
        
    # compute the distance from this fiber to the array of other fibers
    ddx = fiber_array.fiber_array_r - fiber.r
    ddy = fiber_array.fiber_array_a - fiber.a
    ddz = fiber_array.fiber_array_s - fiber.s

    dx = np.square(ddx)
    dy = np.square(ddy)
    dz = np.square(ddz)

    # sum dx dx dz at each point on the fiber and sqrt for threshold
    #distance = np.sqrt(dx + dy + dz)
    distance = dx + dy + dz

    # threshold if requested
    if threshold:
        # set values less than threshold to 0
        distance = distance - threshold*threshold
        idx = np.nonzero(distance < 0)
        distance[idx] = 0

    if distance_method == 'Mean':
        # sum along fiber
        distance = np.sum(np.sqrt(distance), 1)
        # Remove effect of number of points along fiber (mean)
        npts = float(fiber_array.points_per_fiber)
        distance = distance / npts
        # for consistency with other methods we need to square this value
        distance = np.square(distance)
    elif distance_method == 'Hausdorff':
        # take max along fiber
        distance = np.max(distance, 1)
    elif distance_method == 'MeanSquared':
        # sum along fiber
        distance = np.sum(distance, 1)
        # Remove effect of number of points along fiber (mean)
        npts = float(fiber_array.points_per_fiber)
        distance = distance / npts
    elif distance_method == 'StrictSimilarity':
        # for use in laterality
        # this is the product of all similarity values along the fiber
        # not truly a distance but it's easiest to compute here in this function
        # where we have all distances along the fiber
        #print "distance range :", np.min(distance), np.max(distance)
        distance = distance_to_similarity(distance, sigmasq)
        #print "similarity range :", np.min(distance), np.max(distance)        
        distance = np.prod(distance, 1)
        #print "overall similarity range:", np.min(distance), np.max(distance)
    elif distance_method == 'Mean_shape':
        
        # sum along fiber
        distance_square = distance
        distance = np.sqrt(distance_square)

        d = np.sum(distance, 1)
        # Remove effect of number of points along fiber (mean)
        npts = float(fiber_array.points_per_fiber)
        d = np.divide(d, npts)
        # for consistency with other methods we need to square this value
        d = np.square(d)

        distance_endpoints = (distance[:,0] + distance[:,npts-1])/2

        for i in np.linspace(0,np.size(distance,0)-1,np.size(distance,0)):
            for j in np.linspace(0,np.size(distance,1)-1,np.size(distance,1)):
                if distance[i,j] == 0:
                    distance[i,j] = 1
        ddx = np.divide(ddx,distance)
        ddy = np.divide(ddy,distance)
        ddz = np.divide(ddz,distance)
        #print ddx*ddx+ddy*ddy+ddz*ddz
        npts = float(fiber_array.points_per_fiber)
        angles = np.zeros([(np.size(distance))/npts,npts*(npts+1)/2])
        s = 0
        n = np.linspace(0,npts-1,npts)
        for i in n:
            m = np.linspace(0,i,i+1)
            for j in m:
                angles[:,s] = (ddx[:,i]-ddx[:,j])*(ddx[:,i]-ddx[:,j]) + (ddy[:,i]-ddy[:,j])*(ddy[:,i]-ddy[:,j]) + (ddz[:,i]-ddz[:,j])*(ddz[:,i]-ddz[:,j])
                s = s+1
        angles = (np.sqrt(angles))/2
        angle = np.max(angles,1)
        #print angle.max()
        
        distance = 0.5*d + 0.4*d/(0.5+0.5*(1-angle*angle)) + 0.1*distance_endpoints
        
    else:
        print(f"<{os.path.basename(__file__)}> throwing Exception. Unknown input distance method (typo?):", distance_method)
        raise Exception("unknown distance method")
        
    
    return distance

def rectangular_frechet_distances(input_vtk_polydata_m,input_vtk_polydata_n):
    # Computes distance matrix (nxm) for all n+m fibers in input
    # polydata. each fiber in input_polydata_n is compared to each fiber
    # in input_polydata_m.

    number_of_lines_n = input_vtk_polydata_n.GetNumberOfLines()
    number_of_lines_m = input_vtk_polydata_m.GetNumberOfLines()
    all_fibers_n = range(0,number_of_lines_n)
    all_fibers_m = range(0,number_of_lines_m)
    distances = np.zeros([number_of_lines_n,number_of_lines_m])

    input_vtk_polydata_n.GetLines().InitTraversal()
    line1_ptids = vtk.vtkIdList()

    inpoints1 = input_vtk_polydata_n.GetPoints()
    inpoints2 = input_vtk_polydata_m.GetPoints()

    for lidx1 in all_fibers_n:
        input_vtk_polydata_n.GetLines().GetNextCell(line1_ptids)

        input_vtk_polydata_m.GetLines().InitTraversal()
        line2_ptids = vtk.vtkIdList()

        for lidx2 in all_fibers_m:
            input_vtk_polydata_m.GetLines().GetNextCell(line2_ptids)

            distances[lidx1,lidx2] = _frechet_distance_internal_use(inpoints1,inpoints2,line1_ptids,line2_ptids)

    return distances

def pairwise_frechet_distances(input_vtk_polydata_n,input_vtk_polydata_m):
    # Computes distance matrix (nxn) for all n fibers in input polydata.

    number_of_lines_n = input_vtk_polydata_n.GetNumberOfLines()
    number_of_lines_m = input_vtk_polydata_m.GetNumberOfLines()
    all_fibers_n = range(0,number_of_lines_n)
    all_fibers_m = range(0,number_of_lines_m)
    distances = np.zeros([number_of_lines_n,number_of_lines_m])

    input_vtk_polydata_n.GetLines().InitTraversal()
    line1_ptids = vtk.vtkIdList()

    inpoints1 = input_vtk_polydata_n.GetPoints()
    inpoints2 = input_vtk_polydata_m.GetPoints()

    for lidx1 in all_fibers_n:
        input_vtk_polydata_n.GetLines().GetNextCell(line1_ptids)

        input_vtk_polydata_m.GetLines().InitTraversal()
        line2_ptids = vtk.vtkIdList()

        for lidx2 in all_fibers_n:
            input_vtk_polydata_m.GetLines().GetNextCell(line2_ptids)

            if lidx1 == lidx2:
                distances[lidx1,lidx2] = 0
            else:
                if distances[lidx1,lidx2] == 0:
                    distances[lidx1,lidx2] = _frechet_distance_internal_use(inpoints1,inpoints2,line1_ptids,line2_ptids)
                    distances[lidx2,lidx1] = distances[lidx1,lidx2]

    return distances


def Frechet_distances_2(input_vtk_polydata_n, input_vtk_polydata_m):
    # compute Frechet distance from an array of fibers to another array

    number_of_lines_m = input_vtk_polydata_m.GetNumberOfLines()
    number_of_lines_n = input_vtk_polydata_n.GetNumberOfLines()
    all_fibers_m = range(0,number_of_lines_m)
    all_fibers_n = range(0,number_of_lines_n)
    distances = np.zeros([number_of_lines_m,number_of_lines_n])
    
    #input_vtk_polydata2 = input_vtk_polydata
    
    input_vtk_polydata_m.GetLines().InitTraversal()
    line1_ptids = vtk.vtkIdList()
    
    inpoints1 = input_vtk_polydata_m.GetPoints()
    inpoints2 = input_vtk_polydata_n.GetPoints()
    
    for lidx1 in all_fibers_m:
        input_vtk_polydata_m.GetLines().GetNextCell(line1_ptids)
        
        input_vtk_polydata_n.GetLines().InitTraversal()
        line2_ptids = vtk.vtkIdList()
        
        for lidx2 in all_fibers_n:
            input_vtk_polydata_n.GetLines().GetNextCell(line2_ptids)

            distances[lidx1,lidx2] = _frechet_distance_internal_use(inpoints1,inpoints2,line1_ptids,line2_ptids)

    return distances

def _euc_dist(pt1,pt2):
    return math.sqrt((pt2[0]-pt1[0])*(pt2[0]-pt1[0])+(pt2[1]-pt1[1])*(pt2[1]-pt1[1])+(pt2[2]-pt1[2])*(pt2[2]-pt1[2]))

def _c(ca,i,j,P,Q):
    if ca[i,j] > -1:
        return ca[i,j]
    elif i == 0 and j == 0:
        ca[i,j] = _euc_dist(P[0],Q[0])
    elif i > 0 and j == 0:
        ca[i,j] = max(_c(ca,i-1,0,P,Q),_euc_dist(P[i],Q[0]))
    elif i == 0 and j > 0:
        ca[i,j] = max(_c(ca,0,j-1,P,Q),_euc_dist(P[0],Q[j]))
    elif i > 0 and j > 0:
        ca[i,j] = max(min(_c(ca,i-1,j,P,Q),_c(ca,i-1,j-1,P,Q),_c(ca,i,j-1,P,Q)),_euc_dist(P[i],Q[j]))
    else:
        ca[i,j] = float("inf")
    return ca[i,j]

def frechDist(P,Q):
    ca = np.ones((len(P),len(Q)))
    ca = np.multiply(ca,-1)
    return _c(ca,len(P)-1,len(Q)-1,P,Q)    

def _frechet_distance_internal_use(inpoints1,inpoints2,line1_ptids,line2_ptids):
    
    line1_length = line1_ptids.GetNumberOfIds()
    line2_length = line2_ptids.GetNumberOfIds()
    d1 = 0
    d2 = 0
    for i in range(0,9):
        d1 = d1 + _euc_dist(inpoints1.GetPoint(line1_ptids.GetId(i*5)),inpoints2.GetPoint(line2_ptids.GetId(i*5)))
        d2 = d2 + _euc_dist(inpoints1.GetPoint(line1_ptids.GetId(line1_length-1-i*5)),inpoints2.GetPoint(line2_ptids.GetId(line2_length-1-i*5)))
    d1 = d1/9
    d2 = d2/9
    line1_points = np.zeros([int(round(line1_length/80))-1,3])
    line2_points = np.zeros([int(round(line2_length/80))-1,3])
    all_points1 = range(1,int(round((line1_length)/80)))
    all_points2 = range(1,int(round((line2_length)/80)))
    for i1 in all_points1:
        ptidx1 = line1_ptids.GetId(80*i1)
        line1_points[i1-1] = inpoints1.GetPoint(ptidx1)
    for i2 in all_points2:
        ptidx2 = line2_ptids.GetId(80*i2)
        line2_points[i2-1] = inpoints2.GetPoint(ptidx2)
    return (frechDist(line1_points,line2_points)*2+d1+d2)/4

def _fiber_distance_internal_landmarks(fiber, fiber_array, fiber_landmarks, landmarks):
    # compute the distance from this fiber to the array of other fibers
    # where distance is according to differences in landmark distance
    [n_fibers, n_landmarks, dims] = landmarks.shape
    #print "*******************************"
    #print "Using landmarks. N fibers:", n_fibers, "N landmarks:", n_landmarks, "dimensions (3):", dims
    #print "*******************************"
    diffs = np.zeros((n_fibers, n_landmarks))
    for lidx in range(n_landmarks):
        # compute fiber array landmark distances to this landmark
        dx = np.subtract(fiber_array.fiber_array_r.T, landmarks[:,lidx,0]).T
        dy = np.subtract(fiber_array.fiber_array_a.T, landmarks[:,lidx,1]).T
        dz = np.subtract(fiber_array.fiber_array_s.T, landmarks[:,lidx,2]).T
        
        dx = np.power(dx, 2)
        dy = np.power(dy, 2)
        dz = np.power(dz, 2)
        landmark_distance = np.sqrt(dx + dy + dz)
        #print "lm dist array shape", landmark_distance.shape
        
        del dx 
        del dy 
        del dz
        # compute for individual fiber
        dx = fiber.r - fiber_landmarks[lidx,0]
        dy = fiber.a - fiber_landmarks[lidx,1]
        dz = fiber.s - fiber_landmarks[lidx,2]
        dx = np.power(dx, 2)
        dy = np.power(dy, 2)
        dz = np.power(dz, 2)
        fiber_landmark_distance = np.sqrt(dx + dy + dz)
        #print "lm dist fiber shape", fiber_landmark_distance.shape
        del dx 
        del dy 
        del dz
        # difference between LD at all points on fiber
        ld_diff = landmark_distance - fiber_landmark_distance
        # summarize for this landmark
        diffs[:,lidx] = np.mean(np.multiply(ld_diff, ld_diff), 1)
        
    distance = np.sqrt(np.mean(diffs, 1))
    return distance

def total_similarity_for_laterality(fiber, fiber_array, reflect, threshold, sigmasq):
    """ Convenience function for fiber similarity needed for laterality.
    Optionally reflects fiber, returns total similarity to all input
    fibers in fiber_array (these are from one hemisphere). Also
    returns distances array for further analysis.
    """

    if reflect:
        fiber = fiber.get_reflected_fiber()
    similarity = fiber_distance(fiber, fiber_array, threshold=threshold, distance_method='StrictSimilarity', sigmasq=sigmasq)

    # compute fiber similarity total (eg to all fibers in a hemisphere)
    total_similarity = np.sum(similarity)

    # Do not include the self fiber in the computation (subtract 1.0)
    #if ~reflect:
    #    total_similarity = total_similarity - 1
    # DO not do the above. Prevents a completely symmetric brain from having all LIs of 0
    # because the self fiber is missing.
        
    return total_similarity

# this must be a function to allow pickling by Parallel
def total_similarity_and_distances(fiber, fiber_array, reflect, threshold, sigmasq, distance_method='MeanSquared'):
    """ Convenience function for fiber similarity needed for laterality.
    Optionally reflects fiber, returns total similarity to all input
    fibers in fiber_array (these are from one hemisphere). Also
    returns distances array for further analysis.
    """

    if reflect:
        fiber = fiber.get_reflected_fiber()

    distance = fiber_distance(fiber, fiber_array, threshold, distance_method)
    similarity = distance_to_similarity(distance, sigmasq)

    # compute fiber similarity total (eg to all fibers in a hemisphere)
    total_similarity = np.sum(similarity)

    return total_similarity, distance


# this must be a function to allow pickling by Parallel
def total_similarity(fiber, fiber_array, threshold, sigmasq, distance_method='MeanSquared'):
    """ Convenience function that just returns total fiber similarity.
    Returns a single scalar value: total similarity to all input
    fibers in fiber_array (these can be from one hemisphere or the
    whole brain, etc).
    """
    
    distance = fiber_distance(fiber, fiber_array, threshold, distance_method)
    similarity = distance_to_similarity(distance, sigmasq)

    # compute fiber similarity total (eg to all fibers in a hemisphere)
    total_similarity = np.sum(similarity)

    return total_similarity

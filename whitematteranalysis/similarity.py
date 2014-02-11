import numpy

import fibers

def distance_to_similarity(distance, sigmasq=100):

    # compute the similarities using Gaussian kernel
    similarities = numpy.exp(-distance / (sigmasq))

    return similarities

def fiber_distance(fiber, fiber_array, threshold=0, distance_method='MeanSquared', fiber_landmarks=None, landmarks=None, sigmasq=6400):
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
        distance = numpy.maximum(distance_1, distance_2)
    else:
        distance = numpy.minimum(distance_1, distance_2)
    
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
        print "ERROR: Please use distance method Landmarks to compute landmark distances"
        
    # compute the distance from this fiber to the array of other fibers
    dx = fiber_array.fiber_array_r - fiber.r
    dy = fiber_array.fiber_array_a - fiber.a
    dz = fiber_array.fiber_array_s - fiber.s

    dx = numpy.square(dx)
    dy = numpy.square(dy)
    dz = numpy.square(dz)

    # sum dx dx dz at each point on the fiber and sqrt for threshold
    #distance = numpy.sqrt(dx + dy + dz)
    distance = dx + dy + dz

    # threshold if requested
    if threshold:
        # set values less than threshold to 0
        distance = distance - threshold*threshold
        idx = numpy.nonzero(distance < 0)
        distance[idx] = 0

    if distance_method == 'Mean':
        # sum along fiber
        distance = numpy.sum(numpy.sqrt(distance), 1)
        # Remove effect of number of points along fiber (mean)
        npts = float(fiber_array.points_per_fiber)
        distance = distance / npts
        # for consistency with other methods we need to square this value
        distance = numpy.square(distance)
    elif distance_method == 'Hausdorff':
        # take max along fiber
        distance = numpy.max(distance, 1)
    elif distance_method == 'MeanSquared':
        # sum along fiber
        distance = numpy.sum(distance, 1)
        # Remove effect of number of points along fiber (mean)
        npts = float(fiber_array.points_per_fiber)
        distance = distance / npts
    elif distance_method == 'StrictSimilarity':
        # for use in laterality
        # this is the product of all similarity values along the fiber
        # not truly a distance but it's easiest to compute here in this function
        # where we have all distances along the fiber
        #print "distance range :", numpy.min(distance), numpy.max(distance)
        distance = distance_to_similarity(distance, sigmasq)
        #print "similarity range :", numpy.min(distance), numpy.max(distance)        
        distance = numpy.prod(distance, 1)
        #print "overall similarity range:", numpy.min(distance), numpy.max(distance)
        
    else:
        print "<similarity.py> throwing Exception. Unknown input distance method (typo?):", distance_method
        raise Exception("unknown distance method")
        
    
    return distance

def _fiber_distance_internal_landmarks(fiber, fiber_array, fiber_landmarks, landmarks):
    # compute the distance from this fiber to the array of other fibers
    # where distance is according to differences in landmark distance
    [n_fibers, n_landmarks, dims] = landmarks.shape
    #print "*******************************"
    #print "Using landmarks. N fibers:", n_fibers, "N landmarks:", n_landmarks, "dimensions (3):", dims
    #print "*******************************"
    diffs = numpy.zeros((n_fibers, n_landmarks))
    for lidx in range(n_landmarks):
        # compute fiber array landmark distances to this landmark
        dx = numpy.subtract(fiber_array.fiber_array_r.T, landmarks[:,lidx,0]).T
        dy = numpy.subtract(fiber_array.fiber_array_a.T, landmarks[:,lidx,1]).T
        dz = numpy.subtract(fiber_array.fiber_array_s.T, landmarks[:,lidx,2]).T
        
        dx = numpy.power(dx, 2)
        dy = numpy.power(dy, 2)
        dz = numpy.power(dz, 2)
        landmark_distance = numpy.sqrt(dx + dy + dz)
        #print "lm dist array shape", landmark_distance.shape
        
        del dx 
        del dy 
        del dz
        # compute for individual fiber
        dx = fiber.r - fiber_landmarks[lidx,0]
        dy = fiber.a - fiber_landmarks[lidx,1]
        dz = fiber.s - fiber_landmarks[lidx,2]
        dx = numpy.power(dx, 2)
        dy = numpy.power(dy, 2)
        dz = numpy.power(dz, 2)
        fiber_landmark_distance = numpy.sqrt(dx + dy + dz)
        #print "lm dist fiber shape", fiber_landmark_distance.shape
        del dx 
        del dy 
        del dz
        # difference between LD at all points on fiber
        ld_diff = landmark_distance - fiber_landmark_distance
        # summarize for this landmark
        diffs[:,lidx] = numpy.mean(numpy.multiply(ld_diff, ld_diff), 1)
        
    distance = numpy.sqrt(numpy.mean(diffs, 1))
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
    total_similarity = numpy.sum(similarity)

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
    total_similarity = numpy.sum(similarity)

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
    total_similarity = numpy.sum(similarity)

    return total_similarity

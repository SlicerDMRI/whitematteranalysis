import numpy

import fibers

def distance_to_similarity(distance, sigmasq=100):

    # compute the similarities using Gaussian kernel
    similarities = numpy.exp(-numpy.power(distance, 2) / (sigmasq))

    return similarities

def fiber_distance(fiber, fiber_array, threshold=0, distance_method='MeanSquared', fiber_landmarks=None, landmarks=None):
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
    distance_1 = _fiber_distance_internal_use(fiber, fiber_array, threshold, distance_method, fiber_landmarks, landmarks)
    distance_2 = _fiber_distance_internal_use(fiber_equiv, fiber_array, threshold, distance_method, fiber_landmarks, landmarks)

    # choose the lowest distance, corresponding to the optimal fiber
    # representation (either forward or reverse order)
    distance = numpy.minimum(distance_1, distance_2)

    if distance_method == 'Mean':
        # now find the average distance to remove effect of number of
        # points chosen to represent the fiber
        npts = fiber_array.points_per_fiber
        distance = distance / float(npts)
    elif distance_method == 'MeanSquared':
        npts = float(fiber_array.points_per_fiber)
        distance = distance / (npts * npts)
        distance = numpy.sqrt(distance)
        
    return distance

def _fiber_distance_internal_use(fiber, fiber_array, threshold=0, distance_method='MeanSquared', fiber_landmarks=None, landmarks=None):
    """ Compute the total fiber distance from one fiber to an array of
    many fibers.

    This function does not handle equivalent fiber representations,
    for that use fiber_distance, above.
    """

    if landmarks is not None:
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

    # compute the distance from this fiber to the array of other fibers
    dx = fiber_array.fiber_array_r - fiber.r
    dy = fiber_array.fiber_array_a - fiber.a
    dz = fiber_array.fiber_array_s - fiber.s

    dx = numpy.square(dx)
    dy = numpy.square(dy)
    dz = numpy.square(dz)

    # sum dx dx dz at each point on the fiber and sqrt for threshold
    distance = numpy.sqrt(dx + dy + dz)

    # threshold if requested
    if threshold:
        print "FIX thresholded distance computation"
        distance = numpy.maximum(distance - threshold, 0)   

    # note to self: this coding error produces very interesting colors
    # that are completely symmetric across hemispheres.
    #dx = numpy.maximum(numpy.abs(-dx-dxthresh),0)
    #dy = numpy.maximum(numpy.abs(-dy-dythresh),0)
    #dz = numpy.maximum(numpy.abs(-dz-dzthresh),0)

    if distance_method == 'Mean':
        # sum along fiber
        distance = numpy.sum(distance, 1)
    elif distance_method == 'Hausdorff':
        # take max along fiber
        distance = numpy.max(distance, 1)
    elif distance_method == 'MeanSquared':
        # sum along fiber
        distance = numpy.sum(numpy.multiply(distance,distance), 1)
        
    return distance


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

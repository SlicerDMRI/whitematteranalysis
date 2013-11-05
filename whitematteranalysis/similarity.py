import numpy

import fibers

def distance_to_similarity(distance, sigmasq=100):

    # compute the similarities using Gaussian kernel
    similarities = numpy.exp(-numpy.power(distance, 2) / (sigmasq))

    return similarities

def fiber_distance(fiber, fiber_array, threshold=0, distance_method='MeanSquared', fiber_landmarks=None, landmarks=None):
    """

    Find pairwise fiber distance from fiber to all fibers in fiber_array.

    The Mean and MeanSquared distances are the average distance per fiber unit length, to
    remove scaling effects (dependence on number of points chosen for
    fiber parameterization). The Hausdorff distance is the maximum distance between corresponding points.

    input fiber should be class Fiber. fibers should be class FiberArray

    """

    #nfib = fiber_array.number_of_fibers 
    #npts = fiber_array.points_per_fiber

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
            #l0 = numpy.tile(landmarks[:,lidx,0], (1,fiber_array.fiber_array_r.shape[2]))
            #l1 = numpy.tile(landmarks[:,lidx,1], (1,fiber_array.fiber_array_r.shape[2]))
            #l2 = numpy.tile(landmarks[:,lidx,2], (1,fiber_array.fiber_array_r.shape[2]))
            #dx = numpy.subtract(fiber_array.fiber_array_r, l0)
            #dy = numpy.subtract(fiber_array.fiber_array_a, l1)
            #dz = numpy.subtract(fiber_array.fiber_array_s, l2)

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

    #nfib = fiber_array.number_of_fibers 

    # like repmat. copy current line to full matrix size of all lines
    #fiber_r = numpy.tile(fiber.r, (nfib, 1))
    #fiber_a = numpy.tile(fiber.a, (nfib, 1))
    #fiber_s = numpy.tile(fiber.s, (nfib, 1))

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

def computeTotalFiberSimilarity(fiber_array, sigmasq, lidx, numLines, threshold):
    fiber = fiber_array.get_fiber(lidx)
    return total_similarity(fiber, fiber_array, threshold, sigmasq)

# this is not part of the class to allow pickling by Parallel     
def computeTotalFiberSimilarityOLD(fibers, sigmasq, lidx, numLines, threshold):
    # compare reflected fiber to all other fibers

    #print fibers.fiber_array_r[0]

    # like repmat. copy current line to full matrix size of all lines
    currx = -numpy.tile(fibers.fiber_array_r[lidx,:],(numLines,1))
    curry = numpy.tile(fibers.fiber_array_a[lidx,:],(numLines,1))
    currz = numpy.tile(fibers.fiber_array_s[lidx,:],(numLines,1))
    
    # do the same for the reverse order of current line
    # as the fiber can be equivalently represented in either order.
    currx2 = -numpy.tile(fibers.fiber_array_r[lidx,::-1],(numLines,1))
    curry2 = numpy.tile(fibers.fiber_array_a[lidx,::-1],(numLines,1))
    currz2 = numpy.tile(fibers.fiber_array_s[lidx,::-1],(numLines,1))

    # compute the distance (and sum along fiber for total distance)
    # handle forward order of fiber points
    dx = fibers.fiber_array_r - currx
    dy = fibers.fiber_array_a - curry
    dz = fibers.fiber_array_s - currz

    dx = numpy.power(dx, 2)
    dy = numpy.power(dy, 2)
    dz = numpy.power(dz, 2)
    
    # sum dx dx dz at each point and along fiber
    #d  = numpy.sum(dx+dy+dz, 1)
    d  = numpy.sqrt(dx+dy+dz)

    # threshold if requested
    d =  numpy.maximum(d-threshold,0)
    
    # sum along fiber
    d = numpy.sum(d, 1)

    # handle reverse order of fiber points
    dx2 = fibers.fiber_array_r - currx2
    dy2 = fibers.fiber_array_a - curry2
    dz2 = fibers.fiber_array_s - currz2
    
    dx2 = numpy.power(dx2, 2)
    dy2 = numpy.power(dy2, 2)
    dz2 = numpy.power(dz2, 2)

    # sum dx dx dz at each point and along fiber
    #d2  = numpy.sum(dx2+dy2+dz2, 1)
    d2  = numpy.sqrt(dx2+dy2+dz2)

    # threshold if requested
    d2 =  numpy.maximum(d2-threshold,0)
    
    # sum along fiber
    d2 = numpy.sum(d2, 1)

    # use the minimum total squared distance (this must be the one
    # with the best possible point correspondence, either forward or reverse order)
    dfinal = numpy.minimum(d,d2)

    # now find the average squared distance to remove effect of number of 
    # points chosen to represent the fiber
    lineLen = fibers.fiber_array_r.shape[1]
    #dfinal = dfinal/float(lineLen*lineLen)
    dfinal = dfinal/float(lineLen)

    # compute the similarities using Gaussian kernel
    #s  = numpy.exp(-dfinal/(sigmasq))
    s  = numpy.exp(-numpy.power(dfinal,2)/(sigmasq))
     
    # TEST hard threshold for debugging
    #s = (dfinal < 20)

    # compute fiber similarity totals for this hemisphere
    totalSimilarity= numpy.sum( s );

    return totalSimilarity


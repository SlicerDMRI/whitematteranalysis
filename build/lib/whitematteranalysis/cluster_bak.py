""" cluster.py

implementations of fiber clustering

"""

import vtk
import numpy
import colorsys
try:
    import scipy.cluster.vq
    import scipy.cluster.hierarchy
    USE_SCIPY = 1
except ImportError:
    USE_SCIPY = 0
    print "<cluster.py> Failed to import  scipy.cluster, cannot cluster."
    print "<cluster.py> Please install  scipy for this functionality."
try:
    import matplotlib.pyplot as plt
    USE_MATPLOTLIB = 1
except ImportError:
    USE_MATPLOTLIB = 0
    print "<cluster.py> Failed to import  matplotlib.pyplot, cannot plot."
    print "<cluster.py> Please install matplotlib for this functionality."
try:
    from joblib import Parallel, delayed
    USE_PARALLEL = 1
except ImportError:
    USE_PARALLEL = 0
    print "<cluster.py> Failed to import joblib, cannot multiprocess."
    print "<cluster.py> Please install joblib for this functionality."

import fibers
import similarity
import filter
import render


def hierarchical(input_polydata, number_of_clusters=300,
                 threshold=2, fcluster_threshold=0.4,
                 number_of_jobs=3):
    """ This is just a test. This is only a test.  Stay tuned."""

    distance_matrix = \
        _pairwise_distance_matrix(input_polydata, threshold,
                                  number_of_jobs)

    # convert to upper triangular [not needed]
    # (note should compute only this part.)
    #mask_matrix = numpy.ones(distance_matrix.shape)
    #mask_matrix = numpy.triu(mask_matrix)

    #distance_matrix_triu = distance_matrix(numpy.nonzero(mask_matrix))

    z_link = scipy.cluster.hierarchy.linkage(distance_matrix)
    cluster_idx = scipy.cluster.hierarchy.fcluster(z_link, fcluster_threshold)

    cluster_colors = vtk.vtkIntArray()
    cluster_colors.SetName('ClusterNumber')
    output_polydata = input_polydata
    for lidx in range(0, output_polydata.GetNumberOfLines()):
        cluster_colors.InsertNextTuple1(cluster_idx[lidx])

    output_polydata.GetCellData().AddArray(cluster_colors)
    output_polydata.GetCellData().SetActiveScalars('ClusterNumber')

    return output_polydata, cluster_idx

def spectral(input_polydata, number_of_clusters=300,
             number_of_eigenvectors=10, sigma=20, threshold=2,
             number_of_jobs=3, use_nystrom=False, nystrom_mask = None):

    """ Spectral clustering based on pairwise fiber affinity matrix.


    As in O'Donnell and Westin TMI 2007.

    Differences from that implementation: fiber distance is defined
    using constant-length fiber parameterization.

    """

    # test pd has lines first
    num_lines = input_polydata.GetNumberOfLines()
    print "Number of input lines for clustering:", num_lines

    if num_lines == 0:
        print "Cannot cluster with 0 lines."
        return []

    # if computing eigenvectors, etc. from a sample of the data
    if use_nystrom:
        # make sure it's an array for logic operation ~
        nystrom_mask = numpy.array(nystrom_mask)
        # make sure it's large enough
        test = sum(nystrom_mask) > 500
        if not test:
            print "<cluster.py> ERROR: Nystrom data mask is smaller than 500."
            raise AssertionError
        # make sure its size matches the polydata
        test = len(nystrom_mask) == input_polydata.GetNumberOfLines()
        if not test:
            print "<cluster.py> ERROR: Nystrom data mask size does not match polydata number of lines."
            raise AssertionError
        
        polydata_m = filter.mask(input_polydata, nystrom_mask)
        polydata_n = filter.mask(input_polydata, nystrom_mask == False)
        print polydata_m.GetNumberOfLines()
        print polydata_n.GetNumberOfLines()

        similarity_matrix = \
            _pairwise_similarity_matrix(polydata_m, threshold,
                                        sigma, number_of_jobs)
        similarity_matrix_rectangular = \
            _rectangular_similarity_matrix(polydata_n, polydata_m, threshold,
                                           sigma, number_of_jobs)

    else:
        similarity_matrix = \
            _pairwise_similarity_matrix(input_polydata, threshold,
                                    sigma, number_of_jobs)

    # do not normalize by row sums or column sums?
    # would be good to test this. don't need outliers to be closer to good data.

    # See the paper:
    # "Spectral Grouping Using the Nystrom Method"
    # (D^-1/2 W D^-1/2) V = V Lambda
    # embedding_i,j = V_i+i,j./sqrt(D_j,j)
    # norm cuts transform of similarity matrix
    if use_nystrom:
        pinv_A = numpy.linalg.pinv(similarity_matrix)
        row_sum_1 = numpy.sum(similarity_matrix, axis=0) + numpy.sum(similarity_matrix_rectangular.T, axis=0)

        print(row_sum_1.shape)
        print(numpy.sum(similarity_matrix_rectangular.T, axis=0).shape)
        print(pinv_A.shape)
        print(similarity_matrix_rectangular.shape)
        
        row_sum_2 = numpy.sum(similarity_matrix_rectangular, axis=0) + \
            numpy.dot(numpy.dot(numpy.sum(similarity_matrix_rectangular.T, axis=0), pinv_A), similarity_matrix_rectangular)

        print(row_sum_1.shape)
        print(numpy.min(row_sum_1))
        print(numpy.max(row_sum_1))
        
        print(row_sum_2.shape)
        print(numpy.min(row_sum_2))
        print(numpy.max(row_sum_2))
        
        dhat = numpy.divide(1, numpy.sqrt(row_sum_1))
        dhat_2 = numpy.divide(1, numpy.sqrt(row_sum_2))
        affinity_matrix = \
            numpy.multiply(similarity_matrix, numpy.outer(dhat, dhat.T))
        rectangular_affinity_matrix = \
            numpy.multiply(similarity_matrix_rectangular, numpy.outer(dhat, dhat_2.T))
    else:
        row_sum = numpy.sum(similarity_matrix, axis=0)
        dhat = numpy.divide(1, numpy.sqrt(row_sum))
        affinity_matrix = \
            numpy.multiply(similarity_matrix, numpy.outer(dhat, dhat.T))

    # embed
    print 'Embedding, calculating eigenvectors.'
    # this is slow
    e_values, e_vectors = numpy.linalg.eigh(affinity_matrix)

    print "Eigenvalue range should be in -1 to 1:", e_values[0], e_values[-1]
    print "Max eigenvalue should be equal to 1:", e_values[-1]
    print e_values[-1] - 1.0

    # test
    if __debug__:
        # this tests max eigenvalue is equal to 1, within reasonable
        # precision
        test = (e_values[-1] - 1.0) < 1e-5
        if not test:
            print "<cluster.py> ERROR: Max eigenvalue is not 1.0."
            raise AssertionError

    # keep top few eigenvectors
    power = numpy.cumsum(e_values[::-1]) / numpy.sum(e_values)
    print "power: ", power[number_of_eigenvectors]
    embed = e_vectors[:, -number_of_eigenvectors - 1: -1]

    embed = numpy.divide(embed.T, e_vectors[:, -1]).T

    # reverse order of embedding so highest eigenvalue
    # information is first
    embed = embed[:, ::-1]

    print(embed.shape)

    # project onto eigenvectors using similarities for Nystrom method
    if use_nystrom:
        print('Projecting onto selected eigenvectors')
        # (this calculates embedding of a new point based on a weighted
        # sum of the embedded points)
        # create embedding vectors using nystrom approximation to find
        # the approximate top eigenvectors of the matrix
        # L = D^(-1/2) (D - W) D^(-1/2) 
        # See the paper:
        # "Spectral Grouping Using the Nystrom Method"
        # Note they divide embedding by 1st eigenvector rather
        # than sqrt of row sum. (As above in this code)       
        # note basically all this does is adds in the extra measurements
        # by projecting them onto the original eigenvector basis.
        # A=UVU' => U = AUV^-1 => take new rows of extended A (B') and do
        # the same.  U' = [AUV^-1 ; B'UV^-1] = [U ; B'UV^-1]
        # So find embedding for the unsampled fibers.
        print(rectangular_affinity_matrix.T.shape)
        print(e_vectors[:, -number_of_eigenvectors - 1: -1].T.shape)
        print(e_values.shape)
        embed2 = numpy.dot(rectangular_affinity_matrix.T, numpy.divide(e_vectors[:, -number_of_eigenvectors - 1: -1].T, e_values).T)
        print(embed2.shape)
        print(e_vectors.shape)
        print(e_vectors[:, -1].shape)

        # normalize estimated eigenvectors

        # reverse order of embedding so highest eigenvalue
        # information is first
        embed2 = embed2[:, ::-1]
        
        print "EMBED 1 values"
        print(embed.shape)
        print(numpy.min(embed))
        print(numpy.max(embed))
        
        print(embed2.shape)
        print(numpy.min(embed2))
        print(numpy.max(embed2))
        print "ENd embed vals"

        # replace embed variable with whole estimated embedding
        embed = numpy.concatenate((embed, embed2))

 

    # k-means in embed space
    print 'K-means clustering'
    centroids, distortion = scipy.cluster.vq.kmeans(embed, number_of_clusters)
    cluster_idx, dist = scipy.cluster.vq.vq(embed, centroids)

    #print 'View results'
    #plt.figure()
    #not_used = plt.hist(cluster_idx, 200)
    #plt.savefig('cluster_hist.pdf')

    # visualize embedding coordinates as RGB
    color = _embed_to_rgb(embed)

    # set up polydata with clustering output info.
    # for now modify input polydata by adding two arrays
    output_polydata = input_polydata
    output_polydata = \
        _format_output_polydata(output_polydata, cluster_idx, color, embed)

    # to view it, do this
    #ren = render.render(output_polydata)
    #render.save_views(ren)

    return output_polydata, cluster_idx, color, embed, distortion


def view_cluster_number(input_polydata, cluster_number, cluster_indices=None):

    """ Pop up a render window showing the selected cluster.


    Uses cluster_indices to choose corresponding cells in the
    polydata. If no argument cluster_indices is provided, then uses
    the cell data array named ClusterNumber.  One of these inputs must
    be present.

    """

    if cluster_indices == None:
        cluster_indices_vtk = \
            input_polydata.GetCellData().GetArray('ClusterNumber')
        cluster_indices = numpy.zeros(cluster_indices_vtk.GetNumberOfTuples())
        for fidx in range(0, cluster_indices_vtk.GetNumberOfTuples()):
            cluster_indices[fidx] = cluster_indices_vtk.GetTuple(fidx)[0]

    fiber_mask = cluster_indices == cluster_number
    view_polydata = filter.mask(input_polydata, fiber_mask)
    ren = render.render(view_polydata)

    return ren

def _rectangular_distance_matrix(input_polydata_n, input_polydata_m, threshold,
                              number_of_jobs=3):

    """ Internal convenience function available to clustering
    routines.

    Computes distance matrix (nxm) for all n+m fibers in input
    polydata. each fiber in input_polydata_n is compared to each fiber
    in input_polydata_m.


    """

    print 'Convert to array'
    fiber_array_n = fibers.FiberArray()
    fiber_array_n.convert_from_polydata(input_polydata_n)
    fiber_array_m = fibers.FiberArray()
    fiber_array_m.convert_from_polydata(input_polydata_m)

    # pairwise distance matrix
    all_fibers_n = range(0, fiber_array_n.number_of_fibers)

    print 'Compute fiber distances'
    distances = Parallel(n_jobs=number_of_jobs,
                         verbose=1)(
        delayed(similarity.fiber_distance)(
            fiber_array_n.get_fiber(lidx),
            fiber_array_m,
            threshold)
        for lidx in all_fibers_n)

    distances = numpy.array(distances).T

    # remove outliers if desired????

    print "**"
    print distances.shape
    print "**..."

    return distances

def _rectangular_similarity_matrix(input_polydata_n, input_polydata_m, threshold, sigma,
                                number_of_jobs=3):

    """ Internal convenience function available to clustering
    routines.

    Computes similarity matrix (nxn) for all n fibers in input
    polydata.  Calls function _pairwise_distance_matrix first.

    For distance computation and conversion to similarity uses
    parameters threshold and sigma.

    """

    distances = _rectangular_distance_matrix(input_polydata_n, input_polydata_m, threshold,
                                             number_of_jobs)

    print 'Convert to similarity'
    # similarity matrix
    sigmasq = sigma * sigma

    similarity_matrix = similarity.distance_to_similarity(distances, sigmasq)

    print "**"
    print similarity_matrix.shape
    print "**..."

    return similarity_matrix

def _pairwise_distance_matrix(input_polydata, threshold,
                              number_of_jobs=3):

    """ Internal convenience function available to clustering
    routines.

    Computes distance matrix (nxn) for all n fibers in input
    polydata.

    """

    print 'Convert to array'
    fiber_array = fibers.FiberArray()
    fiber_array.convert_from_polydata(input_polydata)

    # pairwise distance matrix
    all_fibers = range(0, fiber_array.number_of_fibers)

    print 'Compute fiber distances'
    distances = Parallel(n_jobs=number_of_jobs,
                         verbose=0)(
        delayed(similarity.fiber_distance)(
            fiber_array.get_fiber(lidx),
            fiber_array,
            threshold)
        for lidx in all_fibers)

    distances = numpy.array(distances)

    # remove outliers if desired????

    return distances

def _pairwise_similarity_matrix(input_polydata, threshold, sigma,
                                number_of_jobs=3):

    """ Internal convenience function available to clustering
    routines.

    Computes similarity matrix (nxn) for all n fibers in input
    polydata.  Calls function _pairwise_distance_matrix first.

    For distance computation and conversion to similarity uses
    parameters threshold and sigma.

    """

    distances = _pairwise_distance_matrix(input_polydata, threshold,
                                          number_of_jobs)

    print 'Convert to similarity'
    # similarity matrix
    sigmasq = sigma * sigma

    similarity_matrix = similarity.distance_to_similarity(distances, sigmasq)

    # sanity check that on-diagonal elements are all 1
    print "This should be 1.0: ", numpy.min(numpy.diag(similarity_matrix))
    print  numpy.min(numpy.diag(similarity_matrix)) == 1.0
    # test
    if __debug__:
        # this tests that on-diagonal elements are all 1
        test = numpy.min(numpy.diag(similarity_matrix)) == 1.0
        if not test:
            print "<cluster.py> ERROR: On-diagonal elements are not all 1.0."
            raise AssertionError

    return similarity_matrix


def _format_output_polydata(output_polydata, cluster_idx, color, embed):
    """ Output polydata with embedding colors, cluster numbers, and
    embedding coordinates.

    Cell data array names are:
    EmbeddingColor
    ClusterNumber
    EmbeddingCoordinate

    """

    embed_colors = vtk.vtkUnsignedCharArray()
    embed_colors.SetNumberOfComponents(3)
    embed_colors.SetName('EmbeddingColor')
    cluster_colors = vtk.vtkIntArray()
    cluster_colors.SetName('ClusterNumber')
    embed_data = vtk.vtkFloatArray()
    embed_data.SetNumberOfComponents(embed.shape[1])
    embed_data.SetName('EmbeddingCoordinate')

    for lidx in range(0, output_polydata.GetNumberOfLines()):
        embed_colors.InsertNextTuple3(
            color[lidx, 0], color[lidx, 1], color[lidx, 2])
        cluster_colors.InsertNextTuple1(cluster_idx[lidx])
        embed_data.InsertNextTupleValue(embed[lidx, :])

    output_polydata.GetCellData().AddArray(embed_data)
    output_polydata.GetCellData().AddArray(cluster_colors)
    output_polydata.GetCellData().AddArray(embed_colors)

    output_polydata.GetCellData().SetActiveScalars('EmbeddingColor')

    return output_polydata


def _embed_to_rgb(embed):
    """

    Get colors from embedding (first 3 coordinates corresponding to
    max eigenvalues). Output minimum color is 0,0,0, max is
    255,255,255, (actually 180 due to brightness constraint) and all
    output values (brightnesses) are equal.

     """

    # first 3 components of embedding
    color = embed[:, 0:3]

    # normalize all colors to length 1
    color_len = numpy.sqrt(numpy.sum(numpy.power(color, 2), 1))
    color = numpy.divide(color.T, color_len).T

    # convert to RGB: use colors from 0 to 255 for unsigned char (no
    # LUT in vtk). Before color components ranged from -1 to +1, now
    # from 0 to 255
    color = 127.5 + (color * 127.5)

    # Avoid dark black fibers that the cluster color is hard to see.
    # Convert to hsv, and make all brightness values the same.  That
    # way shadows show geometry, not cluster number.
    hue = numpy.zeros(color.shape[0])
    sat = numpy.zeros(color.shape[0])
    val = numpy.zeros(color.shape[0])
    for c_idx in range(0, color.shape[0]):
        hue[c_idx], sat[c_idx], val[c_idx] = \
            colorsys.rgb_to_hsv(
            color[c_idx, 0], color[c_idx, 1], color[c_idx, 2])

    val2 = numpy.ones(val.shape) * 180
    for c_idx in range(0, color.shape[0]):
        color[c_idx, 0], color[c_idx, 1], color[c_idx, 2] = \
            colorsys.hsv_to_rgb(hue[c_idx], sat[c_idx], val2[c_idx])

    # To look at this transform of the embed data (sanity check)
    #plt.figure()
    #plt.plot(embed[:, 0], color[:, 0], 'r.')
    #plt.plot(embed[:, 1], color[:, 1], 'g.')
    #plt.plot(embed[:, 2], color[:, 2], 'b.')

    return color

""" cluster.py

implementations of fiber clustering

"""
import os
import pickle

import vtk
import numpy
import colorsys
try:
    import scipy.cluster.vq
    import scipy.cluster.hierarchy
    USE_SCIPY = 1
except ImportError:
    USE_SCIPY = 0
    print "<cluster.py> Failed to import scipy.cluster, cannot cluster."
    print "<cluster.py> Please install scipy for this functionality."
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
import io

class ClusterAtlas:
    """Variables necessary to label a new subject from a spectral cluster atlas."""

    def __init__(self):
        self.pinv_A = None
        self.e_val = None
        self.e_vec = None
        self.row_sum_matrix = None
        self.row_sum_1 = None
        self.e_vec_norm = None
        self.nystrom_polydata = None
        self.threshold = None
        self.number_of_eigenvectors = None
        self.centroids = None
        self.sigma = None

    def save(self, directory, atlas_name):
        # temporarily remove polydata object to enable pickling
        polydata_tmp = self.nystrom_polydata
        self.nystrom_polydata = None
        fname = os.path.join(directory,atlas_name)
        pickle.dump(self,open(fname+'.p','wb'))
        # save the polydata
        io.write_polydata(polydata_tmp, fname+'.vtp')
        # replace the polydata object
        self.nystrom_polydata = polydata_tmp

    def load(self, directory, atlas_name):
        fname = os.path.join(directory,atlas_name)
        atlas = pickle.load(open(fname+'.p','rb'))
        atlas.nystrom_polydata = io.read_polydata(fname+'.vtp')
        return(atlas)
    
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
             number_of_jobs=3, use_nystrom=False, nystrom_mask = None,
             landmarks=None):

    """ Spectral clustering based on pairwise fiber affinity matrix.

    As in O'Donnell and Westin TMI 2007.

    Differences from that implementation: fiber distance is defined
    using fixed-length fiber parameterization.

    """

    # test pd has lines first
    number_fibers = input_polydata.GetNumberOfLines()
    print "<cluster.py> Starting spectral clustering."
    print "<cluster.py> Number of input fibers:", number_fibers
    print "<cluster.py> Number of clusters:", number_of_clusters

    if number_fibers == 0:
        print "<cluster.py> ERROR: Cannot cluster polydata with 0 fibers."
        return

    atlas = ClusterAtlas() 

    # 1) Compute fiber similarities.
    # Nystrom version of the code uses a sample of the data.
    if use_nystrom:
        # make sure it's an array for logic operations
        nystrom_mask = numpy.array(nystrom_mask)
        # make sure it's boolean or 0 and 1
        test = numpy.max(nystrom_mask) == 1.0
        if not test:
            print "<cluster.py> ERROR: Nystrom data mask is may not be Boolean. Max value is not 1.0/True."
            raise AssertionError
        # make sure it's large enough
        test = sum(nystrom_mask) >= 100
        if not test:
            print "<cluster.py> ERROR: Nystrom data mask is smaller than 100."
            raise AssertionError
        # make sure its size matches the polydata input
        test = len(nystrom_mask) == number_fibers
        if not test:
            print "<cluster.py> ERROR: Nystrom data mask size does not match polydata number of lines."
            raise AssertionError
        
        # Separate the Nystrom sample and the rest of the data.
        polydata_m = filter.mask(input_polydata, nystrom_mask)
        atlas.nystrom_polydata = polydata_m
        atlas.threshold = threshold
        polydata_n = filter.mask(input_polydata, nystrom_mask == False)
        sz = polydata_m.GetNumberOfLines()
        print '<cluster.py> Using Nystrom approximation. Subset size:',  sz, '/', number_fibers
        # Determine ordering to get embedding to correspond to original input data.
        reorder_embedding = numpy.concatenate((numpy.where(nystrom_mask)[0], numpy.where(~nystrom_mask)[0]))
        if landmarks is not None:
            landmarks_m = landmarks[nystrom_mask,:,:]
            landmarks_n = landmarks[~nystrom_mask,:,:]
        else:
            landmarks_m = landmarks_n = None

        # Calculate fiber similarities
        A = \
            _pairwise_similarity_matrix(polydata_m, threshold,
                                        sigma, number_of_jobs, landmarks_m)
        B = \
            _rectangular_similarity_matrix(polydata_n, polydata_m, threshold,
                                           sigma, number_of_jobs, landmarks_n, landmarks_m)
        atlas.sigma = sigma
    else:
        # Calculate all fiber similarities
        A = \
            _pairwise_similarity_matrix(input_polydata, threshold,
                                    sigma, number_of_jobs, landmarks)

    # 2) Do Normalized Cuts transform of similarity matrix.
    # See the paper: "Spectral Grouping Using the Nystrom Method"
    # (D^-1/2 W D^-1/2) V = V Lambda
    if use_nystrom:
        # calculate the sum of the rows we know from the full matrix
        atlas.row_sum_1 = numpy.sum(A, axis=0) + numpy.sum(B.T, axis=0)        
        # approximate the sum of the columns
        # weights for approximating row sum use the known columns of B'
        # pseudo inverse of A is needed for atlas (?)
        atlas.pinv_A = numpy.linalg.pinv(A)  

        e_val, e_vec = numpy.linalg.eigh(atlas.pinv_A)
        print "test of non-normalized A pseudoinverse Eigenvalue range:", e_val[0], e_val[-1]  

        # matlab was: atlas.approxRowSumMatrix = sum(B',1)*atlas.pseudoInverseA;
        # this matrix is needed for atlas.
        atlas.row_sum_matrix = numpy.dot(numpy.sum(B.T, axis=0), atlas.pinv_A)
        # row sum estimate for current B part of the matrix
        row_sum_2 = numpy.sum(B, axis=0) + \
            numpy.dot(atlas.row_sum_matrix, B)
        print "row sum check", numpy.min(atlas.row_sum_1), \
            numpy.max(atlas.row_sum_1), numpy.min(row_sum_2), \
            numpy.max(row_sum_2)
        print "check 2:", numpy.min(numpy.sum(B.T, axis=0)), numpy.min(atlas.row_sum_matrix)

        # normalized cuts normalization
        dhat = numpy.sqrt(numpy.divide(1, numpy.concatenate((atlas.row_sum_1, row_sum_2))))

        A = \
            numpy.multiply(A, numpy.outer(dhat[0:sz], dhat[0:sz].T))
        B = \
            numpy.multiply(B, numpy.outer(dhat[0:sz], dhat[sz:].T))
    else:
        # normalized cuts normalization using row (same as column) sums
        row_sum = numpy.sum(A, axis=0)
        dhat = numpy.divide(1, numpy.sqrt(row_sum))
        A = \
            numpy.multiply(A, numpy.outer(dhat, dhat.T))

    # 3) Compute eigenvectors for use in spectral embedding
    print '<cluster.py> Calculating eigenvectors of similarity matrix A...'
    atlas.e_val, atlas.e_vec = numpy.linalg.eigh(A)
    print '<cluster.py> Done calculating eigenvectors.'
    print "<cluster.py> Eigenvalue range:", atlas.e_val[0], atlas.e_val[-1]    
    # Check how well our chosen number of eigenvectors models the data
    power = numpy.cumsum(atlas.e_val[::-1]) / numpy.sum(atlas.e_val)
    print "<cluster.py> Power from chosen number of eigenvectors (", number_of_eigenvectors, ')', power[number_of_eigenvectors]
    print '<cluster.py> Top eigenvalues:', atlas.e_val[::-1][1:number_of_eigenvectors]

    # 4) Compute embedding using eigenvectors
    print('<cluster.py> Compute embedding using eigenvectors.')
    if use_nystrom:
        # Create embedding vectors using nystrom approximation to find
        # the approximate top eigenvectors of the matrix
        # L = D^(-1/2) (D - W) D^(-1/2) 
        # See the paper:
        # "Spectral Grouping Using the Nystrom Method"
        # Basically all this does is adds in the extra measurements
        # by projecting them onto the original eigenvector basis.
        # A=UVU' => U = AUV^-1 => take new rows of extended A (B') and do
        # the same.  U' = [AUV^-1 ; B'UV^-1] = [U ; B'UV^-1]
        # Note they divide embedding by 1st eigenvector rather
        # than sqrt of row sum, as in this code (below).

        # matlab was: % project onto eigenvectors of A:
        # % v' = [v ; B'*v*d^-1 
        # V = [atlas.eigenvectA; B'*atlas.eigenvectA*(diag(1./diag(atlas.eigenvalA)))];
        V = numpy.concatenate((atlas.e_vec, \
                                  numpy.dot(numpy.dot(B.T, atlas.e_vec), \
                                                numpy.diag(numpy.divide(1.0, atlas.e_val)))))

        # normalize estimated eigenvectors to have length of one
        # matlab was:
        # atlas.eigenvectorLengthToNormalize=sqrt(sum(V.*V));
        # V=V./repmat(atlas.eigenvectorLengthToNormalize,length(V),1);
        atlas.e_vec_norm = numpy.sum(numpy.multiply(V, V),0)
        V = numpy.divide(V, atlas.e_vec_norm)

        # Normalize each embedding vector by first eigenvector. Matlab code was:
        # for i = 2:embedLength+1
        #    embedding(:,i-1) = V(:,i)./V(:,1);
        # end
        # This eigenvector corresponds to an eigenvalue of 1, since row sums are 1.
        # The other option from the literature was to use this:
        # embedding_i,j = V_i+i,j./sqrt(D_j,j)
        embed = numpy.zeros((number_fibers, number_of_eigenvectors))
        for i in range(0, number_of_eigenvectors):
            embed[reorder_embedding,i] = numpy.divide(V[:,-(i+2)], V[:,-1])
    else:
        embed = atlas.e_vec[:, -number_of_eigenvectors - 2: -2]
        embed = numpy.divide(embed.T, atlas.e_vec[:, -1]).T
        # reverse order of embedding so highest eigenvalue
        # information is first
        embed = embed[:, ::-1]

    atlas.number_of_eigenvectors = number_of_eigenvectors

    # 5) Find clusters using k-means in embedding space.
    print '<cluster.py> K-means clustering in embedding space.'
    atlas.centroids, distortion = scipy.cluster.vq.kmeans(embed, number_of_clusters)
    cluster_idx, dist = scipy.cluster.vq.vq(embed, atlas.centroids)

    # 6) Output results.
    print '<cluster.py> Done spectral clustering, returning results.'
    # visualize embedding coordinates as RGB
    embed2 = embed
    embed2[numpy.isnan(embed)] = 0.0
    color = _embed_to_rgb(embed2)
    # set up polydata with clustering output info.
    # for now modify input polydata by adding two arrays
    output_polydata = input_polydata
    output_polydata = \
        _format_output_polydata(output_polydata, cluster_idx, color, embed)

    return output_polydata, cluster_idx, color, embed, distortion, atlas


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

def spectral_atlas_label(input_polydata, atlas, number_of_jobs=2):
    """ test

    test


    """
    number_fibers = input_polydata.GetNumberOfLines()
    sz = atlas.nystrom_polydata.GetNumberOfLines()

    # 1) Compute fiber similarities.
    B = \
        _rectangular_similarity_matrix(input_polydata, atlas.nystrom_polydata, 
                                       atlas.threshold, atlas.sigma, number_of_jobs)

    # 2) Do Normalized Cuts transform of similarity matrix.
    # row sum estimate for current B part of the matrix
    row_sum_2 = numpy.sum(B, axis=0) + \
        numpy.dot(atlas.row_sum_matrix, B)
    # normalized cuts normalization
    dhat = numpy.sqrt(numpy.divide(1, numpy.concatenate((atlas.row_sum_1, row_sum_2))))
    B = \
        numpy.multiply(B, numpy.outer(dhat[0:sz], dhat[sz:].T))

    # 3) Compute eigenvectors for use in spectral embedding
    # <done already in atlas creation>

    # 4) Compute embedding using eigenvectors
    V = numpy.dot(numpy.dot(B.T, atlas.e_vec), \
                      numpy.diag(numpy.divide(1.0, atlas.e_val)))
    V = numpy.divide(V, atlas.e_vec_norm)
    embed = numpy.zeros((number_fibers, atlas.number_of_eigenvectors))
    for i in range(0, atlas.number_of_eigenvectors):
        embed[:,i] = numpy.divide(V[:,-(i+2)], V[:,-1])

    # 5) Find clusters using k-means in embedding space.
    # Actually: LABEL using centroids from atlas
    print '<cluster.py>Cluster labeling in embedding space.'
    cluster_idx, dist = scipy.cluster.vq.vq(embed, atlas.centroids)

    # 6) Output results.
    print '<cluster.py> Done spectral clustering, returning results.'
    # visualize embedding coordinates as RGB
    color = _embed_to_rgb(embed)
    # set up polydata with clustering output info.
    # for now modify input polydata by adding two arrays
    output_polydata = input_polydata
    output_polydata = \
        _format_output_polydata(output_polydata, cluster_idx, color, embed)

    return output_polydata, cluster_idx, color, embed

def _rectangular_distance_matrix(input_polydata_n, input_polydata_m, threshold,
                              number_of_jobs=3, landmarks_n=None, landmarks_m=None):

    """ Internal convenience function available to clustering
    routines.

    Computes distance matrix (nxm) for all n+m fibers in input
    polydata. each fiber in input_polydata_n is compared to each fiber
    in input_polydata_m.


    """

    fiber_array_n = fibers.FiberArray()
    fiber_array_n.convert_from_polydata(input_polydata_n)
    fiber_array_m = fibers.FiberArray()
    fiber_array_m.convert_from_polydata(input_polydata_m)

    if landmarks_n is None:
        landmarks_n = numpy.zeros((fiber_array_n.number_of_fibers,3))

    # pairwise distance matrix
    all_fibers_n = range(0, fiber_array_n.number_of_fibers)

    distances = Parallel(n_jobs=number_of_jobs,
                         verbose=1)(
        delayed(similarity.fiber_distance)(
            fiber_array_n.get_fiber(lidx),
            fiber_array_m,
            threshold, distance_method='Hausdorff',
            fiber_landmarks=landmarks_n[lidx,:], 
            landmarks=landmarks_m)
        for lidx in all_fibers_n)

    distances = numpy.array(distances).T

    return distances

def _rectangular_similarity_matrix(input_polydata_n, input_polydata_m, threshold, sigma,
                                number_of_jobs=3, landmarks_n=None, landmarks_m=None):

    """ Internal convenience function available to clustering
    routines.

    Computes similarity matrix (nxn) for all n fibers in input
    polydata.  Calls function _pairwise_distance_matrix first.

    For distance computation and conversion to similarity uses
    parameters threshold and sigma.

    """

    distances = _rectangular_distance_matrix(input_polydata_n, input_polydata_m, threshold,
                                             number_of_jobs, landmarks_n, landmarks_m)

    # similarity matrix
    sigmasq = sigma * sigma

    similarity_matrix = similarity.distance_to_similarity(distances, sigmasq)

    return similarity_matrix

def _pairwise_distance_matrix(input_polydata, threshold,
                              number_of_jobs=3, landmarks=None):

    """ Internal convenience function available to clustering
    routines.

    Computes distance matrix (nxn) for all n fibers in input
    polydata.

    """

    fiber_array = fibers.FiberArray()
    fiber_array.convert_from_polydata(input_polydata)

    # pairwise distance matrix
    all_fibers = range(0, fiber_array.number_of_fibers)

    if landmarks is None:
        landmarks2 = numpy.zeros((fiber_array.number_of_fibers,3))
    else:
        landmarks2 = landmarks

    distances = Parallel(n_jobs=number_of_jobs,
                         verbose=0)(
        delayed(similarity.fiber_distance)(
            fiber_array.get_fiber(lidx),
            fiber_array,
            threshold, distance_method='Hausdorff', 
            fiber_landmarks=landmarks2[lidx,:], 
            landmarks=landmarks)
        for lidx in all_fibers)

    distances = numpy.array(distances)

    # remove outliers if desired????

    return distances

def _pairwise_similarity_matrix(input_polydata, threshold, sigma,
                                number_of_jobs=3, landmarks=None):

    """ Internal convenience function available to clustering
    routines.

    Computes similarity matrix (nxn) for all n fibers in input
    polydata.  Calls function _pairwise_distance_matrix first.

    For distance computation and conversion to similarity uses
    parameters threshold and sigma.

    """

    distances = _pairwise_distance_matrix(input_polydata, threshold,
                                          number_of_jobs, landmarks)

    # similarity matrix
    sigmasq = sigma * sigma

    similarity_matrix = similarity.distance_to_similarity(distances, sigmasq)

    # sanity check that on-diagonal elements are all 1
    #print "This should be 1.0: ", numpy.min(numpy.diag(similarity_matrix))
    #print  numpy.min(numpy.diag(similarity_matrix)) == 1.0
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
        cluster_colors.InsertNextTuple1(int(cluster_idx[lidx]))
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
    print numpy.min(color_len)
    color_len += 0.0001
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

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
import mrml

from pprint import pprint

HAVE_PLT = 1
try:
    import matplotlib
    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except:
    print "<wm_cluster.py> Error importing matplotlib.pyplot package, can't plot quality control data.\n"
    HAVE_PLT = 0
    
# This did not work better. Leave here for future testing if of interest
if 0:
    try:    
        from sklearn.cluster import AffinityPropagation
        from sklearn import metrics
    except ImportError:
        SKLEARN = 0
        print "<cluster.py> Failed to import sklearn, cannot use affinity propagation."
        print "<cluster.py> Please install sklearn for this functionality."
    
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
        self.bilateral = None
        self.distance_method = None
        self.version = 1.0
        self.polydata_filename = None
        # Outlier removal
        self.cluster_outlier_std_threshold = None
        self.cluster_mean_similarity = None
        self.cluster_std_similarity = None
        self.brain_mean_similarity = None
        self.brain_std_similarity = None
        
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

    def load(self, directory, atlas_name, verbose=False):
        if not os.path.isdir(directory):
            print "Error: Atlas directory", directory, "does not exist or is not a directory."
            raise "<cluster.py> I/O error"
        
        fname_base = os.path.join(directory,atlas_name)
        fname_atlas = fname_base+'.p'
        fname_polydata = fname_base+'.vtp'
        
        if not os.path.exists(fname_atlas):
            print "Error: Atlas file", fname_atlas, "does not exist."
            raise "<cluster.py> I/O error"
        if not os.path.exists(fname_polydata):
            print "Error: Atlas file", fname_polydata, "does not exist."
            raise "<cluster.py> I/O error"

        atlas = pickle.load(open(fname_atlas,'rb'))
        atlas.nystrom_polydata = io.read_polydata(fname_polydata)
        print "<cluster.py> Loaded atlas", atlas_name, "from", directory, "."
        print "<cluster.py> Atlas Nystrom polydata sample:", atlas.nystrom_polydata.GetNumberOfLines(), \
            "\nAtlas size:", atlas.pinv_A.shape, "\nAtlas number of eigenvectors:", atlas.number_of_eigenvectors

        atlas.polydata_filename = fname_polydata

        # check for any version issues
        if not hasattr(atlas, 'version'):    
            atlas.version = '0.0'
        # this is an atlas from pre 1.0 version code. Set to default values.
        if not hasattr(atlas, 'distance_method'):    
            atlas.distance_method = 'Hausdorff'
        if not hasattr(atlas, 'bilateral'):    
            atlas.bilateral = 'False'

        if verbose:
            pprint (vars(atlas))

        return(atlas)

def load_atlas(directory, atlas_name):
    """ Convenience function to load a spectral cluster atlas """
    tmp_atlas = ClusterAtlas()
    return tmp_atlas.load(directory, atlas_name)
        
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


def nearPSD(A,epsilon=0):
   n = A.shape[0]
   eigval, eigvec = numpy.linalg.eig(A)
   val = numpy.matrix(numpy.maximum(eigval,epsilon))
   vec = numpy.matrix(eigvec)
   T = 1/(numpy.multiply(vec,vec) * val.T)
   T = numpy.matrix(numpy.sqrt(numpy.diag(numpy.array(T).reshape((n)) )))
   B = T * vec * numpy.diag(numpy.array(numpy.sqrt(val)).reshape((n)))
   out = B*B.T
   return(numpy.asarray(out))
   
def spectral(input_polydata, number_of_clusters=200,
             number_of_eigenvectors=20, sigma=60, threshold=0.0,
             number_of_jobs=3, use_nystrom=False, nystrom_mask=None,
             landmarks=None, distance_method='Mean', normalized_cuts=True,
             outlier_std_threshold = 2.0,
             pos_def_approx=True,
             bilateral=False):

    """ Spectral clustering based on pairwise fiber affinity matrix.

    As in O'Donnell and Westin TMI 2007.

    Differences from that implementation: fiber point correspondences are defined
    using fixed-length fiber parameterization (instead of closest point).

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

    # Store all parameters to this function. They must be identical later to label new data.
    # Below, calculated values will also be stored in the atlas.
    atlas.number_of_eigenvectors = number_of_eigenvectors
    atlas.sigma = sigma
    atlas.threshold = threshold
    atlas.use_nystrom = use_nystrom
    atlas.landmarks = landmarks
    atlas.distance_method = distance_method
    atlas.bilateral = bilateral

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
        polydata_m = filter.mask(input_polydata, nystrom_mask, verbose=False)
        atlas.nystrom_polydata = polydata_m
        polydata_n = filter.mask(input_polydata, nystrom_mask == False, verbose=False)
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
                                        sigma, number_of_jobs, landmarks_m, distance_method, bilateral)
        B = \
            _rectangular_similarity_matrix(polydata_n, polydata_m, threshold,
                                           sigma, number_of_jobs, landmarks_n, landmarks_m, distance_method, bilateral)

        # sanity check
        print "<cluster.py> Range of values in A:", numpy.min(A), numpy.max(A)
        print "<cluster.py> Range of values in B:", numpy.min(B), numpy.max(B)
        
    else:
        # Calculate all fiber similarities
        A = \
            _pairwise_similarity_matrix(input_polydata, threshold,
                                    sigma, number_of_jobs, landmarks, distance_method, bilateral)

        atlas.nystrom_polydata = input_polydata
        # sanity check
        print "<cluster.py> Range of values in A:", numpy.min(A), numpy.max(A)
        
    testval = numpy.max(A-A.T) 
    if not testval == 0.0:
        if testval > 1e-10:
            print "<cluster.py> ERROR: A matrix is not symmetric."
            raise AssertionError
        else:
            print "<cluster.py> Maximum of A - A^T:", testval
        # Ensure that A is symmetric
        A = numpy.divide(A+A.T, 2.0)
        
    testval = numpy.min(A)
    if not testval > 0.0:
        print "<cluster.py> ERROR: A matrix is not positive."
        print "<cluster.py> Minimum value in A: ", testval
        if testval < 0.0:
            raise AssertionError

    # Outliers will have low measured (or estimated) row sums. Detect outliers in A:
    # to turn off for testing: outlier_std_threshold = numpy.inf
    row_sum_A_initial = numpy.sum(A, axis=0) + numpy.sum(B.T, axis=0)
    print "<cluster.py> Initial similarity (row) sum A:", numpy.mean(row_sum_A_initial), numpy.std(row_sum_A_initial), numpy.min(row_sum_A_initial)
    atlas.outlier_std_threshold = outlier_std_threshold
    atlas.row_sum_threshold_for_rejection = numpy.mean(row_sum_A_initial) - outlier_std_threshold*numpy.std(row_sum_A_initial)
    bad_idx = numpy.nonzero(row_sum_A_initial < atlas.row_sum_threshold_for_rejection)[0]
    reject_A = bad_idx
    print "<cluster.py> Rejecting n=", len(bad_idx), "/", sz, "fibers >", outlier_std_threshold, "standard deviations below the mean total fiber similarity"

    A = numpy.delete(A,reject_A,0)
    A = numpy.delete(A,reject_A,1)
    #print A.shape, B.shape
    B = numpy.delete(B,reject_A,0)
    #print A.shape, B.shape, reorder_embedding.shape
                    
    # Ensure that A is positive definite.
    if pos_def_approx:
        e_val, e_vec = numpy.linalg.eigh(A)
        print "<cluster.py> Eigenvalue range of A:", e_val[0], e_val[-1]
        A2 = nearPSD(A)
        e_val, e_vec = numpy.linalg.eigh(A2)
        print "<cluster.py> Eigenvalue range of nearest PSD matrix to A:", e_val[0], e_val[-1]  
        testval = numpy.max(A-A2) 
        if not testval == 0.0:
            print "<cluster.py> A matrix differs by PSD matrix by maximum of:", testval
            if testval > 0.25:
                print "<cluster.py> ERROR: A matrix changed by more than 0.25."
                raise AssertionError
        A = A2
        
    # 2) Do Normalized Cuts transform of similarity matrix.
    # See the paper: "Spectral Grouping Using the Nystrom Method"
    # (D^-1/2 W D^-1/2) V = V Lambda
    if normalized_cuts:
        if use_nystrom:
            # Form of entire affinity matrix: 
            # A   B
            # B^T   C
            # C is not computed.
            # Calculate the sum of the partial rows we've computed:
            atlas.row_sum_1 = numpy.sum(A, axis=0) + numpy.sum(B.T, axis=0)  
            #print "<cluster.py> A size:", A.shape
            #print "<cluster.py> B size:", B.shape
            #print "<cluster.py> A-B matrix row sums range (should be > 0):", numpy.min(atlas.row_sum_1), numpy.max(atlas.row_sum_1)
            
            # Approximate the sum of the rest of the data (including C)
            # These are weighted sums of the columns we did compute
            # where the weight depends on how similar that fiber 
            # was to each path in A.  This uses the dual basis
            # of the columns in A.
            # Approximate the inverse of A for dual basis
            #print "<cluster.py> Using numpy linalg pinv A"
            atlas.pinv_A = numpy.linalg.pinv(A)

            #e_val, e_vec = numpy.linalg.eigh(atlas.pinv_A)
            #print "<cluster.py> test of non-normalized A pseudoinverse Eigenvalue range:", e_val[0], e_val[-1]  

            # row sum formula:
            # dhat = [a_r + b_r; b_c + B^T*A-1*b_r]
            # this matrix is A^-1 * b_r, where b_r are the row sums of B
            # matlab was: atlas.approxRowSumMatrix = sum(B',1)*atlas.pseudoInverseA;
            atlas.row_sum_matrix = numpy.dot(numpy.sum(B.T, axis=0), atlas.pinv_A)
            #test = numpy.sum(B.T, axis=0)
            #print "<cluster.py> B column sums range (should be > 0):", numpy.min(test), numpy.max(test)
            print "<cluster.py> Range of row sum weights:", numpy.min(atlas.row_sum_matrix), numpy.max(atlas.row_sum_matrix)
            #print "<cluster.py> First 10 entries in weight matrix:", atlas.row_sum_matrix[0:10]
            #test = numpy.dot(atlas.row_sum_matrix, B)
            #print "<cluster.py> Test partial sum estimation for B:", numpy.min(test), numpy.max(test)
            #del test
            
            # row sum estimate for current B part of the matrix
            row_sum_2 = numpy.sum(B, axis=0) + \
                numpy.dot(atlas.row_sum_matrix, B)
            print "<cluster.py> Row sum check (min/max, should be > 0) A:", numpy.min(atlas.row_sum_1), numpy.median(atlas.row_sum_1), numpy.max(atlas.row_sum_1),  "B:", numpy.min(row_sum_2), numpy.median(row_sum_2), numpy.max(row_sum_2)

            # reject outliers in B
            bad_idx = numpy.nonzero(row_sum_2 < atlas.row_sum_threshold_for_rejection)[0]
            reject_B = bad_idx
            print "<cluster.py> Rejecting n=", len(bad_idx), "/", B.shape[1], "fibers >", outlier_std_threshold, "standard deviations below the mean total fiber similarity"
            row_sum_2 = numpy.delete(row_sum_2,reject_B)
            B = numpy.delete(B,reject_B,1)

            print "<cluster.py> After outlier rejection A:", A.shape, "B:", B.shape
            print "<cluster.py> Row sum check (min/max, should be > 0) A:", numpy.min(atlas.row_sum_1), numpy.median(atlas.row_sum_1), numpy.max(atlas.row_sum_1),  "B:", numpy.min(row_sum_2), numpy.median(row_sum_2), numpy.max(row_sum_2)

            # Separate the Nystrom sample and the rest of the data after removing outliers
            nystrom_mask_2 = nystrom_mask
            midxA = numpy.nonzero(nystrom_mask_2)[0]
            nystrom_mask_2[midxA[reject_A]] = False
            not_nystrom_mask = nystrom_mask == False
            not_nystrom_mask[midxA[reject_A]] = False
            midxB = numpy.nonzero(not_nystrom_mask)[0]
            not_nystrom_mask[midxB[reject_B]] = False

            polydata_m = filter.mask(input_polydata, nystrom_mask_2, verbose=False)
            atlas.nystrom_polydata = polydata_m
            polydata_n = filter.mask(input_polydata, not_nystrom_mask, verbose=False)
            output_polydata = filter.mask(input_polydata, numpy.add(nystrom_mask_2, not_nystrom_mask),verbose=False)
            sz = polydata_m.GetNumberOfLines()
            number_fibers = output_polydata.GetNumberOfLines()
            print '<cluster.py> Using Nystrom approximation. Subset size (A):',  sz, '/', number_fibers, "B:", polydata_n.GetNumberOfLines()
            # Determine ordering to get embedding to correspond to original input data.
            # reject outliers from masks
            reject_idx = numpy.concatenate((midxA[reject_A],midxB[reject_B]))
            nystrom_mask_2 = numpy.delete(nystrom_mask_2,reject_idx)
            not_nystrom_mask = numpy.delete(not_nystrom_mask,reject_idx)
            #print "hi after mask:", reorder_embedding.shape, numpy.sum(nystrom_mask_2), numpy.sum(not_nystrom_mask)
            reorder_embedding = numpy.concatenate((numpy.where(nystrom_mask_2)[0], numpy.where(not_nystrom_mask)[0]))
            #print "hi after embed reorder calc:", reorder_embedding.shape, numpy.max(reorder_embedding), numpy.min(reorder_embedding)
           
            # in case of negative row sum estimation
            if any(row_sum_2<=0):
                print "<cluster.py> Warning: Consider increasing sigma or using the Mean distance. negative row sum approximations."
                print "Number of negative row sums:", numpy.count_nonzero(row_sum_2<=0)
                #row_sum_2[row_sum_2<0] = 0.1

            # save for testing
            column_sum = numpy.concatenate((numpy.sum(A, axis=1) , numpy.sum(B.T, axis=1)))

            # normalized cuts normalization
            row_sum = numpy.concatenate((atlas.row_sum_1, row_sum_2))
            dhat = numpy.sqrt(numpy.divide(1, row_sum))
            #dhat = numpy.sqrt(numpy.divide(1, numpy.concatenate((atlas.row_sum_1, row_sum_2))))

            A = \
                numpy.multiply(A, numpy.outer(dhat[0:sz], dhat[0:sz].T))
            B = \
                numpy.multiply(B, numpy.outer(dhat[0:sz], dhat[sz:].T))

        else:
            # normalized cuts normalization using row (same as column) sums
            row_sum = numpy.sum(A, axis=0)
            print "<cluster.py> A matrix row sums range (should be > 0):", numpy.min(row_sum), numpy.max(row_sum)
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
        #print "Hi 3:", embed.shape, number_fibers, reorder_embedding.shape, V.shape, A.shape, B.shape
        for i in range(0, number_of_eigenvectors):
            embed[reorder_embedding,i] = numpy.divide(V[:,-(i+2)], V[:,-1])
    else:
        embed = atlas.e_vec[:, -number_of_eigenvectors - 2: -2]
        embed = numpy.divide(embed.T, atlas.e_vec[:, -1]).T
        # reverse order of embedding so highest eigenvalue
        # information is first
        embed = embed[:, ::-1]


    # Default is always k-means. Other code is just left for testing. Did not improve results.
    #centroid_finder = 'AffinityPropagation'
    centroid_finder = 'K-means'
    
    # 5) Find clusters using k-means in embedding space.
    cluster_metric = None
    if centroid_finder == 'K-means':
        print '<cluster.py> K-means clustering in embedding space.'
        centroids, cluster_metric = scipy.cluster.vq.kmeans2(embed, number_of_clusters, minit='points')
        # sort centroids by first eigenvector order
        # centroid_order = numpy.argsort(centroids[:,0])
        # sort centroids according to colormap and save them in this order in atlas
        color = _embed_to_rgb(centroids)
        centroid_order = render.argsort_by_jet_lookup_table(color)
        atlas.centroids = centroids[centroid_order,:]
        cluster_idx, dist = scipy.cluster.vq.vq(embed, atlas.centroids)
        #print "<cluster.py> Distortion metric:", cluster_metric
        if 0:
            # This is extremely slow, but leave code here if ever wanted for testing
            cluster_metric = metrics.silhouette_score(embed, cluster_idx, metric='sqeuclidean')
            print("Silhouette Coefficient: %0.3f" % cluster_metric)
 
    else:
        print "ERROR: Unknown centroid finder", centroid_finder
        ## # This found fewer clusters than we need to represent the anatomy well
        ## # Leave code here in case wanted in future for more testing.
        ## print '<cluster.py> Affinity Propagation clustering in embedding space.'
        ## af = AffinityPropagation(preference=-50).fit(embed)
        ## cluster_centers_indices = af.cluster_centers_indices_
        ## labels = af.labels_
        ## n_clusters_ = len(cluster_centers_indices)
        ## print('Estimated number of clusters: %d' % n_clusters_)
        ## cluster_idx = labels
        ## for k in range(n_clusters_):
        ##     class_members = labels == k
        ##     atlas.centroids = embed[cluster_centers_indices[k]]
        ## # return metrics
        ## if 0:
        ##     # This is extremely slow, but leave code here if ever wanted for testing
        ##     cluster_metric = metrics.silhouette_score(embed, labels, metric='sqeuclidean')
        ##     print("Silhouette Coefficient: %0.3f" % cluster_metric)

    # 6) Output results.
    print '<cluster.py> Done spectral clustering, returning results.'
    # visualize embedding coordinates as RGB
    embed2 = embed
    #embed2[numpy.isnan(embed)] = 0.0
    color = _embed_to_rgb(embed2)
    # set up polydata with clustering output info.
    # for now modify input polydata by adding two arrays
    #output_polydata = input_polydata
    if use_nystrom:
        output_polydata = \
            _format_output_polydata(output_polydata, cluster_idx, color, embed, row_sum[reorder_embedding], column_sum[reorder_embedding])
    else:
        # row and column sums are the same. no need to reorder.
        output_polydata = \
            _format_output_polydata(output_polydata, cluster_idx, color, embed, row_sum, row_sum)
            
    return output_polydata, cluster_idx, color, embed, cluster_metric, atlas, reject_idx 


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
    """ Use an existing atlas to label a new polydata.

    Returns the cluster indices for all the fibers. output_polydata, cluster_numbers, color, embed = wma.cluster.spectral_atlas_label(input_data, atlas)

    """

    number_fibers = input_polydata.GetNumberOfLines()
    sz = atlas.nystrom_polydata.GetNumberOfLines()
    
    # 1) Compute fiber similarities.
    B = \
        _rectangular_similarity_matrix(input_polydata, atlas.nystrom_polydata, 
                                       atlas.threshold, atlas.sigma, number_of_jobs, distance_method=atlas.distance_method,
                                       bilateral=atlas.bilateral)

    # 2) Do Normalized Cuts transform of similarity matrix.
    # row sum estimate for current B part of the matrix
    row_sum_2 = numpy.sum(B, axis=0) + \
        numpy.dot(atlas.row_sum_matrix, B)

    # in case of negative row sum estimation (this should not ever happen)
    if any(row_sum_2<=0):
        print "<cluster.py> Warning: Consider increasing sigma or using the Mean distance. negative row sum approximations."
        print "Number of negative row sums:", numpy.count_nonzero(row_sum_2<=0)
        row_sum_2[row_sum_2<0] = 0.1

             
    # normalized cuts normalization
    row_sum = numpy.concatenate((atlas.row_sum_1, row_sum_2))
    dhat = numpy.sqrt(numpy.divide(1, row_sum))
    #dhat = numpy.sqrt(numpy.divide(1, numpy.concatenate((atlas.row_sum_1, row_sum_2))))
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
    print '<cluster.py> Cluster labeling in embedding space.'
    cluster_idx, dist = scipy.cluster.vq.vq(embed, atlas.centroids)

    # 6) Output results.
    print '<cluster.py> Done labeling clusters, returning results.'
    # visualize embedding coordinates as RGB
    color = _embed_to_rgb(embed)
    # set up polydata with clustering output info.
    # for now modify input polydata by adding two arrays
    output_polydata = input_polydata
    output_polydata = \
        _format_output_polydata(output_polydata, cluster_idx, color, embed, row_sum, numpy.sum(B.T, axis=1))

    return output_polydata, cluster_idx, color, embed

def _rectangular_distance_matrix(input_polydata_n, input_polydata_m, threshold,
                              number_of_jobs=3, landmarks_n=None, landmarks_m=None,
                              distance_method='Hausdorff', bilateral=False):

    """ Internal convenience function available to clustering
    routines.

    Computes distance matrix (nxm) for all n+m fibers in input
    polydata. each fiber in input_polydata_n is compared to each fiber
    in input_polydata_m.


    """

    if distance_method == 'Frechet':

        distances = similarity.rectangular_frechet_distances(input_polydata_n, input_polydata_m)
        distances = numpy.array(distances)

    else:
        
        fiber_array_n = fibers.FiberArray()
        fiber_array_n.convert_from_polydata(input_polydata_n, points_per_fiber=15)
        fiber_array_m = fibers.FiberArray()
        fiber_array_m.convert_from_polydata(input_polydata_m, points_per_fiber=15)

        if landmarks_n is None:
            landmarks_n = numpy.zeros((fiber_array_n.number_of_fibers,3))
    
        # pairwise distance matrix
        all_fibers_n = range(0, fiber_array_n.number_of_fibers)

        distances = Parallel(n_jobs=number_of_jobs,
                             verbose=0)(
            delayed(similarity.fiber_distance)(
                fiber_array_n.get_fiber(lidx),
                fiber_array_m,
                threshold, distance_method=distance_method,
                fiber_landmarks=landmarks_n[lidx,:], 
                landmarks=landmarks_m, bilateral=bilateral)
            for lidx in all_fibers_n)

        distances = numpy.array(distances).T

    return distances

def _rectangular_similarity_matrix(input_polydata_n, input_polydata_m, threshold, sigma,
                                number_of_jobs=3, landmarks_n=None, landmarks_m=None, distance_method='Hausdorff',
                                bilateral=False):

    """ Internal convenience function available to clustering
    routines.

    Computes similarity matrix (nxn) for all n fibers in input
    polydata.  Calls function _pairwise_distance_matrix first.

    For distance computation and conversion to similarity uses
    parameters threshold and sigma.

    """

    distances = _rectangular_distance_matrix(input_polydata_n, input_polydata_m, threshold,
                                             number_of_jobs, landmarks_n, landmarks_m, distance_method, bilateral=bilateral)

    if distance_method == 'StrictSimilarity':
        similarity_matrix = distances
    else:
        # similarity matrix
        sigmasq = sigma * sigma
        similarity_matrix = similarity.distance_to_similarity(distances, sigmasq)

    return similarity_matrix

def _pairwise_distance_matrix(input_polydata, threshold,
                              number_of_jobs=3, landmarks=None, distance_method='Hausdorff',
                              bilateral=False, sigmasq=6400):

    """ Internal convenience function available to clustering
    routines.

    Computes distance matrix (nxn) for all n fibers in input
    polydata.

    """

    if distance_method == 'Frechet':

        distances = similarity.pairwise_frechet_distances(input_polydata,input_polydata)
        distances = numpy.array(distances)

    else:

        fiber_array = fibers.FiberArray()
        fiber_array.convert_from_polydata(input_polydata, points_per_fiber=15)

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
                threshold, distance_method=distance_method, 
                fiber_landmarks=landmarks2[lidx,:], 
                landmarks=landmarks, bilateral=bilateral, sigmasq=sigmasq)
            for lidx in all_fibers)

        distances = numpy.array(distances)

        # remove outliers if desired????

    return distances

def _pairwise_similarity_matrix(input_polydata, threshold, sigma,
                                number_of_jobs=3, landmarks=None, distance_method='Hausdorff',
                                bilateral=False):

    """ Internal convenience function available to clustering
    routines.

    Computes similarity matrix (nxn) for all n fibers in input
    polydata.  Calls function _pairwise_distance_matrix first.

    For distance computation and conversion to similarity uses
    parameters threshold and sigma.

    """

    distances = _pairwise_distance_matrix(input_polydata, threshold,
                                          number_of_jobs, landmarks, distance_method, bilateral=bilateral)
    
    if distance_method == 'StrictSimilarity':
        similarity_matrix = distances
    else:
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
            print" Minimum on-diagonal value:", numpy.min(numpy.diag(similarity_matrix))
            print" Maximum on-diagonal value:", numpy.max(numpy.diag(similarity_matrix))
            print" Mean value:", numpy.mean(numpy.diag(similarity_matrix))
            raise AssertionError

    return similarity_matrix


def _format_output_polydata(output_polydata, cluster_idx, color, embed, estimated_row_sum, computed_column_sum):
    """ Output polydata with embedding colors, cluster numbers, and
    embedding coordinates.

    Cell data array names are:
    EmbeddingColor
    ClusterNumber
    EmbeddingCoordinate

    """
    #print computed_column_sum.shape
     
    embed_colors = vtk.vtkUnsignedCharArray()
    embed_colors.SetNumberOfComponents(3)
    embed_colors.SetName('EmbeddingColor')
    cluster_colors = vtk.vtkIntArray()
    cluster_colors.SetName('ClusterNumber')
    embed_data = vtk.vtkFloatArray()
    embed_data.SetNumberOfComponents(embed.shape[1])
    embed_data.SetName('EmbeddingCoordinate')
    row_sum_data = vtk.vtkFloatArray()
    row_sum_data.SetName('TotalFiberSimilarity')
    col_sum_data = vtk.vtkFloatArray()
    col_sum_data.SetName('MeasuredFiberSimilarity')
    
    for lidx in range(0, output_polydata.GetNumberOfLines()):
        embed_colors.InsertNextTuple3(
            color[lidx, 0], color[lidx, 1], color[lidx, 2])
        cluster_colors.InsertNextTuple1(int(cluster_idx[lidx]))
        embed_data.InsertNextTupleValue(embed[lidx, :])
        row_sum_data.InsertNextTuple1(float(estimated_row_sum[lidx]))
        col_sum_data.InsertNextTuple1(float(computed_column_sum[lidx]))

    output_polydata.GetCellData().AddArray(embed_data)
    output_polydata.GetCellData().AddArray(cluster_colors)
    output_polydata.GetCellData().AddArray(embed_colors)
    output_polydata.GetCellData().AddArray(row_sum_data)
    output_polydata.GetCellData().AddArray(col_sum_data)

    output_polydata.GetCellData().SetActiveScalars('EmbeddingColor')

    return output_polydata


def _embed_to_rgb(embed):
    """

    Get colors from embedding (first 3 coordinates corresponding to
    max eigenvalues). Output minimum color is 0,0,0, max is
    255,255,255, (actually 180 due to brightness constraint) and all
    output values (brightnesses) are equal.

     """

    # make sure nothing is NaN
    testval = numpy.count_nonzero(numpy.isnan(embed))
    if testval:
        print "<cluster.py> Warning: Consider increasing sigma or using the Mean distance. NaN values encountered in embedding."
        print "Number of NaN values:", testval
        embed[numpy.isnan(embed)] = 0.0

    # first 3 components of embedding
    color = embed[:, 0:3]

    # normalize all colors to length 1
    color_len = numpy.sqrt(numpy.sum(numpy.power(color, 2), 1))
    #print numpy.min(color_len)
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


def output_and_quality_control_cluster_atlas(atlas, output_polydata_s, subject_fiber_list, input_polydatas, number_of_subjects, outdir, cluster_numbers_s, color, embed, number_of_fibers_to_display, testing=False, verbose=False, render_images=True):

    """Save the output in our atlas format for automatic labeling of clusters.

    First save the atlas.vtp and atlas.p datasets. This is the data used to
    label a new subject.  Also write the polydata with cluster indices
    saved as cell data. This is a one-file output option for clusters.
    Finally, save some quality control metrics and save the atlas
    clusters as individual polydatas. This is used to set up a mrml
    hierarchy file and to visualize the output in Slicer. This data is
    not used to label a new subject.
    """

    # Write additional output for software testing if code has changed (requires random seed to be constant)
    if testing:
        expected_results_file = os.path.join(outdir, 'test_cluster_atlas_numbers.pkl')
        pickle.dump(cluster_numbers_s, open(expected_results_file, 'wb'))
        expected_results_file = os.path.join(outdir, 'test_cluster_atlas_colors.pkl')
        pickle.dump(color, open(expected_results_file, 'wb'))
        expected_results_file = os.path.join(outdir, 'test_cluster_atlas_embeddings.pkl')
        pickle.dump(embed, open(expected_results_file, 'wb'))
         
    # Save the output in our atlas format for automatic labeling of full brain datasets.
    # This is the data used to label a new subject
    atlas.save(outdir,'atlas')

    # Write the polydata with cluster indices saved as cell data
    fname_output = os.path.join(outdir, 'clustered_whole_brain.vtp')
    io.write_polydata(output_polydata_s, fname_output)

    # output summary file to save information about all subjects
    subjects_qc_fname = os.path.join(outdir, 'input_subjects.txt')
    subjects_qc_file = open(subjects_qc_fname, 'w')
    outstr = "Subject_idx\tSubject_ID\tfilename\n"
    subjects_qc_file.write(outstr)
    idx = 1
    for fname in input_polydatas:
        subject_id = os.path.splitext(os.path.basename(fname))[0]
        outstr =  str(idx) + '\t' + str(subject_id) + '\t' + str(fname) + '\n'
        subjects_qc_file.write(outstr)
        idx += 1
    subjects_qc_file.close()

    # output summary file to save information about all clusters
    clusters_qc_fname = os.path.join(outdir, 'cluster_quality_control.txt')
    clusters_qc_file = open(clusters_qc_fname, 'w')

    # Figure out how many subjects in each cluster (ideally, most subjects in most clusters)
    subjects_per_cluster = list()
    percent_subjects_per_cluster = list()
    fibers_per_cluster = list()
    mean_fiber_len_per_cluster = list()
    std_fiber_len_per_cluster = list()
    mean_fibers_per_subject_per_cluster = list()
    std_fibers_per_subject_per_cluster = list()

    # find out length of each fiber
    fiber_length, step_size = filter.compute_lengths(output_polydata_s)

    # loop over each cluster and compute quality control metrics
    cluster_indices = range(atlas.centroids.shape[0])
    for cidx in cluster_indices:
        cluster_mask = (cluster_numbers_s==cidx) 
        subjects_per_cluster.append(len(set(subject_fiber_list[cluster_mask])))
        fibers_per_subject = list()
        for sidx in range(number_of_subjects):
            fibers_per_subject.append(list(subject_fiber_list[cluster_mask]).count(sidx))
        mean_fibers_per_subject_per_cluster.append(numpy.mean(numpy.array(fibers_per_subject)))
        std_fibers_per_subject_per_cluster.append(numpy.std(numpy.array(fibers_per_subject)))
        mean_fiber_len_per_cluster.append(numpy.mean(fiber_length[cluster_mask]))
        std_fiber_len_per_cluster.append(numpy.std(fiber_length[cluster_mask]))

    percent_subjects_per_cluster = numpy.divide(numpy.array(subjects_per_cluster),float(number_of_subjects))

    # Save output quality control information
    print "<cluster.py> Saving cluster quality control information file."
    clusters_qc_file = open(clusters_qc_fname, 'w')
    print >> clusters_qc_file, 'cluster_idx','\t', 'number_subjects','\t', 'percent_subjects','\t', 'mean_length','\t', 'std_length','\t', 'mean_fibers_per_subject','\t', 'std_fibers_per_subject'
    for cidx in cluster_indices:
        print >> clusters_qc_file, cidx + 1,'\t', subjects_per_cluster[cidx],'\t', percent_subjects_per_cluster[cidx] * 100.0,'\t', \
            mean_fiber_len_per_cluster[cidx],'\t', std_fiber_len_per_cluster[cidx],'\t', \
            mean_fibers_per_subject_per_cluster[cidx],'\t', std_fibers_per_subject_per_cluster[cidx]

    clusters_qc_file.close()

    if HAVE_PLT:
        print "<cluster.py> Saving subjects per cluster histogram."
        fig, ax = plt.subplots()
        counts = numpy.zeros(num_of_subjects+1)
        counts[:numpy.max(subjects_per_cluster)+1] = numpy.bincount(subjects_per_cluster)
        ax.bar(range(num_of_subjects + 1), counts, width=1, align='center')
        ax.set(xlim=[-1, num_of_subjects + 1])
        plt.title('Histogram of Subjects per Cluster')
        plt.xlabel('subjects per cluster')
        plt.ylabel('number of clusters')
        plt.savefig( os.path.join(outdir, 'subjects_per_cluster_hist.pdf'))
        plt.close()
        
    # Save the entire combined atlas as individual clusters for visualization
    # and labeling/naming of structures. This will include all of the data
    # that was clustered to make the atlas.

    # Figure out file name and mean color for each cluster, and write the individual polydatas
    print "<cluster.py> Beginning to save individual clusters as polydata files. TOTAL CLUSTERS:", len(cluster_indices),
    fnames = list()
    cluster_colors = list()
    cluster_sizes = list()
    cluster_fnames = list()
    for c in cluster_indices:
        print c,
        mask = cluster_numbers_s == c
        cluster_size = numpy.sum(mask)
        cluster_sizes.append(cluster_size)
        pd_c = filter.mask(output_polydata_s, mask,verbose=verbose)
        # color by subject so we can see which one it came from
        filter.add_point_data_array(pd_c, subject_fiber_list[mask], "Subject_ID")
        # Save hemisphere information into the polydata
        farray = fibers.FiberArray()
        farray.hemispheres = True
        farray.hemisphere_percent_threshold = 0.90
        farray.convert_from_polydata(pd_c, points_per_fiber=50)
        filter.add_point_data_array(pd_c, farray.fiber_hemisphere, "Hemisphere")
        # The clusters are stored starting with 1, not 0, for user friendliness.
        fname_c = 'cluster_{0:05d}.vtp'.format(c+1)
        # save the filename for writing into the MRML file
        fnames.append(fname_c)
        # prepend the output directory
        fname_c = os.path.join(outdir, fname_c)
        #print fname_c
        io.write_polydata(pd_c, fname_c)
        cluster_fnames.append(fname_c)
        if cluster_size > 0:
            color_c = color[mask,:]
            cluster_colors.append(numpy.mean(color_c,0))
        else:
            cluster_colors.append([0,0,0])
        del pd_c
    print "\n<cluster.py> Finishes saving individual clusters as polydata files."

    # Notify user if some clusters empty
    empty_count = 0
    for sz, fname in zip(cluster_sizes,cluster_fnames):
        if sz == 0:
            print sz, ":", fname
            empty_count += 1
    if empty_count:
        print "<cluster.py> Warning. Empty clusters found:", empty_count

    cluster_sizes = numpy.array(cluster_sizes)
    print "<cluster.py> Mean number of fibers per cluster:", numpy.mean(cluster_sizes), "Range:", numpy.min(cluster_sizes), "..", numpy.max(cluster_sizes)

    # Estimate subsampling ratio to display approximately number_of_fibers_to_display total fibers in 3D Slicer
    number_fibers = len(cluster_numbers_s)
    if number_fibers < number_of_fibers_to_display:
        ratio = 1.0
    else:
        ratio = number_of_fibers_to_display / number_fibers
    print "<cluster.py> Subsampling ratio for display of", number_of_fibers_to_display, "total fibers estimated as:", ratio

    # Write the MRML file into the directory where the polydatas were already stored
    fname = os.path.join(outdir, 'clustered_tracts.mrml')
    mrml.write(fnames, numpy.around(numpy.array(cluster_colors), decimals=3), fname, ratio=ratio)

    # Also write one with 100% of fibers displayed
    fname = os.path.join(outdir, 'clustered_tracts_display_100_percent.mrml')
    mrml.write(fnames, numpy.around(numpy.array(cluster_colors), decimals=3), fname, ratio=1.0)
    
    # View the whole thing in jpg format for quality control
    if render_images:
        print '<cluster.py> Rendering and saving images of cluster atlas.'
        ren = render.render(output_polydata_s, 1000, data_mode='Cell', data_name='EmbeddingColor', verbose=verbose)
        ren.save_views(outdir, verbose=verbose)
        del ren

# -*- coding: utf-8 -*-

""" cluster.py

implementations of fiber clustering

"""
import colorsys
import os
import pickle
from pprint import pprint

import numpy as np
import scipy.cluster.hierarchy
import scipy.cluster.vq
import vtk
from joblib import Parallel, delayed

from whitematteranalysis.utils.opt_pckg import optional_package

from . import fibers, filter, io, mrml, render, similarity

matplotlib, have_mpl, _ = optional_package("matplotlib")
plt, _, _ = optional_package("matplotlib.pyplot")

if have_mpl:
    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use("Agg")
else:
    warnings.warn(matplotlib._msg)
    warnings.warn("Cannot plot quality control data.")


# This did not work better. Leave here for future testing if of interest
if 0:
    try:    
        from sklearn import metrics
        from sklearn.cluster import AffinityPropagation
    except ImportError:
        SKLEARN = 0
        print(f"<{os.path.basename(__file__)}> Failed to import sklearn, cannot use affinity propagation.")
        print(f"<{os.path.basename(__file__)}> Please install sklearn for this functionality.")
    
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
            print(f"Error: Atlas directory {directory} does not exist or is not a directory.")
            raise f"<{os.path.basename(__file__)}> I/O error"
        
        fname_base = os.path.join(directory,atlas_name)
        fname_atlas = f'{fname_base}.p'
        fname_polydata = f'{fname_base}.vtp'
        
        if not os.path.exists(fname_atlas):
            print(f"Error: Atlas file {fname_atlas} does not exist.")
            raise f"<{os.path.basename(__file__)}> I/O error"
        if not os.path.exists(fname_polydata):
            print(f"Error: Atlas file {fname_polydata} does not exist.")
            raise f"<{os.path.basename(__file__)}> I/O error"

        try:
            atlas = pickle.load(open(fname_atlas,'rb'))
        except UnicodeDecodeError:
            atlas = pickle.load(open(fname_atlas,'rb'), encoding="latin1")
        atlas.nystrom_polydata = io.read_polydata(fname_polydata)
        print(f"<{os.path.basename(__file__)}> Loaded atlas {atlas_name} from {directory}.")
        print(f"<{os.path.basename(__file__)}> Atlas Nystrom polydata sample: {atlas.nystrom_polydata.GetNumberOfLines()}\nAtlas size: {atlas.pinv_A.shape}\nAtlas number of eigenvectors: {atlas.number_of_eigenvectors}")

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
    #mask_matrix = np.ones(distance_matrix.shape)
    #mask_matrix = np.triu(mask_matrix)

    #distance_matrix_triu = distance_matrix(np.nonzero(mask_matrix))

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
   eigval, eigvec = np.linalg.eig(A)
   val = np.matrix(np.maximum(eigval,epsilon))
   vec = np.matrix(eigvec)
   T = 1/(np.multiply(vec,vec) * val.T)
   T = np.matrix(np.sqrt(np.diag(np.array(T).reshape(n) )))
   B = T * vec * np.diag(np.array(np.sqrt(val)).reshape(n))
   out = B*B.T
   return(np.asarray(out))
   
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
    print(f"<{os.path.basename(__file__)}> Starting spectral clustering.")
    print(f"<{os.path.basename(__file__)}> Number of input fibers: {number_fibers}")
    print(f"<{os.path.basename(__file__)}> Number of clusters: {number_of_clusters}")

    if number_fibers == 0:
        print(f"<{os.path.basename(__file__)}> ERROR: Cannot cluster polydata with 0 fibers.")
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
        nystrom_mask = np.array(nystrom_mask)
        # make sure it's boolean or 0 and 1
        test = np.max(nystrom_mask) == 1.0
        if not test:
            print(f"<{os.path.basename(__file__)}> ERROR: Nystrom data mask is may not be Boolean. Max value is not 1.0/True.")
            raise AssertionError
        # make sure it's large enough
        test = sum(nystrom_mask) >= 100
        if not test:
            print(f"<{os.path.basename(__file__)}> ERROR: Nystrom data mask is smaller than 100.")
            raise AssertionError
        # make sure its size matches the polydata input
        test = len(nystrom_mask) == number_fibers
        if not test:
            print(f"<{os.path.basename(__file__)}> ERROR: Nystrom data mask size does not match polydata number of lines.")
            raise AssertionError
        
        # Separate the Nystrom sample and the rest of the data.
        polydata_m = filter.mask(input_polydata, nystrom_mask, verbose=False)
        atlas.nystrom_polydata = polydata_m
        polydata_n = filter.mask(input_polydata, nystrom_mask == False, verbose=False)
        sz = polydata_m.GetNumberOfLines()
        print(f'<{os.path.basename(__file__)}> Using Nystrom approximation. Subset size: {sz} / {number_fibers}')
        # Determine ordering to get embedding to correspond to original input data.
        reorder_embedding = np.concatenate((np.where(nystrom_mask)[0], np.where(~nystrom_mask)[0]))
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
        print(f"<{os.path.basename(__file__)}> Range of values in A: {np.min(A)} {np.max(A)}")
        print(f"<{os.path.basename(__file__)}> Range of values in B: {np.min(B)} {np.max(B)}")
        
    else:
        # Calculate all fiber similarities
        A = \
            _pairwise_similarity_matrix(input_polydata, threshold,
                                    sigma, number_of_jobs, landmarks, distance_method, bilateral)

        atlas.nystrom_polydata = input_polydata
        # sanity check
        print(f"<{os.path.basename(__file__)}> Range of values in A: {np.min(A)} {np.max(A)}")
        
    testval = np.max(A-A.T)
    if not testval == 0.0:
        if testval > 1e-10:
            print(f"<{os.path.basename(__file__)}> ERROR: A matrix is not symmetric.")
            raise AssertionError
        else:
            print(f"<{os.path.basename(__file__)}> Maximum of A - A^T: {testval}")
        # Ensure that A is symmetric
        A = np.divide(A+A.T, 2.0)
        
    testval = np.min(A)
    if not testval > 0.0:
        print(f"<{os.path.basename(__file__)}> ERROR: A matrix is not positive.")
        print(f"<{os.path.basename(__file__)}> Minimum value in A: {testval}")
        if testval < 0.0:
            raise AssertionError

    # Outliers will have low measured (or estimated) row sums. Detect outliers in A:
    # to turn off for testing: outlier_std_threshold = np.inf
    row_sum_A_initial = np.sum(A, axis=0) + np.sum(B.T, axis=0)
    print(f"<{os.path.basename(__file__)}> Initial similarity (row) sum A: {np.mean(row_sum_A_initial)} {np.std(row_sum_A_initial)} {np.min(row_sum_A_initial)}")
    atlas.outlier_std_threshold = outlier_std_threshold
    atlas.row_sum_threshold_for_rejection = np.mean(row_sum_A_initial) - outlier_std_threshold*np.std(row_sum_A_initial)
    bad_idx = np.nonzero(row_sum_A_initial < atlas.row_sum_threshold_for_rejection)[0]
    reject_A = bad_idx
    print(f"<{os.path.basename(__file__)}> Rejecting n= {len(bad_idx)} / {sz} fibers > {outlier_std_threshold} standard deviations below the mean total fiber similarity")

    A = np.delete(A,reject_A,0)
    A = np.delete(A,reject_A,1)
    #print A.shape, B.shape
    B = np.delete(B,reject_A,0)
    #print A.shape, B.shape, reorder_embedding.shape
                    
    # Ensure that A is positive definite.
    if pos_def_approx:
        e_val, e_vec = np.linalg.eigh(A)
        print(f"<{os.path.basename(__file__)}> Eigenvalue range of A: {e_val[0]} {e_val[-1]}")
        A2 = nearPSD(A)
        e_val, e_vec = np.linalg.eigh(A2)
        print(f"<{os.path.basename(__file__)}> Eigenvalue range of nearest PSD matrix to A: {e_val[0]} e{_val[-1]}")
        testval = np.max(A-A2)
        if not testval == 0.0:
            print(f"<{os.path.basename(__file__)}> A matrix differs by PSD matrix by maximum of: {testval}")
            if testval > 0.25:
                print(f"<{os.path.basename(__file__)}> ERROR: A matrix changed by more than 0.25.")
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
            atlas.row_sum_1 = np.sum(A, axis=0) + np.sum(B.T, axis=0)
            #print(f"<{os.path.basename(__file__)}> A size: {A.shape}")
            #print(f"<{os.path.basename(__file__)}> B size: {B.shape}")
            #print(f"<{os.path.basename(__file__)}> A-B matrix row sums range (should be > 0): {np.min(atlas.row_sum_1)} {np.max(atlas.row_sum_1)}")

            # Approximate the sum of the rest of the data (including C)
            # These are weighted sums of the columns we did compute
            # where the weight depends on how similar that fiber 
            # was to each path in A.  This uses the dual basis
            # of the columns in A.
            # Approximate the inverse of A for dual basis
            #print(f"<{os.path.basename(__file__)}> Using numpy linalg pinv A")
            atlas.pinv_A = np.linalg.pinv(A)

            #e_val, e_vec = np.linalg.eigh(atlas.pinv_A)
            #print(f"<{os.path.basename(__file__)}> test of non-normalized A pseudoinverse Eigenvalue range: {e_val[0]} {e_val[-1]}")

            # row sum formula:
            # dhat = [a_r + b_r; b_c + B^T*A-1*b_r]
            # this matrix is A^-1 * b_r, where b_r are the row sums of B
            # matlab was: atlas.approxRowSumMatrix = sum(B',1)*atlas.pseudoInverseA;
            atlas.row_sum_matrix = np.dot(np.sum(B.T, axis=0), atlas.pinv_A)
            #test = np.sum(B.T, axis=0)
            #print(f"<{os.path.basename(__file__)}> B column sums range (should be > 0): {np.min(test)} {np.max(test)}")
            print(f"<{os.path.basename(__file__)}> Range of row sum weights: {np.min(atlas.row_sum_matrix)} {np.max(atlas.row_sum_matrix)}")
            #print(f"<{os.path.basename(__file__)}> First 10 entries in weight matrix: {atlas.row_sum_matrix[0:10]}")
            #test = np.dot(atlas.row_sum_matrix, B)
            #print(f"<{os.path.basename(__file__)}> Test partial sum estimation for B: {np.min(test)} {np.max(test)}")
            #del test
            
            # row sum estimate for current B part of the matrix
            row_sum_2 = np.sum(B, axis=0) + \
                np.dot(atlas.row_sum_matrix, B)
            print(f"<{os.path.basename(__file__)}> Row sum check (min/max, should be > 0) A: {np.min(atlas.row_sum_1)} {np.median(atlas.row_sum_1)} {np.max(atlas.row_sum_1)} B: {np.min(row_sum_2)} {np.median(row_sum_2)} {np.max(row_sum_2)}")

            # reject outliers in B
            bad_idx = np.nonzero(row_sum_2 < atlas.row_sum_threshold_for_rejection)[0]
            reject_B = bad_idx
            print(f"<{os.path.basename(__file__)}> Rejecting n= {len(bad_idx)} / {B.shape[1]} fibers > {outlier_std_threshold} standard deviations below the mean total fiber similarity")
            row_sum_2 = np.delete(row_sum_2,reject_B)
            B = np.delete(B,reject_B,1)

            print(f"<{os.path.basename(__file__)}> After outlier rejection A: {A.shape} B: {B.shape}")
            print(f"<{os.path.basename(__file__)}> Row sum check (min/max, should be > 0) A: {np.min(atlas.row_sum_1)} {np.median(atlas.row_sum_1)} {np.max(atlas.row_sum_1)} B: {np.min(row_sum_2)} {np.median(row_sum_2)} {np.max(row_sum_2)}")

            # Separate the Nystrom sample and the rest of the data after removing outliers
            nystrom_mask_2 = nystrom_mask
            midxA = np.nonzero(nystrom_mask_2)[0]
            nystrom_mask_2[midxA[reject_A]] = False
            not_nystrom_mask = nystrom_mask == False
            not_nystrom_mask[midxA[reject_A]] = False
            midxB = np.nonzero(not_nystrom_mask)[0]
            not_nystrom_mask[midxB[reject_B]] = False

            polydata_m = filter.mask(input_polydata, nystrom_mask_2, verbose=False)
            atlas.nystrom_polydata = polydata_m
            polydata_n = filter.mask(input_polydata, not_nystrom_mask, verbose=False)
            output_polydata = filter.mask(input_polydata, np.add(nystrom_mask_2, not_nystrom_mask),verbose=False)
            sz = polydata_m.GetNumberOfLines()
            number_fibers = output_polydata.GetNumberOfLines()
            print(f'<{os.path.basename(__file__)}> Using Nystrom approximation. Subset size (A): {sz} / {number_fibers} B: {polydata_n.GetNumberOfLines()}')
            # Determine ordering to get embedding to correspond to original input data.
            # reject outliers from masks
            reject_idx = np.concatenate((midxA[reject_A],midxB[reject_B]))
            nystrom_mask_2 = np.delete(nystrom_mask_2,reject_idx)
            not_nystrom_mask = np.delete(not_nystrom_mask,reject_idx)
            #print "hi after mask:", reorder_embedding.shape, np.sum(nystrom_mask_2), np.sum(not_nystrom_mask)
            reorder_embedding = np.concatenate((np.where(nystrom_mask_2)[0], np.where(not_nystrom_mask)[0]))
            #print "hi after embed reorder calc:", reorder_embedding.shape, np.max(reorder_embedding), np.min(reorder_embedding)
           
            # in case of negative row sum estimation
            if any(row_sum_2<=0):
                print(f"<{os.path.basename(__file__)}> Warning: Consider increasing sigma or using the Mean distance. negative row sum approximations.")
                print(f"Number of negative row sums: {np.count_nonzero(row_sum_2 <= 0)}")
                #row_sum_2[row_sum_2<0] = 0.1

            # save for testing
            column_sum = np.concatenate((np.sum(A, axis=1) , np.sum(B.T, axis=1)))

            # normalized cuts normalization
            row_sum = np.concatenate((atlas.row_sum_1, row_sum_2))
            dhat = np.sqrt(np.divide(1, row_sum))
            #dhat = np.sqrt(np.divide(1, np.concatenate((atlas.row_sum_1, row_sum_2))))

            A = \
                np.multiply(A, np.outer(dhat[0:sz], dhat[0:sz].T))
            B = \
                np.multiply(B, np.outer(dhat[0:sz], dhat[sz:].T))

        else:
            # normalized cuts normalization using row (same as column) sums
            row_sum = np.sum(A, axis=0)
            print(f"<{os.path.basename(__file__)}> A matrix row sums range (should be > 0): {np.min(row_sum)} {np.max(row_sum)}")
            dhat = np.divide(1, np.sqrt(row_sum))
            A = \
                np.multiply(A, np.outer(dhat, dhat.T))

    # 3) Compute eigenvectors for use in spectral embedding
    print(f'<{os.path.basename(__file__)}> Calculating eigenvectors of similarity matrix A...')
    atlas.e_val, atlas.e_vec = np.linalg.eigh(A)
    print(f'<{os.path.basename(__file__)}> Done calculating eigenvectors.')
    print(f"<{os.path.basename(__file__)}> Eigenvalue range: {atlas.e_val[0]} {atlas.e_val[-1]}")
    # Check how well our chosen number of eigenvectors models the data
    power = np.cumsum(atlas.e_val[::-1]) / np.sum(atlas.e_val)
    print(f"<{os.path.basename(__file__)}> Power from chosen number of eigenvectors ({number_of_eigenvectors}) {power[number_of_eigenvectors]}")
    print(f'<{os.path.basename(__file__)}> Top eigenvalues: {atlas.e_val[::-1][1:number_of_eigenvectors]}')

    # 4) Compute embedding using eigenvectors
    print(f'<{os.path.basename(__file__)}> Compute embedding using eigenvectors.')
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
        V = np.concatenate((atlas.e_vec, \
                                  np.dot(np.dot(B.T, atlas.e_vec), \
                                                np.diag(np.divide(1.0, atlas.e_val)))))

        # normalize estimated eigenvectors to have length of one
        # matlab was:
        # atlas.eigenvectorLengthToNormalize=sqrt(sum(V.*V));
        # V=V./repmat(atlas.eigenvectorLengthToNormalize,length(V),1);
        atlas.e_vec_norm = np.sum(np.multiply(V, V),0)
        V = np.divide(V, atlas.e_vec_norm)

        # Normalize each embedding vector by first eigenvector. Matlab code was:
        # for i = 2:embedLength+1
        #    embedding(:,i-1) = V(:,i)./V(:,1);
        # end
        # This eigenvector corresponds to an eigenvalue of 1, since row sums are 1.
        # The other option from the literature was to use this:
        # embedding_i,j = V_i+i,j./sqrt(D_j,j)
        embed = np.zeros((number_fibers, number_of_eigenvectors))
        #print "Hi 3:", embed.shape, number_fibers, reorder_embedding.shape, V.shape, A.shape, B.shape
        for i in range(0, number_of_eigenvectors):
            embed[reorder_embedding,i] = np.divide(V[:,-(i+2)], V[:,-1])
    else:
        embed = atlas.e_vec[:, -number_of_eigenvectors - 2: -2]
        embed = np.divide(embed.T, atlas.e_vec[:, -1]).T
        # reverse order of embedding so highest eigenvalue
        # information is first
        embed = embed[:, ::-1]


    # Default is always k-means. Other code is just left for testing. Did not improve results.
    #centroid_finder = 'AffinityPropagation'
    centroid_finder = 'K-means'
    
    # 5) Find clusters using k-means in embedding space.
    cluster_metric = None
    if centroid_finder == 'K-means':
        print(f'<{os.path.basename(__file__)}> K-means clustering in embedding space.')
        centroids, cluster_metric = scipy.cluster.vq.kmeans2(embed, number_of_clusters, minit='points')
        # sort centroids by first eigenvector order
        # centroid_order = np.argsort(centroids[:,0])
        # sort centroids according to colormap and save them in this order in atlas
        color = _embed_to_rgb(centroids)
        centroid_order = render.argsort_by_jet_lookup_table(color)
        atlas.centroids = centroids[centroid_order,:]
        cluster_idx, dist = scipy.cluster.vq.vq(embed, atlas.centroids)
        #print(f"<{os.path.basename(__file__)}> Distortion metric: {cluster_metric}")
        if 0:
            # This is extremely slow, but leave code here if ever wanted for testing
            cluster_metric = metrics.silhouette_score(embed, cluster_idx, metric='sqeuclidean')
            print(f"Silhouette Coefficient: {cluster_metric:0.3f}")
 
    else:
        raise NotImplementedError(
            f"Workflow not implemented for centroid finder: {centroid_finder}.")
        ## # This found fewer clusters than we need to represent the anatomy well
        ## # Leave code here in case wanted in future for more testing.
        ## print(f'<{os.path.basename(__file__)}> Affinity Propagation clustering in embedding space.')
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
    print(f'<{os.path.basename(__file__)}> Done spectral clustering, returning results.')
    # visualize embedding coordinates as RGB
    embed2 = embed
    #embed2[np.isnan(embed)] = 0.0
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
        cluster_indices = np.zeros(cluster_indices_vtk.GetNumberOfTuples())
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
    row_sum_2 = np.sum(B, axis=0) + \
        np.dot(atlas.row_sum_matrix, B)

    # in case of negative row sum estimation (this should not ever happen)
    if any(row_sum_2<=0):
        print(f"<{os.path.basename(__file__)}> Warning: Consider increasing sigma or using the Mean distance. negative row sum approximations.")
        print("Number of negative row sums:", np.count_nonzero(row_sum_2<=0))
        row_sum_2[row_sum_2<0] = 0.1

             
    # normalized cuts normalization
    row_sum = np.concatenate((atlas.row_sum_1, row_sum_2))
    dhat = np.sqrt(np.divide(1, row_sum))
    #dhat = np.sqrt(np.divide(1, np.concatenate((atlas.row_sum_1, row_sum_2))))
    B = \
        np.multiply(B, np.outer(dhat[0:sz], dhat[sz:].T))

    # 3) Compute eigenvectors for use in spectral embedding
    # <done already in atlas creation>

    # 4) Compute embedding using eigenvectors
    V = np.dot(np.dot(B.T, atlas.e_vec), \
                      np.diag(np.divide(1.0, atlas.e_val)))
    V = np.divide(V, atlas.e_vec_norm)
    embed = np.zeros((number_fibers, atlas.number_of_eigenvectors))
    for i in range(0, atlas.number_of_eigenvectors):
        embed[:,i] = np.divide(V[:,-(i+2)], V[:,-1])

    # 5) Find clusters using k-means in embedding space.
    # Actually: LABEL using centroids from atlas
    print(f'<{os.path.basename(__file__)}> Cluster labeling in embedding space.')
    cluster_idx, dist = scipy.cluster.vq.vq(embed, atlas.centroids)

    # 6) Output results.
    print(f'<{os.path.basename(__file__)}> Done labeling clusters, returning results.')
    # visualize embedding coordinates as RGB
    color = _embed_to_rgb(embed)
    # set up polydata with clustering output info.
    # for now modify input polydata by adding two arrays
    output_polydata = input_polydata
    output_polydata = \
        _format_output_polydata(output_polydata, cluster_idx, color, embed, row_sum, np.sum(B.T, axis=1))

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
        distances = np.array(distances)

    else:
        
        fiber_array_n = fibers.FiberArray()
        fiber_array_n.convert_from_polydata(input_polydata_n, points_per_fiber=15)
        fiber_array_m = fibers.FiberArray()
        fiber_array_m.convert_from_polydata(input_polydata_m, points_per_fiber=15)

        if landmarks_n is None:
            landmarks_n = np.zeros((fiber_array_n.number_of_fibers,3))
    
        # pairwise distance matrix
        all_fibers_n = list(range(0, fiber_array_n.number_of_fibers))

        distances = Parallel(n_jobs=number_of_jobs,
                             verbose=0)(
            delayed(similarity.fiber_distance)(
                fiber_array_n.get_fiber(lidx),
                fiber_array_m,
                threshold, distance_method=distance_method,
                fiber_landmarks=landmarks_n[lidx,:], 
                landmarks=landmarks_m, bilateral=bilateral)
            for lidx in all_fibers_n)

        distances = np.array(distances).T

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
        distances = np.array(distances)

    else:

        fiber_array = fibers.FiberArray()
        fiber_array.convert_from_polydata(input_polydata, points_per_fiber=15)

        # pairwise distance matrix
        all_fibers = list(range(0, fiber_array.number_of_fibers))

        if landmarks is None:
            landmarks2 = np.zeros((fiber_array.number_of_fibers,3))
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

        distances = np.array(distances)

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
    #print "This should be 1.0: ", np.min(np.diag(similarity_matrix))
    #print  np.min(np.diag(similarity_matrix)) == 1.0
    # test
    if __debug__:
        # this tests that on-diagonal elements are all 1
        test = np.min(np.diag(similarity_matrix)) == 1.0
        if not test:
            print(f"<{os.path.basename(__file__)}> ERROR: On-diagonal elements are not all 1.0.")
            print(f" Minimum on-diagonal value: {np.min(np.diag(similarity_matrix))}")
            print(f" Maximum on-diagonal value: {np.max(np.diag(similarity_matrix))}")
            print(f" Mean value: {np.mean(np.diag(similarity_matrix))}")
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
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 7.1):
          embed_data.InsertNextTypedTuple(embed[lidx, :])
        else:
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
    testval = np.count_nonzero(np.isnan(embed))
    if testval:
        print(f"<{os.path.basename(__file__)}> Warning: Consider increasing sigma or using the Mean distance. NaN values encountered in embedding.")
        print("Number of NaN values:", testval)
        embed[np.isnan(embed)] = 0.0

    # first 3 components of embedding
    color = embed[:, 0:3]

    # normalize all colors to length 1
    color_len = np.sqrt(np.sum(np.power(color, 2), 1))
    #print np.min(color_len)
    color_len += 0.0001
    color = np.divide(color.T, color_len).T

    # convert to RGB: use colors from 0 to 255 for unsigned char (no
    # LUT in vtk). Before color components ranged from -1 to +1, now
    # from 0 to 255
    color = 127.5 + (color * 127.5)

    # Avoid dark black fibers that the cluster color is hard to see.
    # Convert to hsv, and make all brightness values the same.  That
    # way shadows show geometry, not cluster number.
    hue = np.zeros(color.shape[0])
    sat = np.zeros(color.shape[0])
    val = np.zeros(color.shape[0])
    for c_idx in range(0, color.shape[0]):
        hue[c_idx], sat[c_idx], val[c_idx] = \
            colorsys.rgb_to_hsv(
            color[c_idx, 0], color[c_idx, 1], color[c_idx, 2])

    val2 = np.ones(val.shape) * 180
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
        outstr = f'{str(idx)}\t{str(subject_id)}\t{str(fname)}\n'
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
    cluster_indices = list(range(atlas.centroids.shape[0]))
    for cidx in cluster_indices:
        cluster_mask = (cluster_numbers_s==cidx) 
        subjects_per_cluster.append(len(set(subject_fiber_list[cluster_mask])))
        fibers_per_subject = list()
        for sidx in range(number_of_subjects):
            fibers_per_subject.append(list(subject_fiber_list[cluster_mask]).count(sidx))
        mean_fibers_per_subject_per_cluster.append(np.mean(np.array(fibers_per_subject)))
        std_fibers_per_subject_per_cluster.append(np.std(np.array(fibers_per_subject)))
        mean_fiber_len_per_cluster.append(np.mean(fiber_length[cluster_mask]))
        std_fiber_len_per_cluster.append(np.std(fiber_length[cluster_mask]))

    percent_subjects_per_cluster = np.divide(np.array(subjects_per_cluster),float(number_of_subjects))

    # Save output quality control information
    print(f"<{os.path.basename(__file__)}> Saving cluster quality control information file.")
    clusters_qc_file = open(clusters_qc_fname, 'w')
    print('cluster_idx\tnumber_subjects\tpercent_subjects\tmean_length\tstd_length\tmean_fibers_per_subject\tstd_fibers_per_subject', file=clusters_qc_file)
    for cidx in cluster_indices:
        print(f'{cidx + 1}\t{subjects_per_cluster[cidx]}\t{percent_subjects_per_cluster[cidx] * 100.0}\t{mean_fiber_len_per_cluster[cidx]}\t{std_fiber_len_per_cluster[cidx]}\t{mean_fibers_per_subject_per_cluster[cidx]}\t{std_fibers_per_subject_per_cluster[cidx]}', file=clusters_qc_file)

    clusters_qc_file.close()

    if have_mpl:
        print(f"<{os.path.basename(__file__)}> Saving subjects per cluster histogram.")
        fig, ax = plt.subplots()
        counts = np.zeros(number_of_subjects+1)
        counts[:np.max(subjects_per_cluster)+1] = np.bincount(subjects_per_cluster)
        ax.bar(list(range(number_of_subjects + 1)), counts, width=1, align='center')
        ax.set(xlim=[-1, number_of_subjects + 1])
        plt.title('Histogram of Subjects per Cluster')
        plt.xlabel('subjects per cluster')
        plt.ylabel('number of clusters')
        plt.savefig( os.path.join(outdir, 'subjects_per_cluster_hist.pdf'))
        plt.close()
        
    # Save the entire combined atlas as individual clusters for visualization
    # and labeling/naming of structures. This will include all of the data
    # that was clustered to make the atlas.

    # Figure out file name and mean color for each cluster, and write the individual polydatas
    pd_c_list = mask_all_clusters(output_polydata_s, cluster_numbers_s, len(cluster_indices), preserve_point_data=True,
                                  preserve_cell_data=True, verbose=False)
    print(f"<{os.path.basename(__file__)}> Beginning to save individual clusters as polydata files. TOTAL CLUSTERS: {len(cluster_indices)}", end=' ')
    fnames = list()
    cluster_colors = list()
    cluster_sizes = list()
    cluster_fnames = list()
    for c in cluster_indices:
        print(c, end=' ')
        mask = cluster_numbers_s == c
        cluster_size = np.sum(mask)
        cluster_sizes.append(cluster_size)
        #pd_c = filter.mask(output_polydata_s, mask,verbose=verbose)
        pd_c = pd_c_list[c]
        # color by subject so we can see which one it came from
        filter.add_point_data_array(pd_c, subject_fiber_list[mask], "Subject_ID")
        # Save hemisphere information into the polydata
        farray = fibers.FiberArray()
        farray.hemispheres = True
        farray.hemisphere_percent_threshold = 0.90
        farray.convert_from_polydata(pd_c, points_per_fiber=50)
        filter.add_point_data_array(pd_c, farray.fiber_hemisphere, "Hemisphere")
        # The clusters are stored starting with 1, not 0, for user friendliness.
        fname_c = f'cluster_{c+1:05d}.vtp'
        # save the filename for writing into the MRML file
        fnames.append(fname_c)
        # prepend the output directory
        fname_c = os.path.join(outdir, fname_c)
        #print fname_c
        io.write_polydata(pd_c, fname_c)
        cluster_fnames.append(fname_c)
        if cluster_size > 0:
            color_c = color[mask,:]
            cluster_colors.append(np.mean(color_c,0))
        else:
            cluster_colors.append([0,0,0])
        del pd_c
    print(f"\n<{os.path.basename(__file__)}> Finished saving individual clusters as polydata files.")

    # Notify user if some clusters empty
    empty_count = 0
    for sz, fname in zip(cluster_sizes,cluster_fnames):
        if sz == 0:
            print(sz, ":", fname)
            empty_count += 1
    if empty_count:
        print(f"<{os.path.basename(__file__)}> Warning. Empty clusters found: {empty_count}")

    cluster_sizes = np.array(cluster_sizes)
    print(f"<{os.path.basename(__file__)}> Mean number of fibers per cluster: {np.mean(cluster_sizes)} Range: {np.min(cluster_sizes)}..{np.max(cluster_sizes)}")

    # Estimate subsampling ratio to display approximately number_of_fibers_to_display total fibers in 3D Slicer
    number_fibers = len(cluster_numbers_s)
    if number_fibers < number_of_fibers_to_display:
        ratio = 1.0
    else:
        ratio = number_of_fibers_to_display / number_fibers
    print(f"<{os.path.basename(__file__)}> Subsampling ratio for display of {number_of_fibers_to_display} total fibers estimated as: {ratio}")

    # Write the MRML file into the directory where the polydatas were already stored
    fname = os.path.join(outdir, 'clustered_tracts.mrml')
    mrml.write(fnames, np.around(np.array(cluster_colors), decimals=3), fname, ratio=ratio)

    # Also write one with 100% of fibers displayed
    fname = os.path.join(outdir, 'clustered_tracts_display_100_percent.mrml')
    mrml.write(fnames, np.around(np.array(cluster_colors), decimals=3), fname, ratio=1.0)
    
    # View the whole thing in jpg format for quality control
    if render_images:
        print(f'<{os.path.basename(__file__)}> Rendering and saving images of cluster atlas.')
        ren = render.render(output_polydata_s, 1000, data_mode='Cell', data_name='EmbeddingColor', verbose=verbose)
        ren.save_views(outdir, verbose=verbose)
        del ren

def mask_all_clusters(inpd, cluster_numbers_s, number_of_clusters, color=None, preserve_point_data=True, preserve_cell_data=False, verbose=True):

    print(f'<{os.path.basename(__file__)}> Masking all clusters: total {number_of_clusters}')

    inpoints = inpd.GetPoints()
    inpointdata = inpd.GetPointData()
    incelldata = inpd.GetCellData()

    # output and temporary objects
    outpd_list = []
    outlines_list = []
    outpoints_list = []

    for c_idx in range(0, number_of_clusters):
        outpd = vtk.vtkPolyData()
        outlines = vtk.vtkCellArray()
        outpoints = vtk.vtkPoints()

        outpd_list.append(outpd)
        outlines_list.append(outlines)
        outpoints_list.append(outpoints)

        del outpd
        del outlines
        del outpoints

    # check for cell data arrays to keep
    if preserve_cell_data:
        if incelldata.GetNumberOfArrays() > 0:
            cell_data_array_indices = list(range(incelldata.GetNumberOfArrays()))
            for idx in cell_data_array_indices:
                array = incelldata.GetArray(idx)
                dtype = array.GetDataType()

                for c_idx in range(0, number_of_clusters):
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

                    if verbose and c_idx == 0:
                        print(f"Cell data array found: {array.GetName()} {array.GetNumberOfComponents()}")

                    outpd_list[c_idx].GetCellData().AddArray(out_array)

        else:
            preserve_cell_data = False

    # check for point data arrays to keep
    tensor_names = []
    if preserve_point_data:
        if inpointdata.GetNumberOfArrays() > 0:
            point_data_array_indices = list(range(inpointdata.GetNumberOfArrays()))
            for idx in point_data_array_indices:
                array = inpointdata.GetArray(idx)
                for c_idx in range(0, number_of_clusters):
                    out_array = vtk.vtkFloatArray()
                    out_array.SetNumberOfComponents(array.GetNumberOfComponents())
                    out_array.SetName(array.GetName())

                    if verbose and c_idx == 0:
                        print(f"Point data array found: {array.GetName()} {array.GetNumberOfComponents()}")

                    outpd_list[c_idx].GetPointData().AddArray(out_array)

                if array.GetNumberOfComponents() == 9:
                    tensor_names.append(array.GetName())
        else:
            preserve_point_data = False

    # Set up scalars and tensors attributes for correct visualization in Slicer.
    # Slicer works with point data and does not handle cell data well.
    # This set of changes breaks old internal wma default visualization of cell scalars.
    # Changes must be propagated through wma so that render is called with the name of the field to visualize.
    # the new way in wma is like this line, below.
    # ren = wma.render.render(output_polydata_s, 1000, data_mode="Cell", data_name='EmbeddingColor')

    # For Slicer: First set one of the expected tensor arrays as default for vis
    tensors_labeled = False
    for name in tensor_names:
        if name == "tensors":
            for c_idx in range(0, number_of_clusters):
                outpd_list[c_idx].GetPointData().SetTensors(outpd_list[c_idx].GetPointData().GetArray("tensors"))
            tensors_labeled = True
        if name == "Tensors":
            for c_idx in range(0, number_of_clusters):
                outpd_list[c_idx].GetPointData().SetTensors(outpd_list[c_idx].GetPointData().GetArray("Tensors"))
            tensors_labeled = True
        if name == "tensor1":
            for c_idx in range(0, number_of_clusters):
                outpd_list[c_idx].GetPointData().SetTensors(outpd_list[c_idx].GetPointData().GetArray("tensor1"))
            tensors_labeled = True
        if name == "Tensor1":
            for c_idx in range(0, number_of_clusters):
                outpd_list[c_idx].GetPointData().SetTensors(outpd_list[c_idx].GetPointData().GetArray("Tensor1"))
            tensors_labeled = True

    if not tensors_labeled:
        if len(tensor_names) > 0:
            print(f"Data has unexpected tensor name(s). Unable to set active for visualization: {tensor_names}")

    # now set cell data visualization inactive.
    for c_idx in range(0, number_of_clusters):
        outpd_list[c_idx].GetCellData().SetActiveScalars(None)

    # loop over lines
    inpd.GetLines().InitTraversal()
    for c_idx in range(0, number_of_clusters):
        outlines_list[c_idx].InitTraversal()

    for lidx in range(0, inpd.GetNumberOfLines()):
        ptids = vtk.vtkIdList()
        inpd.GetLines().GetNextCell(ptids)

        if cluster_numbers_s[lidx] >= 0:
            c_idx = int(cluster_numbers_s[lidx])

            if verbose:
                if lidx % 100 == 0:
                    print(f"Line: {lidx} / {inpd.GetNumberOfLines()} belonging to cluster {c_idx}")

            # get points for each ptid and add to output polydata
            cellptids = vtk.vtkIdList()

            for pidx in range(0, ptids.GetNumberOfIds()):
                point = inpoints.GetPoint(ptids.GetId(pidx))
                idx_ = outpoints_list[c_idx].InsertNextPoint(point)
                cellptids.InsertNextId(idx_)

                if preserve_point_data:
                    for idx in point_data_array_indices:
                        array = inpointdata.GetArray(idx)
                        outpd_list[c_idx].GetPointData().GetArray(idx).InsertNextTuple(array.GetTuple(ptids.GetId(pidx)))

            outlines_list[c_idx].InsertNextCell(cellptids)

            if preserve_cell_data:
                for idx in cell_data_array_indices:
                    array = incelldata.GetArray(idx)
                    out_array = outpd_list[c_idx].GetCellData().GetArray(idx)
                    out_array.InsertNextTuple(array.GetTuple(lidx))

    for c_idx in range(0, number_of_clusters):
        # put data into output polydata
        outpd_list[c_idx].SetLines(outlines_list[c_idx])
        outpd_list[c_idx].SetPoints(outpoints_list[c_idx])

        if verbose:
            print(f"<{os.path.basename(__file__)}> Cluster {c_idx} have fibers {outpd_list[c_idx].GetNumberOfLines()} / {inpd.GetNumberOfLines()}, points {outpd_list[c_idx].GetNumberOfPoints()}")

    return outpd_list

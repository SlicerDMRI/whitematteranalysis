# goal: to look along tract in neighborhoods.
# maybe just like clustering then doing this though, without the cluster centroids...
# why would random centroids be better? idea is to approximate using all possible ones
# or uniform like voxels without hard boundaries. but...
# maybe need to compare centroids from clustering (with cluster defined boundaries)
# vs evenly spaced with hard or soft boundaries
# vs random.
# maybe this is not all that interesting though.
# it's interesting if the statistics define the segmentation rather than the clustering..... hmm.

# ------------------------------------------------------
# Can re-run the below after changing neighborhood_threshold
# Or after changing group membership
# All slow processing happens above this line
# ------------------------------------------------------
#neighborhood_threshold = 15.0
neighborhood_threshold = 20.0

group_indices =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,\
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

group_indices = numpy.array(group_indices)



from __future__ import print_function
print("Importing whitematteranalysis package.")
for m in ["io", "fibers", "filter", "laterality", "render", "cluster", "relative_distance", "mrml", "congeal_multisubject", "register_two_subjects", "register_two_subjects_nonrigid", "register_two_subjects_nonrigid_bsplines", "congeal_to_atlas", "tract_measurement"]:
  #print "importing ", m
  exec 'import %s'%m


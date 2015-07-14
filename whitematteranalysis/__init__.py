print "Importing whitematteranalysis package."
for m in ["io", "fibers", "filter", "laterality", "render", "cluster", "relative_distance", "mrml", "congeal_multisubject", "register_two_subjects", "register_two_subjects_nonlinear", "congeal_to_atlas"]:
  #print "importing ", m
  exec 'import %s'%m


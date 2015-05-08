print "Importing whitematteranalysis package."
for m in ["io", "fibers", "filter", "laterality", "render", "cluster", "relative_distance", "congeal", "registration_functions", "mrml", "congeal_multisubject", "register_two_subjects"]:
  #print "importing ", m
  exec 'import %s'%m


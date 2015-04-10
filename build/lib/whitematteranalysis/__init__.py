print "Unmodified wma from github"
for m in ["io", "fibers", "filter", "laterality", "render", "cluster", "relative_distance", "congeal", "registration_functions", "mrml"]:
  #print "importing ", m
  exec 'import %s'%m


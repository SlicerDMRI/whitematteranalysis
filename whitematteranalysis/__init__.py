#__all__ = ["io", "fibers", "laterality", "register", "render", "test"]
for m in ["io", "fibers", "filter", "laterality", "midsagalign", "register", "render", "cluster", "relative_distance", "congeal", "registration_functions"]:
  #print "importing ", m
  exec 'import %s'%m


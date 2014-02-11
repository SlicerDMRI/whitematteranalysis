import argparse
import os
import multiprocessing

try:
    import whitematteranalysis as wma
except:
    print "<wm_render_laterality_all.py> Error importing white matter analysis package\n"
    raise

try:
    from joblib import Parallel, delayed
except:
    print "<wm_render_laterality_all.py> Error importing joblib package\n"
    raise


fname_to_read = 'tractography_with_LI.vtk'
        
parser = argparse.ArgumentParser(
    description="Render LI polydatas",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu")

#-----------------
# Parse arguments
#-----------------
parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputDirectory',
    help='Contains laterality results, one directory per subject.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')


args = parser.parse_args()

if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

input_directories = os.listdir(args.inputDirectory)
print input_directories

print "wm_render_laterality_all. Starting white matter outlier computation."
print ""
print "=====render laterality======"

print "input directory:", args.inputDirectory
print 'CPUs detected:', multiprocessing.cpu_count()
if args.numberOfJobs is not None:
    parallel_jobs = args.numberOfJobs
else:
    parallel_jobs = multiprocessing.cpu_count()
print 'Using N jobs:', parallel_jobs

print "=========================="

# =======================================================================
# Above this line is argument parsing. Below this line is the pipeline.
# =======================================================================

def pipeline(indir):
    li_dir = os.path.join(args.inputDirectory,indir)
    print "<wm_render_laterality_all.py> processing directory:", li_dir
    if os.path.isdir(li_dir):
        pd = None
        fname = os.path.join(li_dir,fname_to_read)
        print fname
        if os.path.exists(fname):
            pd = wma.io.read_polydata(fname)
            if pd is not None:
                print "<wm_render_laterality_all.py> Rendering polydata:", fname
                ren = wma.render.render(pd, scalar_range=[-.5,.5], scalar_bar=True, number_of_fibers=1000)
                ren.save_views(li_dir)

# loop over all inputs
Parallel(n_jobs=parallel_jobs, verbose=0)(
        delayed(pipeline)(indir)
        for indir in input_directories)



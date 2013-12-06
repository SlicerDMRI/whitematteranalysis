import whitematteranalysis as wma
import argparse
import os

fname_to_read = 'tractography_with_LI.vtk'
        
parser = argparse.ArgumentParser(
    description="Render LI polydatas",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='Contains laterality results, one directory per subject.')


args = parser.parse_args()

if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

dirs = os.listdir(args.inputDirectory)
print dirs
current_dir = os.getcwd()

for indir in dirs:
    li_dir = os.path.join(args.inputDirectory,indir)
    print li_dir
    if os.path.isdir(li_dir):
        pd = None
        pd = wma.io.read_polydata(os.path.join(li_dir,fname_to_read))
        if pd is not None:
            ren = wma.render.render(pd, scalar_range=[-.5,.5], scalar_bar=True, number_of_fibers=1000)
            ren.save_views(li_dir)
            #os.chdir(current_dir)


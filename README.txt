===========
White Matter Analysis
===========

White Matter Analysis provides clustering and tractography analysis tools.

Installation
=========

1. Install python. You need scipy so: I used to use the Enthought version EPD since it has some pre-built packages needed by some of the code (scipy). However it now is packaged differently and I had trouble building VTK with python wrapping with EPD, on mac. Now I have installed Anaconda python. http://docs.continuum.io/anaconda/index.html 
Perhaps under linux scipy building is easier than on mac, so you may already have it.

2. Install required packages:

2.1 VTK: http://www.vtk.org/Wiki/VTK/Building
http://www.vtk.org/Wiki/VTK/Git/Download
You will need to compile it with python wrapping. VTK_WRAP_PYTHON must be on.
Make sure it finds the version of python that you want to use for this project, at configure time. You may need to toggle t for advanced mode in ccmake. I have something like this when I run:
cd VTK-build
ccmake ../VTK

PYTHON_EXECUTABLE                /Users/lauren/anaconda/bin/python            
 PYTHON_EXTRA_LIBS                                                             
 PYTHON_INCLUDE_DIR               /Users/lauren/anaconda/pkgs/python-2.7.4-1/in
 PYTHON_LIBRARY                   /Users/lauren/anaconda/lib/libpython2.7.dylib
 PYTHON_UTIL_LIBRARY              /usr/lib/libutil.dylib   

Note this requires git and cmake. More information is at vtk.org.
Then install it into your python:
cd VTK-build/Wrapping/Python
python setup.py install

2.2 Install the following python packages. This can be done with easy_install, for example.
joblib
scipy
(scipy is included in anaconda or EPD or other python distributions, possibly)

3. Install WhiteMatterAnalysis in your python

cd WhiteMatterAnalysis
python setup.py install

4. Test the installation
python
>>> import whitematteranalysis as wma

This will produce errors if the required packages are not found.


Test the installation
======================
Copy the script somewhere you can edit it:
WhiteMatterAnalysis/bin/cluster/cluster_atlas_and_label_subjects.py

edit the input and output directories at the top of the file:
indir = '/Users/odonnell/Data/tbi_with_scalars'
outdir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/March2013/TBI/cluster_5'

The input directory should hold several polydatas with whole-brain tractography, already registered. The output directory should be a different directory from the input.

Try to run the script:
python
>> exec file 'path_to_script/cluster_atlas_and_label_subjects.py'




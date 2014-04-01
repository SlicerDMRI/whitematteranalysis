whitematteranalysis
===================

#Synopsis
White Matter Analysis provides clustering and tractography analysis tools.

#Installation
##1. Install python. 
Anaconda is a nice option since it has VTK.
##2. Install required packages:

###2.1 If your python does not already have VTK, install it.
VTK: http://www.vtk.org/Wiki/VTK/Building
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

###2.2 Install the following python packages (dependencies). 
This can be done with easy_install, for example.
joblib
scipy
(scipy is included in anaconda or EPD or other python distributions, possibly)

##3. Install WhiteMatterAnalysis in your python

cd WhiteMatterAnalysis
python setup.py install

##4. Test the installation
python
>>> import whitematteranalysis as wma

This will produce errors if the required packages are not found.
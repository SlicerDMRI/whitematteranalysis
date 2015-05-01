whitematteranalysis
===================

#Synopsis
White Matter Analysis provides clustering and tractography analysis tools.

It implements algorithms from publications listed here:
http://projects.iq.harvard.edu/whitematteranalysis/publications

Also see the github.io page here:
http://ljod.github.io/whitematteranalysis/

#Installation
###1. Install python. 
Anaconda is a nice option since it has VTK and scipy.
First install anaconda, then run: conda install vtk

###2. Install the following python packages (dependencies).

Once you have anaconda installed, run: pip install joblib

####If you decide to use another python that does not already have VTK, compile it.
* VTK: http://www.vtk.org/Wiki/VTK/Building
* http://www.vtk.org/Wiki/VTK/Git/Download

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

###3. Install WhiteMatterAnalysis into your python.

     cd whiteMatterAnalysis
     python setup.py install

###4. Test the installation.
Run python, then:

    >>> import whitematteranalysis as wma

This will produce errors if the required packages are not found.

###5. Please see the wiki for usage instructions.
https://github.com/ljod/whitematteranalysis/wiki

whitematteranalysis
===================

#Synopsis
White Matter Analysis provides clustering and tractography analysis tools.

It implements algorithms from publications listed here:
http://projects.iq.harvard.edu/whitematteranalysis/publications

Also see the github.io page here:
http://slicerdmri.github.io/whitematteranalysis/

Please reference the following papers:

O'Donnell, LJ., and Westin, CF. Automatic
tractography segmentation using a high-dimensional white matter
atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.

O’Donnell LJ, Wells III WM, Golby AJ, Westin CF. Unbiased groupwise registration of white matter tractography. InInternational Conference on Medical Image Computing and Computer-Assisted Intervention 2012 Oct 1 (pp. 123-130). Springer Berlin Heidelberg.

For projects using Slicer please include this text (or similar):

"We performed tractography visualization with anatomical hierarchies in 3D Slicer (http://www.slicer.org) via the SlicerDMRI project (https://github.com/SlicerDMRI), funded by NIH U01 CA199459."

#Installation
###1. Download whitematteranalysis from github. 

      git clone https://github.com/SlicerDMRI/whitematteranalysis.git
      
###2. Install python. 
Anaconda is a nice option since it has VTK and scipy.
First install anaconda from http://continuum.io/downloads, then run: 

      conda install vtk

###3. Install the following python packages (dependencies).

Once you have anaconda installed, run: 

      pip install joblib

Note: If you decide to use another python that does not already have VTK, you can compile VTK.
* VTK: http://www.vtk.org/Wiki/VTK/Building
* http://www.vtk.org/Wiki/VTK/Git/Download

You will need to compile it with python wrapping. VTK_WRAP_PYTHON must be on.
Make sure that at configure time it finds the version of python that you want to use for this project. You may need to toggle t for advanced mode in ccmake. I have something like this when I run:
     cd VTK-build
     ccmake ../VTK

       PYTHON_EXECUTABLE                /Users/lauren/anaconda/bin/python            
       PYTHON_EXTRA_LIBS                                                             
       PYTHON_INCLUDE_DIR               /Users/lauren/anaconda/pkgs/python-2.7.4-1/in
       PYTHON_LIBRARY                   /Users/lauren/anaconda/lib/libpython2.7.dylib
       PYTHON_UTIL_LIBRARY              /usr/lib/libutil.dylib   

Note this requires both git and cmake. More information is at vtk.org.
To install your compiled vtk into your python:
     cd VTK-build/Wrapping/Python
     python setup.py install

###4. Install WhiteMatterAnalysis into your python in the standard way.

     cd whitematteranalysis
     python setup.py install

###7. Please see the wiki for usage instructions.
https://github.com/SlicerDMRI/whitematteranalysis/wiki

whitematteranalysis
===================

# Synopsis
White Matter Analysis provides clustering and tractography analysis tools.

It implements algorithms from publications listed here:
http://projects.iq.harvard.edu/whitematteranalysis/publications

Also see the github.io page here:
http://slicerdmri.github.io/whitematteranalysis/

**Please cite the following papers:**

    O'Donnell, LJ., and Westin, CF. Automatic tractography segmentation
    using a high-dimensional white matter atlas. Medical Imaging,
    IEEE Transactions on 26.11 (2007): 1562-1575.

    O’Donnell LJ, Wells III WM, Golby AJ, Westin CF. Unbiased groupwise registration of white matter tractography.
    In International Conference on Medical Image Computing and Computer-Assisted Intervention 2012 Oct 1 (pp. 123-130).
    Springer Berlin Heidelberg.

**For projects using Slicer and SlicerDMRI please also include the following text (or similar) and citations:**

* How to cite the [Slicer platform](http://wiki.slicer.org/slicerWiki/index.php/CitingSlicer)
* An example of how to cite SlicerDMRI (modify the first part of the sentence according to your use case):

    "We performed diffusion MRI tractography and/or analysis and/or visualization in 3D Slicer (www.slicer.org) via the SlicerDMRI project (dmri.slicer.org) (Norton et al. 2017)."
    
    - [Isaiah Norton, Walid Ibn Essayed, Fan Zhang, Sonia Pujol, Alex Yarmarkovich, Alexandra J. Golby, Gordon Kindlmann, Demian Wassermann, Raul San Jose Estepar, Yogesh Rathi, Steve Pieper, Ron Kikinis, Hans J. Johnson, Carl-Fredrik Westin and Lauren J. O'Donnell. SlicerDMRI: Open Source Diffusion MRI Software for Brain Cancer Research. Cancer Research 77(21), e101-e103, 2017.](http://cancerres.aacrjournals.org/content/77/21/e101)

# Installation
## 1. Install python 2. 
Anaconda is a nice option since it has VTK and scipy, and includes up-to-date pip and setuptools.

  - Download and install Miniconda (**Python 2**) from https://conda.io/miniconda.html

### 1a. (Windows-only): install VTK 6.3 from [C. Gohlke's windows packages](https://www.lfd.uci.edu/~gohlke/pythonlibs/#vtk) (PyPI VTK packages for Windows are only available for Python 3):

      pip install https://download.lfd.uci.edu/pythonlibs/o4uhg4xd/VTK-6.3.0-cp27-cp27m-win_amd64.whl
   
## 2. Install whitematteranalysis with pip: 

      pip install git+https://github.com/SlicerDMRI/whitematteranalysis.git
      


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

## Documentation
* Please see the wiki for usage instructions of whitematteranalysis.

    https://github.com/SlicerDMRI/whitematteranalysis/wiki

* Please see the following page for instructions of applying a pre-provided anatomically curated white matter atlas, along with the computation tools provided in whitematteranalysis, to perform subject-specific tractography parcellation. 

    https://dmri.slicer.org/atlases



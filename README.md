whitematteranalysis
===================

# Synopsis
White Matter Analysis provides clustering and tractography analysis tools.

It implements algorithms from publications listed here:
http://projects.iq.harvard.edu/whitematteranalysis/publications

Also see the github.io page here:
http://slicerdmri.github.io/whitematteranalysis/

# Installation
## 1. Install Python 2. 
Miniconda is a nice option since it includes up-to-date pip and setuptools and all library dependencies (such as VTK and scipy).

  - Download and install Miniconda (**Python 2**) from https://conda.io/miniconda.html

## 2. Install whitematteranalysis with pip:

The following command will use pip to install whitematteranalysis from this source repository and all library dependencies:

      pip install git+https://github.com/SlicerDMRI/whitematteranalysis.git

  (Note: On MacOS, to able to use pip, X-code needs to be installed using `xcode-select --install`.)

Run 'wm_quality_control_tractography.py --help' to test if the installation is successful. 

# Documentation
* Please see the wiki for usage instructions of whitematteranalysis.

    https://github.com/SlicerDMRI/whitematteranalysis/wiki

* Please see the following page for instructions of applying a pre-provided anatomically curated white matter atlas, along with the computation tools provided in whitematteranalysis, to perform subject-specific tractography parcellation. 

    https://dmri.slicer.org/atlases

# References

**Please cite the following papers:**

    O'Donnell, LJ., and Westin, CF. Automatic tractography segmentation
    using a high-dimensional white matter atlas. Medical Imaging,
    IEEE Transactions on 26.11 (2007): 1562-1575.

    O’Donnell LJ, Wells III WM, Golby AJ, Westin CF. Unbiased groupwise registration of white matter tractography.
    In International Conference on Medical Image Computing and Computer-Assisted Intervention 2012 Oct 1 (pp. 123-130).
    Springer Berlin Heidelberg.

    Zhang, F., Wu, Y., Norton, I., Rathi, Y., Makris, N., O'Donnell, LJ. An anatomically curated fiber clustering white matter atlas for consistent white matter tract parcellation across the lifespan. NeuroImage, 2018 (179): 429-447

**For projects using Slicer and SlicerDMRI please also include the following text (or similar) and citations:**

* How to cite the [Slicer platform](http://wiki.slicer.org/slicerWiki/index.php/CitingSlicer)
* An example of how to cite SlicerDMRI (modify the first part of the sentence according to your use case):

    "We performed diffusion MRI tractography and/or analysis and/or visualization in 3D Slicer (www.slicer.org) via the SlicerDMRI project (dmri.slicer.org) (Norton et al. 2017)."
    
    - [Isaiah Norton, Walid Ibn Essayed, Fan Zhang, Sonia Pujol, Alex Yarmarkovich, Alexandra J. Golby, Gordon Kindlmann, Demian Wassermann, Raul San Jose Estepar, Yogesh Rathi, Steve Pieper, Ron Kikinis, Hans J. Johnson, Carl-Fredrik Westin and Lauren J. O'Donnell. SlicerDMRI: Open Source Diffusion MRI Software for Brain Cancer Research. Cancer Research 77(21), e101-e103, 2017.](http://cancerres.aacrjournals.org/content/77/21/e101)

whitematteranalysis
===================

# Synopsis

WhiteMatterAnalysis (WMA) provides fiber clustering and tractography analysis tools.

# WMA Installation
## 1. Install Python 3:

Miniconda is a nice option becuase it includes pip, setuptools, and all library dependencies (such as VTK and scipy).

  - Download and install Miniconda (**Python 3.7**) from https://conda.io/miniconda.html

## 2. Install whitematteranalysis with pip:

The following command will use pip to install whitematteranalysis from this source repository and all library dependencies:

      pip install git+https://github.com/SlicerDMRI/whitematteranalysis.git


  (Note: On MacOS, to able to use pip, X-code needs to be installed using `xcode-select --install`.)

Run `wm_quality_control_tractography.py --help` to test if the installation is successful.

# Documentation

* Please see the following page for instructions of applying a pre-provided anatomically curated white matter atlas, along with the computation tools provided in whitematteranalysis, to perform subject-specific tractography parcellation.

    https://github.com/SlicerDMRI/whitematteranalysis/blob/master/doc/subject-specific-tractography-parcellation.md

* WMA implements algorithms from publications listed here:
http://projects.iq.harvard.edu/whitematteranalysis/publications

# References

**Please cite the following papers for using the code and/or the training data :**
 
    Zhang, F., Wu, Y., Norton, I., Rathi, Y., Makris, N., O'Donnell, LJ. 
    An anatomically curated fiber clustering white matter atlas for consistent white matter tract parcellation across the lifespan. 
    NeuroImage, 2018 (179): 429-447

    O'Donnell LJ, Wells III WM, Golby AJ, Westin CF. 
    Unbiased groupwise registration of white matter tractography.
    In MICCAI, 2012, pp. 123-130.

    O'Donnell, LJ., and Westin, CF. Automatic tractography segmentation
    using a high-dimensional white matter atlas. Medical Imaging,
    IEEE Transactions on 26.11 (2007): 1562-1575.

**For projects using Slicer and SlicerDMRI please also include the following text (or similar) and citations:**

* How to cite the [Slicer platform](http://wiki.slicer.org/slicerWiki/index.php/CitingSlicer)
* An example of how to cite SlicerDMRI (modify the first part of the sentence according to your use case):

    "We performed diffusion MRI tractography and/or analysis and/or visualization in 3D Slicer (www.slicer.org) via the SlicerDMRI project (dmri.slicer.org) (Norton et al. 2017)."
    
    Fan Zhang, Thomas Noh, Parikshit Juvekar, Sarah F Frisken, Laura Rigolo, Isaiah Norton, Tina Kapur, Sonia Pujol, William Wells III, Alex Yarmarkovich, Gordon Kindlmann, Demian Wassermann, Raul San Jose Estepar, Yogesh Rathi, Ron Kikinis, Hans J Johnson, Carl-Fredrik Westin, Steve Pieper, Alexandra J Golby, Lauren J O’Donnell. 
    SlicerDMRI: Diffusion MRI and Tractography Research Software for Brain Cancer Surgery Planning and Visualization. 
    JCO Clinical Cancer Informatics 4, e299-309, 2020.

    Isaiah Norton, Walid Ibn Essayed, Fan Zhang, Sonia Pujol, Alex Yarmarkovich, Alexandra J. Golby, Gordon Kindlmann, Demian Wassermann, Raul San Jose Estepar, Yogesh Rathi, Steve Pieper, Ron Kikinis, Hans J. Johnson, Carl-Fredrik Westin and Lauren J. O'Donnell. 
    SlicerDMRI: Open Source Diffusion MRI Software for Brain Cancer Research. Cancer Research 77(21), e101-e103, 2017.

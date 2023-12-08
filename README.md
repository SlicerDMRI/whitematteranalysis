whitematteranalysis
===================

[![test, package](https://github.com/SlicerDMRI/whitematteranalysis/actions/workflows/test_package.yaml/badge.svg?branch=master)](https://github.com/SlicerDMRI/whitematteranalysis/actions/workflows/test_package.yaml?query=branch%3Amaster)
[![Documentation Status](https://readthedocs.org/projects/whitematteranalysis/badge/?version=latest)](https://whitematteranalysis.readthedocs.io/en/latest/?badge=latest)
[![code format](https://github.com/SlicerDMRI/whitematteranalysis/actions/workflows/check_format.yaml/badge.svg?branch=master)](https://github.com/SlicerDMRI/whitematteranalysis/actions/workflows/check_format.yaml?query=branch%3Amaster)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)


# Synopsis

WhiteMatterAnalysis (WMA) provides fiber clustering and tractography analysis tools.

# WMA Installation and Usage
## 1. Install Python 3:

Miniconda is a nice option because it includes pip, setuptools, and all library dependencies (such as VTK and scipy).

  - Download and install Miniconda from https://conda.io/miniconda.html

## 2. Install whitematteranalysis with pip:

The following command will use pip to install whitematteranalysis from this source repository and all library dependencies:

```shell
$ pip install git+https://github.com/SlicerDMRI/whitematteranalysis.git
```

Note: On macOS, to be able to use `pip`, `X-code` needs to be installed using `$ xcode-select --install`.

Run `$ wm_quality_control_tractography.py --help` to test if the installation is successful.

## 3. Documentation

The `whitematteranalysis` package documentation can be found at
https://whitematteranalysis.readthedocs.io/en/latest/.

* A master shell script `wm_apply_ORG_atlas_to_subject.sh` (see code [here](https://github.com/SlicerDMRI/whitematteranalysis/blob/73a7948751f49d9fda8ec84fb5caeecaeeb92621/bin/wm_apply_ORG_atlas_to_subject.sh#L1-L40)) is provided to apply an anatomically curated white matter atlas ([the ORG atlas](https://dmri.slicer.org/atlases/)), along with the computation tools provided in whitematteranalysis, to perform subject-specific tractography parcellation.

```shell
$ wm_apply_ORG_atlas_to_subject.sh \
  -i input_tractography.vtk \
  -o output_dir \
  -a path_to_atlas/ORG-Atlases-v1.x \
  -s /Applications/Slicer5.2.2.app/Contents/MacOS/Slicer \
  -d 1 \
  -m /Applications/Slicer5.2.2.app/Contents/Extensions-31382/SlicerDMRI/lib/Slicer-5.2/cli-modules/FiberTractMeasurements
```

* A step-by-step tutorial to explain the detailed computational process within the pipeline is provided [here](https://whitematteranalysis.readthedocs.io/en/latest/subject-specific-tractography-parcellation.html).

* The names of the anatomical bundles WMA extracts can be found in the [bundles](https://whitematteranalysis.readthedocs.io/en/latest/bundles.html) page.

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

* How to cite the [Slicer platform](https://slicer.readthedocs.io/en/latest/user_guide/about.html#how-to-cite)
* An example of how to cite SlicerDMRI (modify the first part of the sentence according to your use case):

    "We performed diffusion MRI tractography and/or analysis and/or visualization in 3D Slicer (https://www.slicer.org/) via the SlicerDMRI project (https://dmri.slicer.org/) (Norton et al. 2017)."
    
    Fan Zhang, Thomas Noh, Parikshit Juvekar, Sarah F Frisken, Laura Rigolo, Isaiah Norton, Tina Kapur, Sonia Pujol, William Wells III, Alex Yarmarkovich, Gordon Kindlmann, Demian Wassermann, Raul San Jose Estepar, Yogesh Rathi, Ron Kikinis, Hans J Johnson, Carl-Fredrik Westin, Steve Pieper, Alexandra J Golby, Lauren J O’Donnell. 
    SlicerDMRI: Diffusion MRI and Tractography Research Software for Brain Cancer Surgery Planning and Visualization. 
    JCO Clinical Cancer Informatics 4, e299-309, 2020.

    Isaiah Norton, Walid Ibn Essayed, Fan Zhang, Sonia Pujol, Alex Yarmarkovich, Alexandra J. Golby, Gordon Kindlmann, Demian Wassermann, Raul San Jose Estepar, Yogesh Rathi, Steve Pieper, Ron Kikinis, Hans J. Johnson, Carl-Fredrik Westin and Lauren J. O'Donnell. 
    SlicerDMRI: Open Source Diffusion MRI Software for Brain Cancer Research. Cancer Research 77(21), e101-e103, 2017.

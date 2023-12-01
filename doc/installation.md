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

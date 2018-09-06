import os, glob
from distutils.core import setup
from setuptools import setup, Extension, Command
from setuptools.command.build_ext import build_ext as _build_ext

#    url='http://pypi.python.org/pypi/WhiteMatterAnalysis/',
#    scripts=['bin/test1.py','bin/test2.py'],

# ext_modules work-around: https://stackoverflow.com/a/38057196
# get_include work-around: https://stackoverflow.com/a/21621689

import sys
if sys.platform == 'win32':
    # force setuptools to use the MSVC compiler environment
    # otherwise it will check version for VC9, and error out
    os.environ['MSSdk'] = '1'
    os.environ['DISTUTILS_USE_SDK'] = '1'


class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())


setup(
    name='WhiteMatterAnalysis',
    version='0.2.0',
    author='Lauren O\'Donnell',
    author_email='odonnell@bwh.harvard.edu',
    packages=['whitematteranalysis'],
    license='LICENSE.txt',
    description='Processing of whole-brain streamline tractography.',
    long_description=open('README.md').read(),

    install_requires = ['setuptools>=18.0', 'numpy', 'scipy', 'cython', 'joblib', 'statsmodels', 'xlrd'],    
    
    ext_modules = [
        Extension('whitematteranalysis.fibers', sources=['whitematteranalysis/fibers.pyx']),
        Extension('whitematteranalysis.similarity', sources=['whitematteranalysis/similarity.pyx']),
        ],
    cmdclass={'build_ext':build_ext},

    scripts = [ 
        'bin/picktracts_converter.py',
        'bin/harden_transform_with_slicer.py',
        'bin/wm_cluster_atlas.py',
        'bin/wm_cluster_from_atlas.py',
        'bin/wm_extract_cluster.py',
        'bin/wm_laterality_all.py',
        'bin/wm_preprocess_all.py',
        'bin/wm_quality_control_tractography.py',
        'bin/wm_register_multisubject_faster.py',
        'bin/wm_register_to_atlas_new.py',
        'bin/wm_separate_clusters_by_hemisphere.py',
        'bin/wm_flatten_length_distribution.py',
        'bin/wm_create_mrml_file.py',
        'bin/wm_cluster_remove_outliers.py',
        'bin/wm_measure_all_subjects.py',
        'bin/wm_quality_control_cluster_measurements.py',
        'bin/wm_statistics.py',
        'bin/wm_statistics_export_data.py',
        'bin/wm_quality_control_after_clustering.py',
        'bin/wm_harden_transform.py',
        'bin/wm_append_clusters.py',
        'bin/wm_remove_data_along_tracts.py',
        'bin/wm_measure_endpoint_overlap.py',
        'bin/wm_extract_clusters_by_endpoints.py',
        'bin/wm_assess_cluster_location.py',
        'bin/wm_vtp2vtk.py',
        'bin/wm_download_anatomically_curated_atlas.py'
    ]
)

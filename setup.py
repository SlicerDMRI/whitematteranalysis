import os, glob
from distutils.core import setup

import setuptools
from setuptools import setup, Extension, Command
from setuptools.command.build_ext import build_ext as _build_ext
from pkg_resources import resource_filename

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

########################################################################
# Begin attribution section
########################################################################
#
# License: BSD
#       Created: August 16, 2012
#       Author:  Francesc Alted - francesc@blosc.org
#
########################################################################
class LazyCommandClass(dict):
    """
    Lazy command class that defers operations requiring Cython and numpy until
    they've actually been downloaded and installed by setup_requires.
    """
    def __contains__(self, key):
        return (
            key == 'build_ext'
            or super().__contains__(key)
        )

    def __setitem__(self, key, value):
        if key == 'build_ext':
            raise AssertionError("build_ext overridden!")
        super().__setitem__(key, value)

    def __getitem__(self, key):
        if key != 'build_ext':
            return super().__getitem__(key)

        from Cython.Distutils import build_ext as cython_build_ext

        class build_ext(cython_build_ext):
            """
            Custom build_ext command that lazily adds numpy's include_dir to
            extensions.
            """
            def build_extensions(self):
                """
                Lazily append numpy's include directory to Extension includes.
                This is done here rather than at module scope because setup.py
                may be run before numpy has been installed, in which case
                importing numpy and calling `numpy.get_include()` will fail.
                """
                numpy_incl = resource_filename('numpy', 'core/include')
                for ext in self.extensions:
                    ext.include_dirs.append(numpy_incl)

                # This explicitly calls the superclass method rather than the
                # usual super() invocation because distutils' build_class, of
                # which Cython's build_ext is a subclass, is an old-style class
                # in Python 2, which doesn't support `super`.
                cython_build_ext.build_extensions(self)
        return build_ext
########################################################################
# End attribution section
########################################################################


setup_requires = ['cython', 'numpy']
setup(
    name='WhiteMatterAnalysis',
    version='0.3.0',
    author='Fan Zhang and Lauren O\'Donnell',
    author_email='fzhang@bwh.harvard.edu; odonnell@bwh.harvard.edu',
    packages=['whitematteranalysis'],
    license='LICENSE.txt',
    description='Processing of whole-brain streamline tractography.',
    long_description=open('README.md').read(),
  
    setup_requires = setup_requires,
    install_requires = setup_requires + ['setuptools', 'scipy', 'vtk',
                        'joblib', 'statsmodels', 'xlrd', 'matplotlib', 'nibabel'],
    
    ext_modules = [
        Extension('whitematteranalysis.fibers', sources=['whitematteranalysis/fibers.pyx']),
        Extension('whitematteranalysis.similarity', sources=['whitematteranalysis/similarity.pyx']),
        ],
    cmdclass=LazyCommandClass(),

    scripts = [ 
        'bin/harden_transform_with_slicer.py',
        'bin/wm_append_clusters.py',
        'bin/wm_append_clusters_to_anatomical_tracts.py',
        'bin/wm_apply_ORG_atlas_to_subject.sh',
        'bin/wm_assess_cluster_location_by_hemisphere.py',
        'bin/wm_cluster_volumetric_measurements.py',
        'bin/wm_cluster_atlas.py',
        'bin/wm_cluster_from_atlas.py',
        'bin/wm_cluster_remove_outliers.py',
        'bin/wm_compare_vtks.py',
        'bin/wm_create_mrml_file.py',
        'bin/wm_diffusion_measurements.py',
        'bin/wm_download_anatomically_curated_atlas.py',
        'bin/wm_harden_transform.py',
        'bin/wm_preprocess_all.py',
        'bin/wm_quality_control_tractography.py',
        'bin/wm_quality_control_after_clustering.py',
        'bin/wm_quality_control_cluster_measurements.py',
        'bin/wm_quality_control_tract_overlap.py',
        'bin/wm_register_multisubject_faster.py',
        'bin/wm_register_to_atlas_new.py',
        'bin/wm_remove_data_along_tracts.py',
        'bin/wm_separate_clusters_by_hemisphere.py',
        'bin/wm_tract_to_volume.py',
        'bin/wm_vtp2vtk.py',
        'bin/wm_change_nrrd_dir.py',
        'testing/test_run.py'
    ]
)

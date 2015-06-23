from distutils.core import setup
#    url='http://pypi.python.org/pypi/WhiteMatterAnalysis/',
#    scripts=['bin/test1.py','bin/test2.py'],
from Cython.Build import cythonize
import numpy

setup(
    name='WhiteMatterAnalysis',
    version='0.1.0',
    author='Lauren O\'Donnell',
    author_email='odonnell@bwh.harvard.edu',
    packages=['whitematteranalysis'],
    license='LICENSE.txt',
    description='Processing of whole-brain streamline tractography.',
    long_description=open('README.md').read(),
    ext_modules = cythonize("whitematteranalysis/*.pyx"),
    scripts = [ 
        'bin/convert_ascii_vtk_to_binary.py',
        'bin/picktracts_converter.py',
        'bin/wm_append_all.py',
        'bin/wm_cluster_atlas.py',
        'bin/wm_cluster_from_atlas.py',
        'bin/wm_cluster_subject.py',
        'bin/wm_extract_cluster.py',
        'bin/wm_laterality_all.py',
        'bin/wm_mask_lines.py',
        'bin/wm_measure_all_clusters.py',
        'bin/wm_preprocess_all.py',
        'bin/wm_quality_control_tractography.py',
        'bin/wm_register_multisubject_faster.py',
        'bin/wm_register_multisubject.py',
        'bin/wm_register_to_atlas_new.py',
        'bin/wm_register_to_atlas.py',
        'bin/wm_separate_clusters_by_hemisphere.py',
        'bin/wm_flatten_length_distribution.py'
    ],
    include_dirs=[numpy.get_include()]
)


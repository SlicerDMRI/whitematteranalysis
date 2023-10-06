import os
from setuptools import setup
from pkg_resources import resource_filename

#    url='http://pypi.python.org/pypi/whitematteranalysis/',
#    scripts=['bin/test1.py','bin/test2.py'],

# ext_modules work-around: https://stackoverflow.com/a/38057196
# get_include work-around: https://stackoverflow.com/a/21621689

import sys
if sys.platform == 'win32':
    # force setuptools to use the MSVC compiler environment
    # otherwise it will check version for VC9, and error out
    os.environ['MSSdk'] = '1'
    os.environ['DISTUTILS_USE_SDK'] = '1'


with open("requirements.txt") as f:
    required_dependencies = f.read().splitlines()
    external_dependencies = []
    for dependency in required_dependencies:
        if dependency[0:2] == "-e":
            repo_name = dependency.split("=")[-1]
            repo_url = dependency[3:]
            external_dependencies.append("f{repo_name} @ {repo_url}")
        else:
            external_dependencies.append(dependency)


setup_requires = ['numpy>=1.20.0']
setup(
    name='whitematteranalysis',
    version='0.3.0',
    author='Fan Zhang and Lauren O\'Donnell',
    author_email='fzhang@bwh.harvard.edu; odonnell@bwh.harvard.edu',
    packages=['whitematteranalysis'],
    license='LICENSE.txt',
    description='Processing of whole-brain streamline tractography.',
    long_description=open('README.md').read(),
    setup_requires=setup_requires,
    install_requires=setup_requires + external_dependencies,
    scripts=[
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

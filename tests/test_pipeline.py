#!/usr/bin/env python

from __future__ import print_function

import sys, os
from subprocess import check_output

test_out_dir = sys.argv[1]

if not os.path.isdir(test_out_dir):
    os.makedirs(test_out_dir)

print("#### Step 1. tractography quality control:")

check_output([
    'wm_quality_control_tractography.py', 'test_data_small', os.path.join(test_out_dir, 'S1_QC')
    ])

print("#### Step 2. groupwise tract registration:")

check_output([
    'wm_register_multisubject_faster.py', '-l', '20', '-f', '200', '-j', '2',
    '-midsag_symmetric', '-mode', 'affine', 'test_data_small',
    os.path.join(test_out_dir, 'S2_GroupRegistration/')
    ])

print("#### Step 3. atlas generation:")

check_output(['wm_cluster_atlas.py', '-l', '20', '-f', '200',
    '-nystrom_sample', '250', '-k', '10', '-j', '2',
    os.path.join(test_out_dir, 'S2_GroupRegistration/output_tractography/'),
    os.path.join(test_out_dir, 'S3_Atlas')
    ])

print("#### Step 4. subject registration:")

check_output(['wm_register_to_atlas_new.py', '-l', '40', '-mode', 'affine',
    os.path.join(test_data_small, 'brain_0001.vtk'),
    os.path.join(test_out_dir, 'S2_GroupRegistration/registration_atlas.vtk'),
    os.path.join(test_out_dir, 'S4_RegisteredSubject')
    ])

print("#### Step 5. subject clustering:")

check_output(['wm_cluster_from_atlas.py', '-j', '4', '-l', '40',
    os.path.join(test_out_dir, 'S4_RegisteredSubject/brain_0001/output_tractography/brain_0001_reg.vtk'),
    os.path.join(test_out_dir, 'S3_Atlas/iteration_00002/initial_clusters/'),
    os.path.join(test_out_dir, 'S5_ClusteredSubject')
    ])

print("#### Step 6. subject outlier removal:")

check_output([
    'wm_cluster_remove_outliers.py',
    os.path.join(test_out_dir, 'S5_ClusteredSubject/brain_0001_reg'),
    os.path.join(test_out_dir, 'S3_Atlas/iteration_00002/initial_clusters/'),
    os.path.join(test_out_dir, 'S6_OutlierRemovedSubject')
    ])

print("#### Step 7. subject hemisphere separation:")

check_output([
    'wm_separate_clusters_by_hemisphere.py', '-pthresh', '0.6', '-atlasMRML',
    os.path.join(test_out_dir, 'S3_Atlas/iteration_00002/initial_clusters/clustered_tracts_display_100_percent.mrml'),
    os.path.join(test_out_dir, 'S6_OutlierRemovedSubject/brain_0001_reg_outlier_removed'),
    os.path.join(test_out_dir, 'S7_HemisphereSeparatedSubject')
    ])
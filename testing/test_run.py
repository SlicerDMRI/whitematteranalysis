#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import wm_quality_control_tractography
import wm_register_to_atlas_new
import wm_cluster_from_atlas
import wm_cluster_remove_outliers
import wm_separate_clusters_by_hemisphere
import numpy as np
import random

np.random.seed(0)
random.seed(0)

sys.argv = ['wm_quality_control_tractography.py', 'input_tractography/', 'qc_output/']
wm_quality_control_tractography.main()

sys.argv =['wm_register_to_atlas_new.py', '-l', '40', '-mode', 'affine', 'input_tractography/bundle.vtk', 'ORG-RegAtlas-100HCP/registration_atlas.vtk', 'registered_subject_output/']
wm_register_to_atlas_new.main()

sys.argv =['wm_cluster_from_atlas.py', '-l', '40', 'registered_subject_output/bundle/output_tractography/bundle_reg.vtk', 'ORG-800FC-100HCP/', 'subject_cluster_output/']
wm_cluster_from_atlas.main()

sys.argv =['wm_cluster_remove_outliers.py', 'subject_cluster_output/bundle_reg', 'ORG-800FC-100HCP/', 'subject_cluster_outlier_removed_output/']
wm_cluster_remove_outliers.main()

sys.argv = ['wm_separate_clusters_by_hemisphere.py', 'subject_cluster_outlier_removed_output/', 'subject_cluster_separated_output/']
wm_separate_clusters_by_hemisphere.main()

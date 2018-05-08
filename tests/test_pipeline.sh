#set -x # show commands as-run
#set -e # fail on error

test_out_dir=$1 # need to specify an output folder

mkdir ${test_out_dir}

echo "#### Step 1. tractography quality control:"

wm_quality_control_tractography.py test_data_small ${test_out_dir}/S1_QC

echo "#### Step 2. groupwise tract registration:"

wm_register_multisubject_faster.py -l 20 -f 200 -j 2 -midsag_symmetric -mode affine test_data_small/ ${test_out_dir}/S2_GroupRegistration/

echo "#### Step 3. atlas generation:"

wm_cluster_atlas.py -l 20 -f 200 -nystrom_sample 250 -k 10 -j 2 ${test_out_dir}/S2_GroupRegistration/output_tractography/ ${test_out_dir}/S3_Atlas

echo "#### Step 4. subject registration:"

wm_register_to_atlas_new.py -l 40 -mode affine test_data_small/brain_0001.vtk ${test_out_dir}/S2_GroupRegistration/registration_atlas.vtk ${test_out_dir}/S4_RegisteredSubject

echo "#### Step 5. subject clustering:"

wm_cluster_from_atlas.py -j 4 -l 40  ${test_out_dir}/S4_RegisteredSubject/brain_0001/output_tractography/brain_0001_reg.vtk ${test_out_dir}/S3_Atlas/iteration_00002/initial_clusters/ ${test_out_dir}/S5_ClusteredSubject

echo "#### Step 6. subject outlier removal:"

wm_cluster_remove_outliers.py ${test_out_dir}/S5_ClusteredSubject/brain_0001_reg ${test_out_dir}/S3_Atlas/iteration_00002/initial_clusters/ ${test_out_dir}/S6_OutlierRemovedSubject

echo "#### Step 7. subject hemisphere separation:"

wm_separate_clusters_by_hemisphere.py -pthresh 0.6 -atlasMRML ${test_out_dir}/S3_Atlas/iteration_00002/initial_clusters/clustered_tracts_display_100_percent.mrml ${test_out_dir}/S6_OutlierRemovedSubject/brain_0001_reg_outlier_removed ${test_out_dir}/S7_HemisphereSeparatedSubject

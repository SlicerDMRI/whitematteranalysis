
echo "<test> Running clustering test:"

# turn off outlier rejection and positive definite approx for consistency with test results from before
python -m cProfile -o profile.txt ../bin/wm_cluster_atlas.py test_data_small output/test_cluster_atlas -nystrom_sample 100 -advanced_only_testing_on -advanced_only_random_seed 10 -cluster_outlier_std 1000 -iter 1 -advanced_only_force_pos_def_off -l 0 -f 500

echo "<test> cProfile run statistics:"
python view_profile_stats.py  profile.txt


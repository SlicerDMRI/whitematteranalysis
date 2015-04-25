
echo "<test> Running clustering test:"


python -m cProfile -o profile.txt ../bin/wm_cluster_atlas.py test_data_small output/test_cluster_atlas -nystrom_sample 100 -advanced_only_testing_on -advanced_only_random_seed 10

echo "<test> cProfile run statistics:"
python view_profile_stats.py


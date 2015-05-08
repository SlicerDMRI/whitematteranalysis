echo "<test> Testing registration functionality."

echo "<test> MULTISUBJECT REGISTRATION."

echo "<test> python -m cProfile -o profile.txt ../bin/wm_register_multisubject_faster.py test_data_small output/test_reg_multi_faster"
python -m cProfile -o profile.txt ../bin/wm_register_multisubject_faster.py test_data_small output/test_reg_multi_faster  -norender

echo "<test> cProfile run statistics:"
python view_profile_stats.py profile.txt 


#echo "<test> REGISTRATION TO ATLAS."
#python -m cProfile -o profile.txt ../bin/wm_register_multisubject.py test_data_2 output/test_reg_atlas



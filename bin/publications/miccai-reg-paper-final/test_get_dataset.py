# commands to copy the dataset from pnl to local disk brage at LMI

cmd = 'scp -rp ringo:///projects/schiz/3Tdata/{0}/diff/*dwi-filt-Ed.* .'

# fill in the subject ids
subjects = ['case01193', 'case01194', 'case01202', 'case01215', 'case01231', 'case01234', 'case01242', 'case01243', 'case01244', 'case01247', 'case01249', 'case01254', 'case01261', 'case01262', 'case01263', 'case01264', 'case01271', 'case01272', 'case01276', 'case01277', 'case01278', 'case01280', 'case01281', 'case01282', 'case01288', 'case01302', 'case01307', 'case01309', 'case01316', 'case01324', 'case01327', 'case01328', 'case01331', 'case01337', 'case01339', 'case01342', 'case01345', 'case01347', 'case01349', 'case01353', 'case01355', 'case01356', 'case01357', 'case01358', 'case01361', 'case01363', 'case01374', 'case01376', 'case01383', 'case01387', 'case01395']

subjects2 = '{case01193,case01194,case01202,case01215,case01231,case01234,case01242,case01243,case01244,case01247,case01249,case01254,case01261,case01262,case01263,case01264,case01271,case01272,case01276,case01277,case01278,case01280,case01281,case01282,case01288,case01302,case01307,case01309,case01316,case01324,case01327,case01328,case01331,case01337,case01339,case01342,case01345,case01347,case01349,case01353,case01355,case01356,case01357,case01358,case01361,case01363,case01374,case01376,case01383,case01387,case01395}'

for subject in subjects:
    print cmd.format(subject)
    
print cmd.format(subjects2)

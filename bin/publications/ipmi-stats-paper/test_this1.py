if USE_PARALLEL:
    retval = Parallel(n_jobs=n_jobs, verbose=1)(
        delayed(do_permutation)(
            pidx, number_of_permutations, group_indices_real, number_of_subjects, number_of_fibers_per_subject, distances, distances_sq, neighborhood_sigma, statistic_threshold, distance_threshold, kl_threshold, bc_threshold
            ) for pidx in range(number_of_permutations))


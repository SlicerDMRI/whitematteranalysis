#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_separate_clusters_by_hemisphere.py", "--help"])
    assert ret.success


# def test_execution(script_runner):
    # in_dirname = "subject_cluster_outlier_removed_output"
    # out_dirname = "subject_cluster_separated_output"
    # ret = script_runner.run(
    #    "wm_separate_clusters_by_hemisphere.py",
    #    in_dirname,
    #    out_dirname,
    # )
    # assert ret.success

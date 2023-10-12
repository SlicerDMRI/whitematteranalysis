#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_cluster_remove_outliers.py", "--help"])
    assert ret.success


# def test_execution(script_runner):
    # in_tractography_fname = "subject_cluster_output/bundle_reg"
    # in_atlas_fname = "ORG-800FC-100HCP"
    # out_dirname = "ORG-800FC-100HCP"
    # ret = script_runner.run(
    #    "wm_cluster_remove_outliers.py",
    #    in_tractography_fname,
    #    in_atlas_fname,
    #    out_dirname,
    # )
    # assert ret.success

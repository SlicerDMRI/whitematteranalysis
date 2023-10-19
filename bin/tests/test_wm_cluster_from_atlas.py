#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_cluster_from_atlas.py", "--help"])
    assert ret.success


# def test_execution(script_runner):
    # length = 40
    # in_tractography_fname = (
    #    "registered_subject_output/bundle/output_tractography/bundle_reg.vtk"
    # )
    # in_atlas_dirname = "ORG-800FC-100HCP"
    # out_dirname = "subject_cluster_output"
    # ret = script_runner.run(
    #    "wm_cluster_from_atlas.py",
    #    "-l",
    #    length,
    #    in_tractography_fname,
    #    in_atlas_dirname,
    #    out_dirname,
    # )
    # assert ret.success

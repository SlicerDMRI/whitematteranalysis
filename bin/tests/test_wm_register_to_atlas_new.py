#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_register_to_atlas_new.py", "--help"])
    assert ret.success


# def test_execution(script_runner):
    # length = 40
    # mode = "affine"
    # in_tractography_fname = "input_tractography/bundle.vtk"
    # in_atlas_fname = "ORG-RegAtlas-100HCP/registration_atlas.vtk"
    # out_dirname = "registered_subject_output"
    # ret = script_runner.run(
    #     "wm_register_to_atlas_new.py",
    #     "-l",
    #     length,
    #     "-mode",
    #     mode,
    #     in_tractography_fname,
    #     in_atlas_fname,
    #     out_dirname,
    # )
    # assert ret.success

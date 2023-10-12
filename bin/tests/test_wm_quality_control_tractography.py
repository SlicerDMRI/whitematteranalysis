#!/usr/bin/env python
# -*- coding: utf-8 -*-


def test_help_option(script_runner):
    ret = script_runner.run(["wm_quality_control_tractography.py", "--help"])
    assert ret.success


# def test_execution(script_runner):
    # in_dirname = "input_tractography"
    # out_dirname = "qc_output"
    # ret = script_runner.run(
    #    "wm_quality_control_tractography.py", in_dirname, out_dirname,
    # )
    # assert ret.success

#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(
        ["wm_append_diffusion_measures_across_subjects.py", "--help"])
    assert ret.success

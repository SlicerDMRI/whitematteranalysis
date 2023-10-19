#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(
        ["wm_compute_bundle_feature_population_math.py", "--help"])
    assert ret.success

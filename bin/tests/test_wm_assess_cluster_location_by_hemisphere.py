#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(
        ["wm_assess_cluster_location_by_hemisphere.py", "--help"])
    assert ret.success

#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_quality_control_after_clustering.py", "--help"])
    assert ret.success

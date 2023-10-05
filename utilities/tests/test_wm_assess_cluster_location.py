#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_assess_cluster_location.py", "--help"])
    assert ret.success

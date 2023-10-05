#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_measure_all_clusters.py", "--help"])
    assert ret.success

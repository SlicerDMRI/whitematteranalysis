#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(
        ["utilities/wm_measure_endpoint_overlap.py", "--help"])
    assert ret.success

#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_fix_UKF_trace.py", "--help"])
    assert ret.success

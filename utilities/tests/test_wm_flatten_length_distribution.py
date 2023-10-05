#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_flatten_length_distribution.py", "--help"])
    assert ret.success

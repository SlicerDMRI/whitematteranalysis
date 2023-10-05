#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_register_multisubject_faster.py", "--help"])
    assert ret.success

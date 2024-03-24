#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["bin/wm_preprocess_all.py", "--help"])
    assert ret.success

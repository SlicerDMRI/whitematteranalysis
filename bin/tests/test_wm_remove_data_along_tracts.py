#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_remove_data_along_tracts.py", "--help"])
    assert ret.success

#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(
        ["wm_fix_hemisphere_loc_name_in_tractogram.py", "--help"])
    assert ret.success

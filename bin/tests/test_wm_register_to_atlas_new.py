#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_register_to_atlas_new.py", "--help"])
    assert ret.success

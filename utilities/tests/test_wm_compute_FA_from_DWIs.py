#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_compute_FA_from_DWIs.py", "--help"])
    assert ret.success

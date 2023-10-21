#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_query_atlas_bundles_names.py", "--help"])
    assert ret.success

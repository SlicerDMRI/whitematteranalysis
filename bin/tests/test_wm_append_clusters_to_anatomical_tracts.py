#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(
        ["wm_append_clusters_to_anatomical_tracts.py", "--help"])
    assert ret.success

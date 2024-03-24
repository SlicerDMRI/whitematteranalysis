#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["bin/harden_transform_with_slicer.py", "--help"])
    assert ret.success

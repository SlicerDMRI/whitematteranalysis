#!/usr/bin/env python
# -*- coding: utf-8 -*-

def test_help_option(script_runner):
    ret = script_runner.run(["wm_diffusion_measurements.py", "--help"])
    assert ret.success

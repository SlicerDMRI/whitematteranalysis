#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

from whitematteranalysis.utils.tripwire import (TripWire, TripWireError,
                                                is_tripwire)


def test_is_tripwire():
    assert not is_tripwire(object())
    assert is_tripwire(TripWire("some message"))


def test_tripwire():
    # Test tripwire object
    not_a_package_name = TripWire("Do not have not_a_package")
    unittest.TestCase.assertRaises(
        TripWireError, getattr, not_a_package_name, "do_something")
    unittest.TestCase.assertRaises(TripWireError, not_a_package_name)
    # Check AttributeError can be checked too
    try:
        not_a_package_name.__wrapped__
    except TripWireError as err:
        assert_true(isinstance(err, AttributeError))
    else:
        raise RuntimeError("No error raised, but expected.")

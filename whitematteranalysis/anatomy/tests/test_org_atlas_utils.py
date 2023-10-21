#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd

from whitematteranalysis.anatomy.org_atlas_utils import (
    ORGAtlasBundleFileHeading, get_association_bundles, get_bundle_long_name,
    get_bundle_short_name, get_cerebellar_bundles,
    get_commissural_augmented_bundles, get_commissural_bundles,
    get_hemispheric_bundles, get_hemispheric_mono_bundles,
    get_projection_bundles, get_striatal_bundles, get_superficial_bundles)
from whitematteranalysis.data.atlas.utils import (ORGAtlasVersion,
                                                  get_local_atlas_bundle_fname)


def _check_bundle_count(method_name, expected_bndl_count):

    version = ORGAtlasVersion.__members__.values()

    for ver, expected_count in zip(version, expected_bndl_count):
        fname = get_local_atlas_bundle_fname(ver)
        df = pd.read_csv(fname, sep=",")
        cer = method_name(df)
        obtained_count = len(cer)

        assert obtained_count == expected_count


def _check_bundle_name(method_name, version, bundle_id, expected_name):

    fname = get_local_atlas_bundle_fname(version)
    df = pd.read_csv(fname, sep=",")

    obtained_name = method_name(df, bundle_id)

    assert obtained_name == expected_name


def test_wma_get_association_bundles():

    expected_bndl_count = [20, 20, 20, 20, 20, 20]
    _check_bundle_count(get_association_bundles, expected_bndl_count)


def test_wma_get_bundle_long_name():

    version = ORGAtlasVersion.V1_2

    bundle_id = [201]
    method_name = get_bundle_long_name
    expected_name = ["cortico - ponto - cerebellar left"]
    _check_bundle_name(method_name, version, bundle_id, expected_name)

    fname = get_local_atlas_bundle_fname(version)
    df = pd.read_csv(fname, sep=",")
    stl = get_striatal_bundles(df)
    bundle_id = stl[ORGAtlasBundleFileHeading.ID.value].values

    expected_name = ["external capsule left", "external capsule right"]
    _check_bundle_name(method_name, version, bundle_id, expected_name)


def test_wma_get_bundle_short_name():

    version = ORGAtlasVersion.V1_2

    bundle_id = [201]
    method_name = get_bundle_short_name
    expected_name = ["CPC_left"]
    _check_bundle_name(method_name, version, bundle_id, expected_name)

    fname = get_local_atlas_bundle_fname(version)
    df = pd.read_csv(fname, sep=",")
    stl = get_striatal_bundles(df)
    bundle_id = stl[ORGAtlasBundleFileHeading.ID.value].values

    expected_name = ["EC_left", "EC_right"]
    _check_bundle_name(method_name, version, bundle_id, expected_name)


def test_wma_get_cerebellar_bundles():

    expected_bndl_count = [9, 9, 9, 9, 9, 9, 9]
    _check_bundle_count(get_cerebellar_bundles, expected_bndl_count)


def test_wma_get_commissural_bundles():

    expected_bndl_count = [7, 7, 7, 7, 7, 7]
    _check_bundle_count(get_commissural_bundles, expected_bndl_count)


def test_wma_get_commissural_augmented_bundles():

    expected_count = 8

    version = ORGAtlasVersion.V1_2
    fname = get_local_atlas_bundle_fname(version)
    df = pd.read_csv(fname, sep=",")

    com_aug = get_commissural_augmented_bundles(df)
    obtained_count = len(com_aug)

    assert obtained_count == expected_count


def test_wma_get_hemispheric_bundles():

    expected_count = 66

    version = ORGAtlasVersion.V1_2
    fname = get_local_atlas_bundle_fname(version)
    df = pd.read_csv(fname, sep=",")

    com_aug = get_hemispheric_bundles(df)
    obtained_count = len(com_aug)

    assert obtained_count == expected_count


def test_wma_get_hemispheric_mono_bundles():

    expected_count = 33

    version = ORGAtlasVersion.V1_2
    fname = get_local_atlas_bundle_fname(version)
    df = pd.read_csv(fname, sep=",")

    com_aug = get_hemispheric_mono_bundles(df)
    obtained_count = len(com_aug)

    assert obtained_count == expected_count


def test_wma_get_striatal_bundles():

    expected_bndl_count = [2, 2, 2, 2, 2, 2]
    _check_bundle_count(get_striatal_bundles, expected_bndl_count)


def test_wma_get_projection_bundles():

    expected_bndl_count = [20, 20, 20, 22, 20, 22]
    _check_bundle_count(get_projection_bundles, expected_bndl_count)


def test_wma_get_superficial_bundles():

    expected_bndl_count = [16, 16, 16, 16, 16, 16]
    _check_bundle_count(get_superficial_bundles, expected_bndl_count)

# -*- coding: utf-8 -*-

import enum

import pandas as pd

from whitematteranalysis.fileio.utils import (SlicerFileExtension,
                                              VTKFileExtension, find_filenames)

atlas_bundle_file_prefix = "T_"
nonatlas_bundle_names = ["T_FalsePositive"]


class Hemisphere(enum.Enum):
    LEFT = "left"
    RIGHT = "right"


class AtlasAvailability(enum.Enum):
    DEFAULT = "default"
    SEPARATED = "separated"


class ORGAtlasBundleType(enum.Enum):
    ASS = "Association"
    CER = "Cerebellar"
    COM = "Commissural"
    PRO = "Projection"
    STR = "Striatal"
    SUP = "Superficial"


class ORGAtlasBundleFileHeading(enum.Enum):
    ID = "id"
    SHORT_NAME = "short_name"
    LONG_NAME = "long_name"


class ORGAtlasBundleRange(enum.Enum):
    ASO = (101, 199)
    CER = (201, 299)
    COM = (301, 399)
    PRO = (401, 499)
    STR = (501, 599)
    SUP = (601, 699)

    @staticmethod
    def get_range_low(name):
        return ORGAtlasBundleRange(name).value[0]

    @staticmethod
    def get_range_high(name):
        return ORGAtlasBundleRange(name).value[1]


def get_bundles(df, bundle_id):
    return df[df[ORGAtlasBundleFileHeading.ID.value].isin(bundle_id)]


def get_bundle_feature(df, bundle_id, column):
    return get_bundles(df, bundle_id)[column].values.tolist()


def get_bundle_long_name(df, bundle_id):
    return get_bundle_feature(df, bundle_id, ORGAtlasBundleFileHeading.LONG_NAME.value)


def get_bundle_short_name(df, bundle_id):
    return get_bundle_feature(df, bundle_id, ORGAtlasBundleFileHeading.SHORT_NAME.value)


def get_bundle_type(df, bundle_type):
    return df[
        df[
            ORGAtlasBundleFileHeading.ID.value].between(
            ORGAtlasBundleRange.get_range_low(bundle_type),
            ORGAtlasBundleRange.get_range_high(bundle_type)
        )
    ]


def get_association_bundles(df):
    return get_bundle_type(df, ORGAtlasBundleRange.ASO)


def get_cerebellar_bundles(df):
    return get_bundle_type(df, ORGAtlasBundleRange.CER)


def get_commissural_bundles(df):
    return get_bundle_type(df, ORGAtlasBundleRange.COM)


def get_commissural_augmented_bundles(df):
    com = get_commissural_bundles(df)
    # Add the MCP
    mcp = get_bundles(df, [209])

    return pd.concat([com, mcp], ignore_index=True)


def get_hemispheric_bundles(df):
    aso = get_association_bundles(df)
    cer = get_cerebellar_bundles(df)
    # Remove the MCP
    cer = cer.drop(cer[cer[ORGAtlasBundleFileHeading.ID.value] == 209].index)
    stl = get_striatal_bundles(df)
    pro = get_projection_bundles(df)
    sup = get_superficial_bundles(df)

    return pd.concat([aso, cer, stl, pro, sup], ignore_index=True)


def get_hemispheric_mono_bundles(df):
    hem_mono = get_hemispheric_bundles(df)

    # Remove all right values
    hem_mono = hem_mono.drop(hem_mono.loc[hem_mono[ORGAtlasBundleFileHeading.SHORT_NAME.value].str.contains(Hemisphere.RIGHT.value)].index)

    # Drop the left substring
    short_name_substr = "_" + Hemisphere.LEFT.value
    column = ORGAtlasBundleFileHeading.SHORT_NAME.value
    hem_mono[column] = hem_mono[column].str.replace(short_name_substr, "")
    long_name_substr = " " + Hemisphere.LEFT.value
    column = ORGAtlasBundleFileHeading.LONG_NAME.value
    hem_mono[column] = hem_mono[column].str.replace(long_name_substr, "")

    # Strip the ID as it makes no longer sense
    hem_mono.drop(ORGAtlasBundleFileHeading.ID.value, axis=1, inplace=True)

    return hem_mono


def get_projection_bundles(df):
    return get_bundle_type(df, ORGAtlasBundleRange.PRO)


def get_striatal_bundles(df):
    return get_bundle_type(df, ORGAtlasBundleRange.STR)


def get_superficial_bundles(df):
    return get_bundle_type(df, ORGAtlasBundleRange.SUP)


def query_bundle_names_from_anatomical_tracts(path):
    bundle_names = find_filenames(path, VTKFileExtension.VTP, stem=True)
    return bundle_names


def query_bundle_names_from_scene_files(path):
    bundle_names = find_filenames(path, SlicerFileExtension.MRML, stem=True)

    # Remove non-atlas filenames if present
    [bundle_names.remove(elem) for elem in nonatlas_bundle_names]
    bundle_names = [elem for elem in bundle_names if elem.startswith(atlas_bundle_file_prefix)]

    return bundle_names


def add_org_atlas_prefix(bundles):
    return [atlas_bundle_file_prefix + elem for elem in bundles]


def query_bundle_names(path, atlas_availability):

    if atlas_availability == AtlasAvailability.DEFAULT:
        return query_bundle_names_from_scene_files(path)
    elif atlas_availability == AtlasAvailability.SEPARATED:
        return query_bundle_names_from_anatomical_tracts(path)
    else:
        raise ValueError(
            f"Unsupported value:\n"
            f"Found: {atlas_availability}; Available: {AtlasAvailability.__members__}")

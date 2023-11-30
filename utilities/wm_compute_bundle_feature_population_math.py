#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Compute bundle feature population math: computes the
- difference
- stats
- sum

for each bundle feature across participants according to the values contained in
the provided CSV files.

For the difference case, the first CSV file given is taken as the reference; for
the sum case, only the streamline count addition is performed.
"""

import argparse
import enum
import os
from pathlib import Path

import pandas as pd

relative_label = "relative"
sum_label = "sum"
underscore = "_"

# ToDo
# Think about the effects of dealing with NANs in the subtractions, stats, etc.


class BundlePopulationMathOperation(enum.Enum):
    DIFFERENCE = "difference"  # for each subject, with respect to a ref
    STATS = "stats"  # across all subjects
    SUM = "sum"  # for each subject, only the streamline count


class WMABundleFeatureDataMap(enum.Enum):
    BundleName = ("Name", str)
    PointCount = ("Num_Points", int)
    StreamlineCount = ("Num_Fibers", int)
    MeanLength = ("Mean_Length", float)
    EstimatedUncertaintyMean = ("EstimatedUncertainty.Mean", float)
    FA1Mean = ("FA1.Mean", float)
    FA2Mean = ("FA2.Mean", float)
    FWMean = ("FreeWater.Mean", float)
    # ToDo
    # HemisphereLocataion.Mean can contain integers, and cannot be cast directly
    # into integers; would need to do .astype(float).astype("Int64")
    HemisphereLocationMean = ("HemisphereLocation.Mean", float)
    ClusterIdxMean = ("cluster_idx.Mean", float)
    Trace1Mean = ("trace1.Mean", float)
    Trace2Mean = ("trace2.Mean", str)

    @staticmethod
    def get_name(name):
        return WMABundleFeatureDataMap(name).value[0]

    @staticmethod
    def get_type(name):
        return WMABundleFeatureDataMap(name).value[1]

    @staticmethod
    def get_type_map():
        type_map = dict(
            map(
                lambda x: x.value, WMABundleFeatureDataMap._member_map_.values()
            )
        )
        # Drop the name
        type_map.pop(
            WMABundleFeatureDataMap.get_name(WMABundleFeatureDataMap.BundleName)
        )

        return type_map


def compute_bundle_population_feature_diff(df_list):
    return [df_list[0].subtract(elem, fill_value=0) for elem in df_list[1:]]


def compute_bundle_population_feature_diff_relative(df_diff, def_ref):
    # rel_change = (new_value – ref_value) / ref_value * 100
    return [df.divide(def_ref).multiply(100) for df in df_diff]


def compute_bundle_population_feature_sum(df_list):
    column_name = WMABundleFeatureDataMap.get_name(WMABundleFeatureDataMap.StreamlineCount)
    return [pd.DataFrame([df[column_name].sum()], columns=[column_name]) for df in df_list]


def compute_bundle_population_feature_stats(df_list):
    df = pd.concat(df_list)
    column_name = WMABundleFeatureDataMap.get_name(WMABundleFeatureDataMap.BundleName)
    return df.groupby(column_name).describe()


def perform_bundle_population_operation(operation, df_list):

    if operation == BundlePopulationMathOperation.DIFFERENCE:
        return compute_bundle_population_feature_diff(df_list)
    elif operation == BundlePopulationMathOperation.STATS:
        return compute_bundle_population_feature_stats(df_list)
    elif operation == BundlePopulationMathOperation.SUM:
        return compute_bundle_population_feature_sum(df_list)
    else:
        raise NotImplementedError(
           f"Unsupported operation:\nFound: {operation}\n"
           f"Available: {list(BundlePopulationMathOperation.__members__)}"
        )


def cast_feature_data(df):

    # Cast all columns to the appropriate types
    type_map = WMABundleFeatureDataMap.get_type_map()
    df = df.astype(type_map)

    return df


def clean_up_feature_data(df):

    # Remove whitespaces from column names
    df.columns = df.columns.str.rstrip(" ")
    df.columns = df.columns.str.lstrip(" ")

    # Strip the path from the "Name" column
    column_name = WMABundleFeatureDataMap.get_name(WMABundleFeatureDataMap.BundleName)
    df[column_name] = pd.Series([str(Path(path)).replace(str(Path(path).parent) + os.sep, "") for path in df[column_name]])
    df[column_name] = pd.Series([str(Path(path)).replace(str(Path(path).suffix), "") for path in df[column_name]])
    return df.set_index(column_name)


def process_feature_data(df):

    df = clean_up_feature_data(df)
    return cast_feature_data(df)


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "operation",
        help="Population feature math operation.",
        type=BundlePopulationMathOperation,
        choices=list(BundlePopulationMathOperation),
    )
    parser.add_argument(
        "bundle_names_fname", help="Bundle names filename (*.txt).", type=Path
    )
    parser.add_argument(
        "--out_fnames",
        nargs="+",
        help="Output TSV filenames (*.tsv).",
        type=Path,
    )
    parser.add_argument(
        "--in_feature_fnames", nargs="+", help="Input feature filenames (*.csv).", type=Path
    )

    return parser


def _parse_args(parser):

    args = parser.parse_args()

    return args


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    with open(args.bundle_names_fname) as file:
        bndl_names = [line.rstrip() for line in file]

    # for bndl_name in bndl_names:

    df_list = []
    for fname in args.in_feature_fnames:

        df = pd.read_csv(fname)

        # Prepare the data
        df = process_feature_data(df)

        # ToDo
        # Keep only the requested bundles

        df_list.append(df)

    df = perform_bundle_population_operation(args.operation, df_list)

    sep = "\t"

    # Compute the relative difference if applicable
    if args.operation == BundlePopulationMathOperation.DIFFERENCE:
        df_rel = compute_bundle_population_feature_diff_relative(df, df_list[0])

        for _df, _df_rel, fname in zip(df, df_rel, args.out_fnames):
            path = Path(fname).parent
            stem = Path(fname).stem + underscore + relative_label
            ext = Path(fname).suffix
            _fname = Path(path, stem).with_suffix(ext)
            _df_rel.to_csv(_fname, sep=sep)

            _df.to_csv(fname, sep=sep)

    elif args.operation == BundlePopulationMathOperation.STATS:
        df.to_csv(args.out_fnames[0], sep=sep)
    elif args.operation == BundlePopulationMathOperation.SUM:
        for _df, fname in zip(df, args.out_fnames):
            _df.to_csv(fname, sep=sep, index=False)

    else:
        raise NotImplementedError(
           f"Unsupported operation:\nFound: {args.operation}\n"
           f"Available: {list(BundlePopulationMathOperation.__members__)}"
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Query the atlas bundle names.

e.g. python wm_query_atlas_bundles.py \
    /mnt/data/atlas/ORG-Atlas-1.4/ORG-800FC-100HCP/ \
    /mnt/data/atlas/ORG-Atlas-1.4/bundle_names_default.txt \
    DEFAULT

 python wm_query_atlas_bundles.py \
    /mnt/data/atlas/ORG-Atlas-1.4/ORG-800FC-100HCP-separated/AnatomicalTracts/ \
    /mnt/data/atlas/ORG-Atlas-1.4/bundle_names_separated.txt \
    SEPARATED
"""

import argparse
from pathlib import Path

from whitematteranalysis.anatomy.org_atlas_utils import (AtlasAvailability,
                                                         query_bundle_names)
from whitematteranalysis.fileio.utils import save2txt


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("in_path", type=Path, help="Input path.")
    parser.add_argument(
        "out_fname",
        help="Output filename (*.txt).",
        type=Path,
    )
    parser.add_argument(
        "atlas_availability",
        choices=AtlasAvailability._member_names_,
        help="Atlas availability.",
    )

    return parser


def _parse_args(parser):

    args = parser.parse_args()

    return args


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    bundle_names = query_bundle_names(
        args.in_path,
        AtlasAvailability.__getitem__(args.atlas_availability),
    )
    save2txt(args.out_fname, bundle_names)


if __name__ == "__main__":
    main()

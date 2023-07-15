# -*- coding: utf-8 -*-

import enum
import json
import pathlib
import sys

if sys.version_info < (3, 10):
    from importlib_resources import files
else:
    from importlib.resources import files


class ORGAtlasVersion(enum.Enum):
    V1_1 = "V1_1"
    V1_1_1 = "V1_1_1"
    V1_2 = "V1_2"
    V1_3_A = "V1_3_A"
    V1_3_B = "V1_3_B"
    V1_4 = "V1_4"


def get_local_atlas_bundle_fname(_atlas):

    path = files(__package__).joinpath("org_atlas_version.json")

    with path.open() as f:
        org_atlas_version = json.load(f)

    return pathlib.Path(str(path)).parent.joinpath(org_atlas_version[_atlas.value])

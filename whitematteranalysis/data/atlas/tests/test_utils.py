#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from whitematteranalysis.data.atlas.utils import (ORGAtlasVersion,
                                                  get_local_atlas_bundle_fname)


def test_get_local_atlas_bundle_fname():

    version = ORGAtlasVersion.__members__.values()
    assert [
        os.path.exists(get_local_atlas_bundle_fname(ver)) for ver in version
    ]

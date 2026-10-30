"""
This package provides support for sky background subtraction and equalization
(matching), including a HEALPix-based multi-exposure offset solver.
"""

import logging

from .healpix_match import (
    match_healpix_skyvals,
    pipeline_level_from_epsilon,
)
from .offset_io import (
    load_offset_table,
    lookup_offset,
    normalize_exposure_key,
    resolve_offsets_path_from_asn,
    write_offset_table,
    write_skyval_parquet,
)
from .skymatch_step import SkyMatchStep

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = [
    "SkyMatchStep",
    "load_offset_table",
    "lookup_offset",
    "match_healpix_skyvals",
    "normalize_exposure_key",
    "pipeline_level_from_epsilon",
    "resolve_offsets_path_from_asn",
    "write_offset_table",
    "write_skyval_parquet",
]

"""
Parquet I/O for HEALPix skyval samples and skymatch offset tables.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from urllib.parse import urlparse

import numpy as np
from astropy.table import Table

log = logging.getLogger(__name__)

__all__ = [
    "ASN_OFFSETS_KEY",
    "load_offset_table",
    "load_skyval_parquet",
    "lookup_offset",
    "normalize_exposure_key",
    "offsets_path_from_asn",
    "write_offset_table",
    "write_skyval_parquet",
]

# Top-level ASN key for the offsets parquet path (one file per L3 ASN).
ASN_OFFSETS_KEY = "skymatch_offsets_file"

# Strip a single known product suffix + extension for offset-table keys.
_SUFFIX_RE = re.compile(
    r"(?i)(_segm|_cal|_skymatch|_cat|_i2d|_crf)?(\.asdf|\.parquet|\.h5|\.fits)?$"
)


def normalize_exposure_key(name: str | Path) -> str:
    """
    Normalize an exposure filename/stem for offset table lookup.

    Strips the final path component's extension and a single trailing Roman
    product suffix so ``..._cal.asdf``, ``..._segm.asdf``, and bare stems map
    to the same key.
    """
    stem = Path(str(name)).name
    stem = _SUFFIX_RE.sub("", stem)
    if "." in stem:
        stem = Path(stem).stem
    return stem


def write_skyval_parquet(
    outfile: str | Path,
    *,
    pixel: np.ndarray,
    value: np.ndarray,
    unc: np.ndarray,
    coverage: np.ndarray | None = None,
    filename: str | None = None,
    conversion_megajanskys: float | None = None,
) -> Path:
    """
    Write one exposure's HEALPix skyval samples to a parquet file.

    Columns: ``healpix17``, ``data``, ``err`` (DN/s). Optional constant
    ``conversion_megajanskys`` (DN/s → MJy/sr) is stored as a column.
    Coarse coverage is stored in a sidecar ``*_cov.parquet`` when provided.

    ``outfile`` is written exactly as given (caller supplies ``.parquet``).
    """
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    n = len(np.asarray(pixel).reshape(-1))
    table = Table(
        {
            "healpix17": np.asarray(pixel, dtype=np.int64).reshape(-1),
            "data": np.asarray(value, dtype=np.float32).reshape(-1),
            "err": np.asarray(unc, dtype=np.float32).reshape(-1),
        }
    )
    if filename is not None:
        table["filename"] = np.full(n, normalize_exposure_key(filename))
    if conversion_megajanskys is not None:
        table["conversion_megajanskys"] = np.full(
            n, float(conversion_megajanskys), dtype=np.float64
        )
    table.write(outfile, format="parquet", overwrite=True)

    if coverage is not None:
        cov_path = outfile.with_name(outfile.stem + "_cov.parquet")
        cov_table = Table(
            {"healpix11": np.asarray(coverage, dtype=np.int64).reshape(-1)}
        )
        cov_table.write(cov_path, format="parquet", overwrite=True)
    return outfile


def load_skyval_parquet(path: str | Path) -> dict[str, np.ndarray | float]:
    """
    Load skyval samples from a parquet file.

    Requires columns ``healpix17``, ``data``, ``err``.
    """
    path_str = str(path)
    table = Table.read(path_str, format="parquet")
    for col in ("healpix17", "data", "err"):
        if col not in table.colnames:
            raise ValueError(f"{path_str} missing required column {col!r}")

    sample: dict[str, np.ndarray | float] = {
        "pixel": np.asarray(table["healpix17"], dtype=np.int64),
        "value": np.asarray(table["data"], dtype=np.float64),
        "unc": np.asarray(table["err"], dtype=np.float64),
    }
    if "conversion_megajanskys" in table.colnames:
        conv = np.asarray(table["conversion_megajanskys"], dtype=np.float64)
        sample["conversion_megajanskys"] = float(conv.flat[0])

    cov_path = Path(path_str).with_name(Path(path_str).stem + "_cov.parquet")
    if not urlparse(path_str).scheme and cov_path.exists():
        cov = Table.read(str(cov_path), format="parquet")
        sample["coverage"] = np.asarray(cov["healpix11"], dtype=np.int64)
    return sample


def write_offset_table(
    outfile: str | Path,
    filenames: Sequence[str],
    levels: Sequence[float],
) -> Path:
    """
    Write a skymatch offsets parquet table.

    Columns: ``filename``, ``level`` (values for ``meta.background.level``).
    ``outfile`` is written exactly as given.
    """
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    keys = [normalize_exposure_key(f) for f in filenames]
    level_arr = np.asarray(levels, dtype=np.float64)
    if level_arr.shape[0] != len(keys):
        raise ValueError("filenames and levels must have the same length")

    table = Table({"filename": keys, "level": level_arr})
    table.write(outfile, format="parquet", overwrite=True)
    log.info("Wrote skymatch levels for %d exposures to %s", len(keys), outfile)
    return outfile


def load_offset_table(path: str | Path) -> dict[str, float]:
    """
    Load offsets parquet into ``{normalized_filename: level}``.

    Requires columns ``filename`` and ``level``.
    """
    path_str = str(path)
    table = Table.read(path_str, format="parquet")
    if "filename" not in table.colnames:
        raise ValueError(f"{path_str} missing filename column")
    if "level" not in table.colnames:
        raise ValueError(f"{path_str} missing level column")

    out: dict[str, float] = {}
    for name, level in zip(table["filename"], table["level"], strict=True):
        out[normalize_exposure_key(str(name))] = float(level)
    return out


def lookup_offset(filename: str | Path, table: Mapping[str, float]) -> float:
    """Look up pipeline level for an exposure name; raise if missing."""
    key = normalize_exposure_key(filename)
    if key not in table:
        raise KeyError(
            f"No skymatch offset for {filename!r} (key={key!r}); have {len(table)} entries"
        )
    return table[key]


def offsets_path_from_asn(asn: Mapping) -> str | None:
    """
    Return the top-level ``skymatch_offsets_file`` string from an ASN, if set.

    Path joining / URL handling is left to ``ModelLibrary._member_to_filename``
    (same rules as association science members).
    """
    raw = asn.get(ASN_OFFSETS_KEY)
    if not raw:
        return None
    return str(raw)

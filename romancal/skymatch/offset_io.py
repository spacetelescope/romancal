"""
Parquet I/O for HEALPix skyval samples and skymatch offset tables.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from astropy.table import Table

log = logging.getLogger(__name__)

__all__ = [
    "ASN_OFFSETS_KEY",
    "load_offset_table",
    "load_skyval_parquet",
    "lookup_offset",
    "normalize_exposure_key",
    "resolve_offsets_path_from_asn",
    "write_offset_table",
    "write_skyval_parquet",
]

# Product-level ASN key for the offsets parquet path.
ASN_OFFSETS_KEY = "skymatch_offsets_file"

_SUFFIX_RE = re.compile(
    r"(?i)(_segm|_cal|_skymatch|_cat|_i2d|_crf)?(\.asdf|\.parquet|\.h5|\.fits)?$"
)


def normalize_exposure_key(name: str | Path) -> str:
    """
    Normalize an exposure filename/stem for offset table lookup.

    Strips directories and common Roman product suffixes so ``..._cal.asdf``,
    ``..._segm.asdf``, and bare stems map to the same key.
    """
    stem = Path(str(name)).name
    prev = None
    while prev != stem:
        prev = stem
        stem = _SUFFIX_RE.sub("", stem)
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
) -> Path:
    """
    Write one exposure's HEALPix skyval samples to a parquet file.

    Columns: ``healpix17``, ``data``, ``err``. Coarse coverage is stored in a
    sidecar ``*_cov.parquet`` when provided.
    """
    outfile = Path(outfile)
    if outfile.suffix.lower() != ".parquet":
        outfile = outfile.with_suffix(".parquet")
    outfile.parent.mkdir(parents=True, exist_ok=True)

    table = Table(
        {
            "healpix17": np.asarray(pixel, dtype=np.int64),
            "data": np.asarray(value, dtype=np.float32),
            "err": np.asarray(unc, dtype=np.float32),
        }
    )
    if filename is not None:
        table["filename"] = np.full(len(table), normalize_exposure_key(filename))
    table.write(outfile, format="parquet", overwrite=True)

    if coverage is not None:
        cov_path = outfile.with_name(outfile.stem + "_cov.parquet")
        cov_table = Table(
            {"healpix11": np.asarray(coverage, dtype=np.int64).reshape(-1)}
        )
        cov_table.write(cov_path, format="parquet", overwrite=True)
    return outfile


def load_skyval_parquet(path: str | Path) -> dict[str, np.ndarray]:
    """
    Load skyval samples from a parquet file.

    Returns
    -------
    sample : dict
        Keys ``pixel``, ``value``, ``unc``, and optional ``coverage``.
    """
    path = Path(path)
    table = Table.read(path, format="parquet")
    if "healpix17" in table.colnames:
        pixel = np.asarray(table["healpix17"], dtype=np.int64)
    elif "pixel" in table.colnames:
        pixel = np.asarray(table["pixel"], dtype=np.int64)
    else:
        raise ValueError(f"{path} missing healpix17/pixel column")

    if "data" in table.colnames:
        value = np.asarray(table["data"], dtype=np.float64)
    elif "value" in table.colnames:
        value = np.asarray(table["value"], dtype=np.float64)
    else:
        raise ValueError(f"{path} missing data/value column")

    if "err" in table.colnames:
        unc = np.asarray(table["err"], dtype=np.float64)
    elif "unc" in table.colnames:
        unc = np.asarray(table["unc"], dtype=np.float64)
    else:
        raise ValueError(f"{path} missing err/unc column")

    sample: dict[str, np.ndarray] = {"pixel": pixel, "value": value, "unc": unc}

    cov_path = path.with_name(path.stem + "_cov.parquet")
    if cov_path.exists():
        cov = Table.read(cov_path, format="parquet")
        col = "healpix11" if "healpix11" in cov.colnames else cov.colnames[0]
        sample["coverage"] = np.asarray(cov[col], dtype=np.int64)
    return sample


def write_offset_table(
    outfile: str | Path,
    filenames: Sequence[str],
    offsets: Sequence[float],
    *,
    as_pipeline_level: bool = False,
) -> Path:
    """
    Write a skymatch offsets parquet table.

    Parameters
    ----------
    outfile
        Destination ``.parquet`` path.
    filenames
        Exposure keys (will be normalized).
    offsets
        Solver epsilon (additive) or pipeline levels if ``as_pipeline_level``.
    as_pipeline_level
        If True, store values as ``level`` (to subtract). If False (default),
        store solver ``epsilon`` (to add to skyvals) and also ``level=-epsilon``.
    """
    outfile = Path(outfile)
    if outfile.suffix.lower() != ".parquet":
        outfile = outfile.with_suffix(".parquet")
    outfile.parent.mkdir(parents=True, exist_ok=True)

    keys = [normalize_exposure_key(f) for f in filenames]
    offs = np.asarray(offsets, dtype=np.float64)
    if offs.shape[0] != len(keys):
        raise ValueError("filenames and offsets must have the same length")

    if as_pipeline_level:
        level = offs
        epsilon = -offs
    else:
        epsilon = offs
        level = -offs

    table = Table(
        {
            "filename": keys,
            "epsilon": np.asarray(epsilon, dtype=np.float64),
            "level": np.asarray(level, dtype=np.float64),
        }
    )
    table.write(outfile, format="parquet", overwrite=True)
    log.info("Wrote skymatch offsets for %d exposures to %s", len(keys), outfile)
    return outfile


def load_offset_table(path: str | Path) -> dict[str, float]:
    """
    Load offsets parquet into ``{normalized_filename: pipeline level}``.

    Accepts columns:

    * ``level`` — used directly as ``meta.background.level``
    * else ``epsilon`` / ``offset`` / ``offsets`` — converted with ``level = -epsilon``
    * ``filename`` / ``file`` / ``files`` / ``name`` for keys
    """
    path = Path(path)
    table = Table.read(path, format="parquet")

    name_col = next(
        (c for c in ("filename", "file", "files", "name") if c in table.colnames),
        None,
    )
    if name_col is None:
        raise ValueError(f"{path} missing filename column")

    if "level" in table.colnames:
        levels = np.asarray(table["level"], dtype=np.float64)
    else:
        eps_col = next(
            (c for c in ("epsilon", "offset", "offsets") if c in table.colnames),
            None,
        )
        if eps_col is None:
            raise ValueError(f"{path} missing level/epsilon/offset column")
        levels = -np.asarray(table[eps_col], dtype=np.float64)

    out: dict[str, float] = {}
    for name, level in zip(table[name_col], levels, strict=True):
        key = normalize_exposure_key(str(name))
        out[key] = float(level)
    return out


def lookup_offset(filename: str | Path, table: Mapping[str, float]) -> float:
    """Look up pipeline level for an exposure name; raise if missing."""
    key = normalize_exposure_key(filename)
    if key in table:
        return table[key]
    raw = Path(str(filename)).stem
    if raw in table:
        return table[raw]
    raise KeyError(
        f"No skymatch offset for {filename!r} (key={key!r}); have {len(table)} entries"
    )


def resolve_offsets_path_from_asn(
    asn: Mapping,
    *,
    asn_path: str | Path | None = None,
) -> Path | None:
    """
    Resolve ``skymatch_offsets_file`` from an ASN dict if present.

    Checks product-level then top-level. Relative paths are resolved against
    the ASN file's directory when ``asn_path`` is given.
    """
    candidates: list[str] = []
    products = asn.get("products") or []
    if products:
        prod = products[0]
        if isinstance(prod, Mapping) and prod.get(ASN_OFFSETS_KEY):
            candidates.append(str(prod[ASN_OFFSETS_KEY]))
    if asn.get(ASN_OFFSETS_KEY):
        candidates.append(str(asn[ASN_OFFSETS_KEY]))
    if not candidates:
        return None

    path = Path(candidates[0])
    if path.is_absolute():
        return path
    if asn_path is not None:
        return (Path(asn_path).parent / path).resolve()
    return path.resolve()

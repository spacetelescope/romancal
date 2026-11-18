#!/usr/bin/env python
"""
Match HEALPix skyval parquet files and write a skymatch levels table.

Example
-------
$ roman_match_skyvals_parquet -o offsets.parquet skyvals/*.parquet
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from romancal.skymatch.healpix_match import match_healpix_skyvals
from romancal.skymatch.offset_io import (
    load_skyval_parquet,
    normalize_exposure_key,
    write_offset_table,
)


def _expand(inputs: list[str]) -> list[Path]:
    out: list[Path] = []
    for item in inputs:
        p = Path(item)
        if p.is_dir():
            out.extend(sorted(p.glob("*.parquet")))
        elif any(c in item for c in "*?["):
            parent = p.parent if str(p.parent) != "" else Path(".")
            out.extend(sorted(parent.glob(p.name)))
        else:
            out.append(p)
    # drop coverage sidecars
    out = [p for p in out if not p.stem.endswith("_cov")]
    seen = set()
    uniq = []
    for p in out:
        r = p.resolve()
        if r not in seen:
            seen.add(r)
            uniq.append(p)
    return uniq


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Solve skymatch levels from skyval parquet files."
    )
    parser.add_argument(
        "inputs", nargs="+", help="Skyval parquet files, dirs, or globs"
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output offsets parquet path (filename and level columns)",
    )
    parser.add_argument("-v", "--verbose", action="count", default=0)
    args = parser.parse_args(argv)

    level = [logging.WARNING, logging.INFO, logging.DEBUG][min(args.verbose, 2)]
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")

    paths = _expand(args.inputs)
    if not paths:
        raise SystemExit("No skyval parquet inputs found")

    samples = []
    names = []
    for path in paths:
        sample = load_skyval_parquet(path)
        samples.append(sample)
        names.append(normalize_exposure_key(path.stem))
        logging.info("Loaded %s", path)

    levels, names_out = match_healpix_skyvals(samples, names=names)
    out = write_offset_table(args.output, names_out, levels)
    print(f"Wrote levels for {len(names_out)} exposures to {out}")


if __name__ == "__main__":
    main()

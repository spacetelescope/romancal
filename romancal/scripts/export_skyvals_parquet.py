#!/usr/bin/env python
"""
Export L2 segmentation skyvals to parquet and optionally run healpix sky match.

Example
-------
$ roman_export_skyvals_parquet -o ./skyvals *_segm.asdf
$ roman_export_skyvals_parquet -o ./skyvals --match-offsets offsets.parquet *_segm.asdf
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
from roman_datamodels import datamodels

from romancal.skymatch.healpix_match import match_healpix_skyvals
from romancal.skymatch.offset_io import (
    normalize_exposure_key,
    write_offset_table,
    write_skyval_parquet,
)


def _expand(inputs: list[str]) -> list[Path]:
    out: list[Path] = []
    for item in inputs:
        p = Path(item)
        if p.is_dir():
            out.extend(sorted(p.glob("*_segm.asdf")))
        elif any(c in item for c in "*?["):
            parent = p.parent if str(p.parent) != "" else Path(".")
            out.extend(sorted(parent.glob(p.name)))
        else:
            out.append(p)
    # unique
    seen = set()
    uniq = []
    for p in out:
        r = p.resolve()
        if r not in seen:
            seen.add(r)
            uniq.append(p)
    return uniq


def _stem_from_segm(path: Path) -> str:
    stem = path.stem
    if stem.endswith("_segm"):
        stem = stem[: -len("_segm")]
    return stem


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Export segmentation skyvals to parquet "
            "(healpix17/data/err) and optionally solve skymatch offsets."
        )
    )
    parser.add_argument("inputs", nargs="+", help="Segm ASDF files, dirs, or globs")
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Directory for per-exposure skyval parquet files",
    )
    parser.add_argument(
        "--match-offsets",
        default="",
        help="If set, run healpix match and write this offsets parquet path",
    )
    parser.add_argument("-v", "--verbose", action="count", default=0)
    args = parser.parse_args(argv)

    level = [logging.WARNING, logging.INFO, logging.DEBUG][min(args.verbose, 2)]
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")

    paths = _expand(args.inputs)
    if not paths:
        raise SystemExit("No input segmentation files found")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    samples = []
    names = []
    written = []
    for path in paths:
        with datamodels.open(path) as model:
            if "skyvals" not in model:
                raise SystemExit(
                    f"{path} has no skyvals; run SourceCatalogStep with compute_skyvals=True"
                )
            sky = np.asarray(model.skyvals)
            cov = np.asarray(model.healpix11_cov) if "healpix11_cov" in model else None
        stem = _stem_from_segm(path)
        out = write_skyval_parquet(
            outdir / f"{stem}.parquet",
            pixel=sky["healpix17"],
            value=sky["data"],
            unc=sky["err"],
            coverage=cov,
            filename=stem,
        )
        written.append(out)
        samples.append(
            {
                "pixel": np.asarray(sky["healpix17"], dtype=np.int64),
                "value": np.asarray(sky["data"], dtype=np.float64),
                "unc": np.asarray(sky["err"], dtype=np.float64),
                "coverage": cov,
            }
        )
        names.append(normalize_exposure_key(stem))
        logging.info("Wrote %s", out)

    print(f"Wrote {len(written)} skyval parquet file(s) to {outdir}")
    for p in written:
        print(p)

    if args.match_offsets:
        epsilon, names_out = match_healpix_skyvals(samples, names=names)
        off_path = write_offset_table(args.match_offsets, names_out, epsilon)
        print(f"Wrote offsets to {off_path}")


if __name__ == "__main__":
    main()

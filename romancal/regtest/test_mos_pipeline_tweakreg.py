"""High Level Pipeline (HLP) regression tests on tweakreg'd Level 2 files.

Purpose
-------
romancal#1236: the existing HLP / mosaic regtests (`test_mos_pipeline.py`)
run on plain ``*_cal.asdf`` members produced by ELP. Those products can differ
materially from the same exposures after absolute/relative WCS alignment via
tweakreg. This module exercises:

1. ``TweakRegStep`` on the multi-member L3 association cal set
2. ``MosaicPipeline`` (``roman_mos``) on an ASN whose members are the
   tweakreg'd outputs

Assertions cover pipeline contract (steps complete, mosaic product type,
side products, tweakreg markers) and a required ``compare_asdf`` against a
dedicated tweakreg-fed coadd truth on Artifactory.

Artifactory inputs (read-only download via ``rtdata``)
-----------------------------------------------------
* ``WFI/image/L3_regtest_asn.json`` and its ``*_cal.asdf`` members
* Optional matching ``*_cat.parquet`` catalogs when already published
* Truth: ``truth/WFI/image/r0000101001001001001_f158_tweakreg_coadd.asdf``
  (MosaicModel coadd from HLP run on tweakreg'd L2 — not the cal-only HLP truth)

Catalogs that are not on Artifactory are produced in the test jail with
``SourceCatalogStep`` from the downloaded cal files (no Artifactory upload
from this module).

Environment
-----------
Requires bigdata access (``TEST_BIGDATA`` / ci_watson) like other regtests.
This module does not upload/okify truths; missing truth fails with an explicit
message so the suite is not mistaken for a science regression.
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path

import pytest
import roman_datamodels as rdm
from ci_watson.artifactory_helpers import BigdataError, get_bigdata
from roman_datamodels.datamodels import ImageModel, MosaicModel

from romancal.pipeline.mosaic_pipeline import MosaicPipeline
from romancal.source_catalog import SourceCatalogStep
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf

log = logging.getLogger(__name__)

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]

# Same L3 ASN family as test_mos_pipeline.py (cal members).
L3_ASN_REMOTE = "WFI/image/L3_regtest_asn.json"
# Distinct product stem so the tweakreg-fed coadd does not collide with the
# cal-only HLP truth (r0000101001001001001_f158_coadd.asdf).
L3_PRODUCT_NAME = "r0000101001001001001_f158_tweakreg"
TWEAKREG_SUFFIX = "tweakregstep"
# Required dedicated truth: MosaicModel coadd from HLP on tweakreg'd L2.
TWEAKREG_COADD_OUTPUT = f"{L3_PRODUCT_NAME}_coadd.asdf"
TWEAKREG_COADD_TRUTH = f"truth/WFI/image/{TWEAKREG_COADD_OUTPUT}"


def _cal_to_catalog_name(cal_name: str) -> str:
    """Map ``*_cal.asdf`` member name to sibling source catalog filename."""
    if cal_name.endswith("_cal.asdf"):
        return cal_name.replace("_cal.asdf", "_cat.parquet")
    stem = Path(cal_name).stem
    # Fallback: strip last underscore suffix and append _cat.parquet
    if "_" in stem:
        return stem.rsplit("_", 1)[0] + "_cat.parquet"
    return stem + "_cat.parquet"


def _cal_to_tweakreg_name(cal_name: str) -> str:
    """Map cal member basename to TweakRegStep output basename."""
    path = Path(cal_name)
    stem = path.stem
    if stem.endswith("_cal"):
        stem = stem[: -len("_cal")]
    return f"{stem}_{TWEAKREG_SUFFIX}.asdf"


def _ensure_catalogs_for_cal_members(rtdata, cal_members: list[str]) -> None:
    """Ensure each cal has a sibling ``*_cat.parquet`` for TweakRegStep.

    Prefer an existing Artifactory catalog when present. Otherwise run
    ``SourceCatalogStep`` on the local cal (L3 ASN members often ship without
    published catalogs — only the cals are required for HLP).
    """
    asn_dir = os.path.dirname(L3_ASN_REMOTE)
    for cal_name in cal_members:
        cat_name = _cal_to_catalog_name(cal_name)
        if os.path.isfile(cat_name):
            continue

        remote_cat = os.path.join(asn_dir, cat_name)
        try:
            # Do not use get_data: it overwrites rtdata.input/asn bookkeeping.
            get_bigdata(
                rtdata._inputs_root,
                rtdata._env,
                remote_cat,
                docopy=rtdata.docopy,
            )
        except BigdataError:
            log.info(
                "Catalog %s not on bigdata; generating with SourceCatalogStep from %s",
                cat_name,
                cal_name,
            )
            SourceCatalogStep.call(
                cal_name,
                save_results=True,
                suffix="cat",
            )

        if not os.path.isfile(cat_name):
            raise FileNotFoundError(
                f"Failed to obtain or generate source catalog for tweakreg: {cat_name}"
            )


def _attach_catalogs_to_cal_files(member_names: list[str]) -> None:
    """Point each on-disk cal model at its ``*_cat.parquet`` for TweakRegStep.

    TweakRegStep requires ``meta.source_catalog.tweakreg_catalog`` or
    ``meta.source_catalog.tweakreg_catalog_name``. ELP-produced cals used as
    HLP inputs may not carry a usable path after Artifactory download, so we
    set the catalog name explicitly before running tweakreg.
    """
    for cal_name in member_names:
        cat_name = _cal_to_catalog_name(cal_name)
        if not os.path.isfile(cat_name):
            raise FileNotFoundError(
                f"Source catalog for tweakreg not found next to cal: {cat_name}"
            )
        with rdm.open(cal_name) as model:
            if "source_catalog" not in model.meta:
                model.meta["source_catalog"] = {}
            # Prefer file path so catalogs are not embedded into the asdf.
            model.meta.source_catalog["tweakreg_catalog_name"] = os.path.abspath(
                cat_name
            )
            # Drop any stale in-memory catalog reference from prior processing.
            model.meta.source_catalog.pop("tweakreg_catalog", None)
            model.to_asdf(cal_name)


def _write_tweakreg_asn(cal_member_names: list[str], asn_path: str) -> str:
    """Write an L3 ASN whose science members are tweakreg'd L2 products."""
    members = []
    for cal_name in cal_member_names:
        tweak_name = _cal_to_tweakreg_name(cal_name)
        if not os.path.isfile(tweak_name):
            raise FileNotFoundError(f"Expected tweakreg output missing: {tweak_name}")
        members.append({"expname": tweak_name, "exptype": "science"})

    asn = {
        "asn_type": "None",
        "asn_rule": "DMS_ELPP_Base",
        "version_id": None,
        "code_version": "regtest",
        "target": "None",
        "degraded_status": "No known degraded exposures in association.",
        "program": "noprogram",
        "constraints": "No constraints",
        "asn_id": "a3001_tweakreg",
        "asn_pool": "none",
        "skycell_wcs_info": "none",
        "products": [
            {
                "name": L3_PRODUCT_NAME,
                "members": members,
            }
        ],
    }
    with open(asn_path, "w") as fh:
        json.dump(asn, fh, indent=2)
    return asn_path


@pytest.fixture(scope="module")
def run_mos_tweakreg(rtdata_module, resource_tracker):
    """Run tweakreg on L3 cal members, then HLP/mosaic on tweakreg'd files.

    Purpose: build the tweakreg → HLP path described in romancal#1236 without
    modifying the cal-only ``test_mos_pipeline`` coverage.
    """
    rtdata = rtdata_module

    # Download ASN + cal members (read-only from Artifactory / bigdata).
    rtdata.get_asn(L3_ASN_REMOTE)
    asn = rtdata.asn
    cal_members = [
        m["expname"]
        for product in asn["products"]
        for m in product["members"]
        if m.get("exptype", "science") == "science"
    ]
    assert cal_members, "L3 ASN contains no science members"

    # Catalogs are optional on Artifactory for L3 cal members; generate any
    # missing ones locally so CI does not depend on unpublished parquet files.
    _ensure_catalogs_for_cal_members(rtdata, cal_members)
    _attach_catalogs_to_cal_files(cal_members)

    # Run tweakreg on the multi-member cal set (single association-like call
    # via a temporary ASN of cal files so relative/absolute matching sees
    # the full group — same members as the L3 ASN).
    cal_asn_path = "L3_cal_for_tweakreg_asn.json"
    with open(cal_asn_path, "w") as fh:
        json.dump(
            {
                "asn_type": "None",
                "asn_rule": "DMS_ELPP_Base",
                "version_id": None,
                "code_version": "regtest",
                "target": "None",
                "degraded_status": "No known degraded exposures in association.",
                "program": "noprogram",
                "constraints": "No constraints",
                "asn_id": "a3001_cal_tweakreg",
                "asn_pool": "none",
                "skycell_wcs_info": "none",
                "products": [
                    {
                        "name": f"{L3_PRODUCT_NAME}_cal_tweakreg",
                        "members": [
                            {"expname": name, "exptype": "science"}
                            for name in cal_members
                        ],
                    }
                ],
            },
            fh,
            indent=2,
        )

    tweakreg_args = [
        "romancal.step.TweakRegStep",
        cal_asn_path,
        f"--suffix={TWEAKREG_SUFFIX}",
        "--save_results=True",
    ]
    with resource_tracker.track():
        RomanStep.from_cmdline(tweakreg_args)

    # Sanity: each member produced a tweakreg'd ImageModel with COMPLETE step.
    for cal_name in cal_members:
        tweak_name = _cal_to_tweakreg_name(cal_name)
        assert os.path.isfile(tweak_name), f"missing tweakreg product {tweak_name}"
        with rdm.open(tweak_name) as model:
            assert isinstance(model, ImageModel)
            assert model.meta.cal_step.tweakreg == "COMPLETE"
            assert "v2v3corr" in model.meta.wcs.available_frames

    # HLP on tweakreg'd members.
    tweak_asn_path = "L3_tweakreg_asn.json"
    _write_tweakreg_asn(cal_members, tweak_asn_path)

    rtdata.input = tweak_asn_path
    rtdata.output = TWEAKREG_COADD_OUTPUT

    mos_args = [
        "roman_mos",
        tweak_asn_path,
        "--save_results=True",
    ]
    with resource_tracker.track():
        MosaicPipeline.from_cmdline(mos_args)

    # Required truth download (same pattern as test_mos_pipeline.run_mos).
    # Fail with an explicit message if the dedicated tweakreg-HLP truth has
    # not been published yet — distinct from a science compare_asdf mismatch.
    try:
        rtdata.get_truth(TWEAKREG_COADD_TRUTH)
    except BigdataError as exc:
        pytest.fail(
            "Tweakreg-HLP truth file is not available yet on bigdata/Artifactory. "
            "This failure means the baseline coadd has not been published; it is "
            "not a MosaicPipeline/tweakreg science regression.\n"
            f"Expected remote path: {TWEAKREG_COADD_TRUTH}\n"
            f"Local pipeline output to promote after review: {TWEAKREG_COADD_OUTPUT}\n"
            "Generate by running this module's fixture path, then upload/okify the "
            "coadd ASDF to that truth path (do not overwrite the cal-only HLP truth "
            f"r0000101001001001001_f158_coadd.asdf).\n"
            f"Original error: {exc}"
        )

    # Stash member list for downstream tests (no Artifactory upload).
    rtdata.tweakreg_members = [_cal_to_tweakreg_name(n) for n in cal_members]
    return rtdata


@pytest.fixture(scope="module")
def output_filename(run_mos_tweakreg):
    return run_mos_tweakreg.output


@pytest.fixture(scope="module")
def output_model(output_filename):
    with rdm.open(output_filename) as model:
        yield model


@pytest.fixture(scope="module")
def truth_filename(run_mos_tweakreg):
    return run_mos_tweakreg.truth


def test_log_tracked_resources(log_tracked_resources, run_mos_tweakreg):
    """Purpose: emit resource-tracker summary for the tweakreg+HLP fixture."""
    log_tracked_resources()


def test_tweakreg_members_used(run_mos_tweakreg):
    """Purpose: HLP inputs are tweakreg'd L2 files, not plain cal products."""
    for name in run_mos_tweakreg.tweakreg_members:
        assert name.endswith(f"_{TWEAKREG_SUFFIX}.asdf")
        assert os.path.isfile(name)
        with rdm.open(name) as model:
            assert model.meta.cal_step.tweakreg == "COMPLETE"


def test_output_matches_truth(output_filename, truth_filename, ignore_asdf_paths):
    """Purpose: tweakreg-fed HLP coadd matches the dedicated Artifactory truth."""
    diff = compare_asdf(output_filename, truth_filename, **ignore_asdf_paths)
    assert diff.identical, diff.report()


def test_output_is_mosaic(output_model):
    """Purpose: HLP on tweakreg'd inputs still yields a MosaicModel coadd."""
    assert isinstance(output_model, MosaicModel)


@pytest.mark.parametrize(
    "step_name",
    (
        "skymatch",
        "outlier_detection",
        "resample",
    ),
)
def test_hlp_steps_ran(output_model, step_name):
    """Purpose: core HLP calibration steps complete on the tweakreg-fed path."""
    assert getattr(output_model.meta.cal_step, step_name) == "COMPLETE"


@pytest.mark.parametrize("suffix", ("cat.parquet", "segm.asdf"))
def test_side_products_exist(output_filename, suffix):
    """Purpose: mosaic source catalog and segmentation products are written."""
    expected_filename = output_filename.rsplit("_", 1)[0] + f"_{suffix}"
    assert os.path.isfile(expected_filename)

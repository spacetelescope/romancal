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

Assertions focus on pipeline contract (steps complete, mosaic product type,
side products exist, tweakreg markers present). A full ``compare_asdf`` truth
lock is optional and only runs when a dedicated truth file is already present
on Artifactory — this test never uploads results to Artifactory.

Artifactory inputs (read-only download via ``rtdata``)
-----------------------------------------------------
* ``WFI/image/L3_regtest_asn.json`` and its ``*_cal.asdf`` members
* Matching ``*_cat.parquet`` catalogs for each cal (same basename stem)

Environment
-----------
Requires bigdata access (``TEST_BIGDATA`` / ci_watson) like other regtests.
No Artifactory write/okify upload is performed by this module.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest
import roman_datamodels as rdm
from ci_watson.artifactory_helpers import get_bigdata
from roman_datamodels.datamodels import ImageModel, MosaicModel

from romancal.pipeline.mosaic_pipeline import MosaicPipeline
from romancal.stpipe import RomanStep

from .regtestdata import compare_asdf

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]

# Same L3 ASN family as test_mos_pipeline.py (cal members).
L3_ASN_REMOTE = "WFI/image/L3_regtest_asn.json"
# Product name from L3_regtest_asn.json — mosaic output basename stem.
L3_PRODUCT_NAME = "r0000101001001001001_f158"
TWEAKREG_SUFFIX = "tweakregstep"
# Optional dedicated truth for the tweakreg-fed coadd (download-only if present).
TWEAKREG_COADD_TRUTH = f"truth/WFI/image/{L3_PRODUCT_NAME}_tweakreg_coadd.asdf"


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


def _attach_catalogs_to_cal_files(member_names: list[str]) -> None:
    """Point each on-disk cal model at its downloaded ``*_cat.parquet``.

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

    # Fetch sibling source catalogs required by TweakRegStep.
    # get_asn leaves input_remote as the ASN relative path (e.g. WFI/image/...json).
    asn_dir = os.path.dirname(L3_ASN_REMOTE)
    for cal_name in cal_members:
        cat_name = _cal_to_catalog_name(cal_name)
        # Use get_bigdata path only; do not go through get_data so we do not
        # overwrite rtdata.input/asn bookkeeping mid-fixture.
        get_bigdata(
            rtdata._inputs_root,
            rtdata._env,
            os.path.join(asn_dir, cat_name),
            docopy=rtdata.docopy,
        )

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

    output = f"{L3_PRODUCT_NAME}_coadd.asdf"
    rtdata.input = tweak_asn_path
    rtdata.output = output
    # Record intended truth location for local failure bookkeeping only.
    # Do not download or upload Artifactory truth from this fixture.
    rtdata.truth_remote = TWEAKREG_COADD_TRUTH

    mos_args = [
        "roman_mos",
        tweak_asn_path,
        "--save_results=True",
    ]
    with resource_tracker.track():
        MosaicPipeline.from_cmdline(mos_args)

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


def test_optional_truth_compare(run_mos_tweakreg, ignore_asdf_paths):
    """Purpose: if a dedicated tweakreg-HLP truth exists, compare coadd to it.

    Skips when the optional truth is not on bigdata. Never uploads or okifies
    truth to Artifactory from this test.
    """
    rtdata = run_mos_tweakreg
    try:
        rtdata.get_truth(TWEAKREG_COADD_TRUTH)
    except Exception as exc:
        # Missing optional truth (BigdataError/IO) should not fail the suite.
        pytest.skip(f"Optional tweakreg HLP truth not available: {exc}")

    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()

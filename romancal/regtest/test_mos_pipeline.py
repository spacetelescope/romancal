""" Roman tests for the High Level Pipeline """

import json
import os
from pathlib import Path

import asdf
import pytest
import roman_datamodels as rdm
from astropy.units import Quantity

from romancal.associations.asn_from_list import asn_from_list
from romancal.pipeline.mosaic_pipeline import MosaicPipeline

from ..associations.association_io import json as asn_json
from .regtestdata import compare_asdf

# mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.soctests]



class RegtestFileModifier:
    # TODO: remove this entire class once the units
    #  have been removed from the regtest files

    def __init__(self, rtdata):
        self.rtdata = rtdata
        self.updated_asn_fname = None
        self.truth_parent = Path(rtdata.truth).parent
        self.input_parent = Path(rtdata.input).parent
        self.truth_relative_path = Path(self.truth_parent).relative_to(
            self.input_parent
        )
        self.truth_path = self.truth_relative_path / f"{Path(self.rtdata.truth).name}"

    @staticmethod
    def create_unitless_file(input_filename: str, output_filename: str) -> None:
        with asdf.config_context() as cfg:
            cfg.validate_on_read = False
            cfg.validate_on_save = False
            af = asdf.open(input_filename)

            for attr in af.tree["roman"]:
                item = getattr(af.tree["roman"], attr)
                if isinstance(item, Quantity):
                    setattr(af.tree["roman"], attr, item.value)

            for attr in af.tree["roman"].meta.photometry:
                item = getattr(af.tree["roman"].meta.photometry, attr)
                if isinstance(item, Quantity):
                    setattr(af.tree["roman"].meta.photometry, attr, item.value)

            af.write_to(output_filename)

    def create_new_asn_file(self, output_filename_list: list):
        updated_asn = asn_from_list(
            output_filename_list,
            product_name=f"{self.rtdata.asn['products'][0]['name']}_no_units",
        )
        updated_asn["target"] = "none"

        current_asn_fname = Path(self.rtdata.input)
        self.updated_asn_fname = (
            f"{current_asn_fname.stem}_no_units{current_asn_fname.suffix}"
        )

        _, serialized_updated_asn = asn_json.dump(updated_asn)
        with open(self.updated_asn_fname, "w") as f:
            json.dump(
                json.loads(serialized_updated_asn), f, indent=4, separators=(",", ": ")
            )

    def update_rtdata(self):
        rtdata_root_path = Path(self.rtdata.input).parent
        self.rtdata.input = f"{rtdata_root_path}/{Path(self.updated_asn_fname)}"
        # r0099101001001001001_F158_visit_no_units_i2d.asdf
        self.rtdata.output = f"{rtdata_root_path}/{Path(self.rtdata.output.split('_i2d')[0]).stem}_no_units_i2d{Path(self.rtdata.output).suffix}"

    def prepare_regtest_input_files(self):
        input_filenames = [
            x["expname"] for x in self.rtdata.asn["products"][0]["members"]
        ]
        input_filenames.append(str(self.truth_path))
        output_filename_list = []
        # include truth file
        for input_filename in input_filenames:
            fname = Path(input_filename)
            if str(fname).startswith(str(self.truth_relative_path)):
                output_filename = Path(
                    f"{str(fname).split('_i2d.asdf')[0]}_no_units_i2d{fname.suffix}"
                )
                self.rtdata.truth = str(self.truth_parent / output_filename.name)
            else:
                output_filename = f"{fname.stem}_no_units{fname.suffix}"
                output_filename_list.append(output_filename)
            self.create_unitless_file(input_filename, output_filename)

        self.create_new_asn_file(output_filename_list)
        self.update_rtdata()



@pytest.fixture(scope="module")
def run_mos(rtdata_module):
    rtdata = rtdata_module

    rtdata.get_asn("WFI/image/L3_regtest_asn.json")

    # Test Pipeline
    output = "r0099101001001001001_F158_visit_i2d.asdf"
    rtdata.output = output

    rtdata.get_truth(f"truth/WFI/image/{output}")

    fixer = RegtestFileModifier(rtdata)
    fixer.prepare_regtest_input_files()

    args = [
        "roman_mos",
        rtdata.input,
    ]
    MosaicPipeline.from_cmdline(args)

    # rtdata.get_truth(f"truth/WFI/image/{output}")
    return rtdata


@pytest.fixture(scope="module")
def output_filename(run_mos):
    return run_mos.output


@pytest.fixture(scope="module")
def output_model(output_filename):
    with rdm.open(output_filename) as model:
        yield model


@pytest.fixture(scope="module")
def truth_filename(run_mos):
    return run_mos.truth


@pytest.fixture(scope="module")
def thumbnail_filename(output_filename):
    thumbnail_filename = output_filename.rsplit("_", 1)[0] + "_thumb.png"
    preview_cmd = f"stpreview to {output_filename} {thumbnail_filename} 256 256 roman"
    os.system(preview_cmd)  # nosec
    return thumbnail_filename


@pytest.fixture(scope="module")
def preview_filename(output_filename):
    preview_filename = output_filename.rsplit("_", 1)[0] + "_preview.png"
    preview_cmd = f"stpreview to {output_filename} {preview_filename} 1080 1080 roman"
    os.system(preview_cmd)  # nosec
    return preview_filename


def test_output_matches_truth(output_filename, truth_filename, ignore_asdf_paths):
    # DMS356
    diff = compare_asdf(output_filename, truth_filename, **ignore_asdf_paths)
    assert diff.identical, diff.report()


def test_thumbnail_exists(thumbnail_filename):
    # DMS356
    assert os.path.isfile(thumbnail_filename)


def test_preview_exists(preview_filename):
    # DMS356
    assert os.path.isfile(preview_filename)


@pytest.mark.parametrize("suffix", ("cat", "segm"))
def test_file_exists(output_filename, suffix):
    # DMS374 for catalog and segm
    expected_filename = output_filename.rsplit("_", 1)[0] + f"_{suffix}.asdf"
    assert os.path.isfile(expected_filename)


def test_output_is_mosaic(output_model):
    # DMS356
    assert isinstance(output_model, rdm.datamodels.MosaicModel)


@pytest.mark.parametrize(
    "step_name",
    (
        "skymatch",
        "outlier_detection",
        "resample",
    ),
)
def test_steps_ran(output_model, step_name):
    # DMS356
    # DMS400 for skymatch
    # DMS86 for outlier_detection and resample
    assert getattr(output_model.meta.cal_step, step_name) == "COMPLETE"


def test_added_background(output_model):
    # DMS400
    assert hasattr(output_model.meta.individual_image_meta, "background")


def test_added_background_level(output_model):
    # DMS400
    assert any(output_model.meta.individual_image_meta.background["level"] != 0)

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


def passfail(bool_expr):
    """set pass fail"""
    if bool_expr:
        return "Pass"
    return "Fail"


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


@pytest.mark.bigdata
@pytest.mark.soctests
def test_level3_mos_pipeline(rtdata, ignore_asdf_paths):
    """Tests for level 3 processing requirements DMS356"""
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

    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()

    # Generate thumbnail image
    input_file = "r0099101001001001001_F158_visit_i2d.asdf"
    thumbnail_file = "r0099101001001001001_F158_visit_thumb.png"

    preview_cmd = f"stpreview to {input_file} {thumbnail_file} 256 256 roman"
    os.system(preview_cmd)  # nosec

    # Generate preview image
    input_file = "r0099101001001001001_F158_visit_i2d.asdf"
    preview_file = "r0099101001001001001_F158_visit_preview.png"
    preview_cmd = f"stpreview to {input_file} {preview_file} 1080 1080 roman"
    os.system(preview_cmd)  # nosec

    # expected catalog and segmentation files
    catalog_file = "r0099101001001001001_F158_visit_cat.asdf"
    segm_file = "r0099101001001001001_F158_visit_segm.asdf"

    # Perform DMS tests
    # Initial prep
    model = rdm.open(rtdata.output)
    pipeline = MosaicPipeline()

    # DMS356 result is an ImageModel
    pipeline.log.info(
        "DMS356 MSG: Testing that result is a Level 3 mosaic model......."
        + passfail(isinstance(model, rdm.datamodels.MosaicModel))
    )

    # DMS356 Test that skymatch step is complete
    pipeline.log.info(
        "Status of the step:             skymatch    "
        + str(model.meta.cal_step.skymatch)
    )
    # DMS356 Test that the thumbnail image exists
    pipeline.log.info(
        "Status of the step:             thumbnail image    "
        + passfail(os.path.isfile(thumbnail_file))
    )
    # DMS356 Test that the preview image exists
    pipeline.log.info(
        "Status of the step:             preview image    "
        + passfail(os.path.isfile(preview_file))
    )
    # DMS374 Test that the output catalog exists
    pipeline.log.info(
        "Check that the catalog file exists   " + passfail(os.path.isfile(catalog_file))
    )
    # DMS374 Test that the segmentation file exists
    pipeline.log.info(
        "Check that the degmentation file exists   "
        + passfail(os.path.isfile(segm_file))
    )
    pipeline.log.info(
        "DMS400 MSG: Testing completion of skymatch in the Level 3  output......."
        + passfail(model.meta.cal_step.skymatch == "COMPLETE")
    )
    assert model.meta.cal_step.skymatch == "COMPLETE"
    pipeline.log.info(
        "Status of the step:             skymatch    "
        + str(model.meta.cal_step.skymatch)
    )
    pipeline.log.info(
        "DMS400 MSG: SkyMatchStep added meta.background? :"
        f'  {hasattr(model.meta.individual_image_meta, "background")}'
    )
    assert hasattr(model.meta.individual_image_meta, "background")

    pipeline.log.info(
        "DMS400 MSG: SkyMatchStep populated meta.background.level? :"
        f"  {any(model.meta.individual_image_meta.background['level'] != 0)}"
    )
    assert any(model.meta.individual_image_meta.background["level"] != 0)
    pipeline.log.info(
        "DMS86 MSG: Testing completion of outlier detection in the Level 3 image output......."
        + passfail(model.meta.cal_step.outlier_detection == "COMPLETE")
    )
    assert model.meta.cal_step.outlier_detection == "COMPLETE"
    pipeline.log.info(
        "Status of the step:             outlier_detection    "
        + str(model.meta.cal_step.outlier_detection)
    )
    assert model.meta.cal_step.resample == "COMPLETE"
    pipeline.log.info(
        "Status of the step:             resample          "
        + str(model.meta.cal_step.resample)
    )
    pipeline.log.info(
        "DMS86 MSG: Testing completion of resample in the Level 3 image output......."
        + passfail(model.meta.cal_step.resample == "COMPLETE")
    )


@pytest.mark.bigdata
@pytest.mark.soctests
def test_hlp_mosaic_pipeline(rtdata, ignore_asdf_paths):
    """Tests for level 3 mosaic requirements DMS373"""
    rtdata.get_asn("WFI/image/L3_mosaic_asn.json")

    # Test Pipeline
    output = "r0099101001001001001_r274dp63x31y81_prompt_F158_i2d.asdf"
    rtdata.output = output
    args = [
        "roman_mos",
        rtdata.input,
    ]
    MosaicPipeline.from_cmdline(args)
    rtdata.get_truth(f"truth/WFI/image/{output}")
    pipeline = MosaicPipeline()
    diff = compare_asdf(rtdata.output, rtdata.truth, **ignore_asdf_paths)
    assert diff.identical, diff.report()

    model = rdm.open(rtdata.output, lazy_load=False)

    pipeline.log.info(
        "DMS373 MSG: Testing the creation of a Level 3 mosaic image resampled to a skycell"
        + passfail(model.meta.cal_step.resample == "COMPLETE")
    )

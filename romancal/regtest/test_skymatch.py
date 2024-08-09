import pytest
from metrics_logger.decorators import metrics_logger

from romancal.datamodels.library import ModelLibrary
from romancal.pipeline.mosaic_pipeline import MosaicPipeline


@metrics_logger()
@pytest.mark.bigdata
def test_l3_output_contains_background_info(rtdata):
    """Test for the presence of meta.background in L3 files."""

    rtdata.get_asn("WFI/image/L3_regtest_asn.json")
    rtdata.output = f"{rtdata.asn['products'][0]['name']}_i2d.asdf"

    args = [
        "roman_mos",
        rtdata.input,
    ]
    MosaicPipeline.from_cmdline(args)

    step = MosaicPipeline()

    # instantiate a ModelLibrary object
    library = ModelLibrary([rtdata.output])
    with library:
        for index, model in enumerate(library):
            # Perform DMS tests
            step.log.info(
                "DMSXXX MSG: SkyMatchStep added meta.background? :"
                f'  {hasattr(model.meta.individual_image_meta, "background")}'
            )
            assert hasattr(model.meta.individual_image_meta, "background")

            step.log.info(
                "DMSXXX MSG: SkyMatchStep populated meta.background? :"
                f"  {all(v is not None for v in model.meta.individual_image_meta.background.values())}"
            )
            assert all(
                v is not None
                for v in model.meta.individual_image_meta.background.values()
            )

            library.shelve(model, index)

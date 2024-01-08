import json
from io import StringIO

import numpy as np
import pytest
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.models import RotationSequence3D, Scale, Shift
from gwcs import coordinate_frames as cf
from gwcs import wcs
from gwcs.geometry import CartesianToSpherical, SphericalToCartesian
from roman_datamodels import datamodels as rdm
from roman_datamodels import maker_utils

from romancal.datamodels import ModelContainer
from romancal.outlier_detection import OutlierDetectionStep, outlier_detection
from pathlib import Path
import os


def _create_tel2sky_model(input_dm):
    """
    Create the transform from telescope to sky.

    The transform is defined with a reference point in a Frame
    associated with the telescope (V2, V3) in arcsec, the corresponding
    reference point on sky (RA_REF, DEC_REF) in deg, and the position angle
    at the center of the aperture, ROLL_REF in deg.

    Parameters
    ----------
    input_dm : roman_datamodels.ImageModel
        A Roman image datamodel.

    Returns
    -------
    astropy.modeling.core.CompoundModel
        An astropy compound model with the transformations to map
        the telescope reference system onto the world reference systems.
    """
    v2_ref = input_dm.meta.wcsinfo.v2_ref / 3600
    v3_ref = input_dm.meta.wcsinfo.v3_ref / 3600
    roll_ref = input_dm.meta.wcsinfo.roll_ref
    ra_ref = input_dm.meta.wcsinfo.ra_ref
    dec_ref = input_dm.meta.wcsinfo.dec_ref

    angles = np.array([v2_ref, -v3_ref, roll_ref, dec_ref, -ra_ref])
    axes = "zyxyz"
    rot = RotationSequence3D(angles, axes_order=axes)

    # The sky rotation expects values in deg.
    # This should be removed when models work with quantities.
    model = (
        (Scale(0.1 / 3600) & Scale(0.1 / 3600))
        | SphericalToCartesian(wrap_lon_at=180)
        | rot
        | CartesianToSpherical(wrap_lon_at=360)
    )
    model.name = "v23tosky"
    return model


def _create_wcs(input_dm, shift_1=0, shift_2=0):
    """
    Create a basic WCS object (with optional shift) and
    append it to the input_dm.meta attribute.

    The WCS object will have a pipeline with the required
    steps to validate against the TweakReg pipeline.

    Parameters
    ----------
    input_dm : roman_datamodels.ImageModel
        A Roman image datamodel.
    shift_1 : int, optional
        Shift to be applied in the direction of the first axis, by default 0
    shift_2 : int, optional
        Shift to be applied in the direction of the second axis, by default 0
    """

    shape = input_dm.data.shape

    # create necessary transformations
    distortion = Shift(-shift_1) & Shift(-shift_2)
    distortion.bounding_box = ((-0.5, shape[-1] + 0.5), (-0.5, shape[-2] + 0.5))
    tel2sky = _create_tel2sky_model(input_dm)

    # create required frames
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(
        name="v2v3",
        axes_order=(0, 1),
        axes_names=("v2", "v3"),
        unit=(u.arcsec, u.arcsec),
    )
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name="world")

    # create pipeline
    pipeline = [
        wcs.Step(detector, distortion),
        wcs.Step(v2v3, tel2sky),
        wcs.Step(world, None),
    ]

    wcs_obj = wcs.WCS(pipeline)

    input_dm.meta["wcs"] = wcs_obj


@pytest.fixture
def base_image():
    """
    Create a base image with a realistic WCSInfo and a WCS.

    Notes
    -----
    The size of the image needs to be relatively large in order for
    the source catalog step to find a reasonable number of sources in the image.

    shift_1 and shift_2 (units in pixel) are used to shift the WCS projection plane.
    """

    def _base_image(shift_1=0, shift_2=0):
        l2 = maker_utils.mk_level2_image(shape=(100, 100))
        l2_im = rdm.ImageModel(l2)
        _create_wcs(l2_im)
        l2_im.meta.wcsinfo.vparity = -1
        return l2_im

    return _base_image


def create_asn_file(tmp_path, members_mapping=None):
    asn_content = """
        {
            "asn_type": "None",
            "asn_rule": "DMS_ELPP_Base",
            "version_id": null,
            "code_version": "0.9.1.dev28+ge987cc9.d20230106",
            "degraded_status": "No known degraded exposures in association.",
            "program": "noprogram",
            "constraints": "No constraints",
            "asn_id": "a3001",
            "target": "none",
            "asn_pool": "test_pool_name",
            "products": [
                {
                    "name": "files.asdf",
                    "members": [
                        {
                            "expname": "img_1.asdf",
                            "exptype": "science"
                        },
                        {
                            "expname": "img_2.asdf",
                            "exptype": "science"
                        }
                    ]
                }
            ]
        }
    """
    if members_mapping is not None:
        asn_dict = json.loads(asn_content)
        asn_dict["products"][0]["members"] = []
        for x in members_mapping:
            asn_dict["products"][0]["members"].append(x)
        asn_content = json.dumps(asn_dict)

    asn_file_path = str(tmp_path / "sample_asn.json")
    asn_file = StringIO()
    asn_file.write(asn_content)
    with open(asn_file_path, mode="w") as f:
        print(asn_file.getvalue(), file=f)

    return asn_file_path


@pytest.mark.parametrize(
    "input_models",
    [
        list(),
        [""],
        "",
        "testing",
        {},
    ],
)
def test_outlier_raises_error_on_invalid_input_models(input_models):
    """Test that OutlierDetection raises an error if input is not a ModelContainer."""

    step = OutlierDetectionStep()

    with pytest.raises(Exception) as exec_info:
        step.process(input_models)

    assert exec_info.type == TypeError


def test_outlier_skips_step_on_invalid_number_of_elements_in_input(base_image):
    img = base_image()

    step = OutlierDetectionStep()
    step.input_models = ModelContainer([img])

    step._run_process()

    assert step.input_models[0].meta.cal_step.outlier_detection == "SKIPPED"
    assert step.skip


def test_outlier_skips_step_on_exposure_type_different_from_wfi_image(base_image):
    img_1 = base_image()
    img_1.meta.exposure.type = "WFI_PRISM"
    img_2 = base_image()
    img_2.meta.exposure.type = "WFI_PRISM"

    step = OutlierDetectionStep()
    step.input_models = ModelContainer([img_1, img_2])

    step._run_process()

    assert step.input_models[0].meta.cal_step.outlier_detection == "SKIPPED"
    assert step.skip


def test_outlier_valid_input_datamodels(tmp_path, base_image):
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    img_2 = base_image()
    img_2.meta.filename = "img_2.asdf"

    step = OutlierDetectionStep()
    step.input_models = ModelContainer([img_1, img_2])

    step._run_process()

    assert step.skip is False
    assert img_1.meta.cal_step.outlier_detection == "COMPLETE"
    assert img_2.meta.cal_step.outlier_detection == "COMPLETE"


def test_outlier_valid_input_asn(tmp_path, base_image):
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    img_1.save(tmp_path / "img_1.asdf")
    img_2 = base_image()
    img_1.meta.filename = "img_2.asdf"
    img_2.save(tmp_path / "img_2.asdf")

    asn_filepath = create_asn_file(tmp_path)

    step = OutlierDetectionStep()
    step.input_models = ModelContainer(asn_filepath)

    step._run_process()

    assert step.skip is False
    assert all(
        x.meta.cal_step.outlier_detection == "COMPLETE" for x in step.input_models
    )


@pytest.mark.parametrize(
    "pars",
    [
        {
            "weight_type": "exptime",
            "pixfrac": 1.0,
            "kernel": "square",
            "fillval": "INDEF",
            "nlow": 0,
            "nhigh": 0,
            "maskpt": 0.7,
            "grow": 1,
            "snr": "4.0 3.0",
            "scale": "0.5 0.4",
            "backg": 0.0,
            "kernel_size": "7 7",
            "save_intermediate_results": False,
            "resample_data": True,
            "good_bits": 0,
            "allowed_memory": None,
            "in_memory": False,
            "make_output_path": None,
            "resample_suffix": "i2d",
        },
        {
            "weight_type": "exptime",
            "save_intermediate_results": True,
            "make_output_path": None,
            "resample_suffix": "some_other_suffix",
        },
    ],
)
def test_outlier_init_default_parameters(pars, base_image):
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    input_models = ModelContainer([img_1])

    step = outlier_detection.OutlierDetection(input_models, **pars)

    assert step.input_models == input_models
    assert step.outlierpars == pars
    assert step.make_output_path == pars["make_output_path"]
    assert step.resample_suffix == f"_outlier_{pars['resample_suffix']}.asdf"


def test_outlier_do_detection(base_image):
    img_1 = base_image()
    img_1.meta.filename = "img_1.asdf"
    img_2 = base_image()
    img_2.meta.filename = "img_2.asdf"
    input_models = ModelContainer([img_1, img_2])

    pars = {
        "weight_type": "exptime",
        "pixfrac": 1.0,
        "kernel": "square",
        "fillval": "INDEF",
        "nlow": 0,
        "nhigh": 0,
        "maskpt": 0.7,
        "grow": 1,
        "snr": "4.0 3.0",
        "scale": "0.5 0.4",
        "backg": 0.0,
        "kernel_size": "7 7",
        "save_intermediate_results": False,
        "resample_data": True,
        "good_bits": 0,
        "allowed_memory": None,
        "in_memory": False,
        "make_output_path": OutlierDetectionStep().make_output_path,
        "resample_suffix": "i2d",
    }

    blot_path_1 = Path(os.getcwd()) / img_1.meta.filename.replace(".asdf", "_blot.asdf")
    blot_path_2 = Path(os.getcwd()) / img_2.meta.filename.replace(".asdf", "_blot.asdf")
    median_path = Path(os.getcwd()) / "img_median.asdf"
    output_path = Path(os.getcwd()) / f"img_outlier_{pars['resample_suffix']}.asdf"

    outlier_files_path = [
        blot_path_1,
        blot_path_2,
        median_path,
        output_path,
    ]

    detection_step = outlier_detection.OutlierDetection
    step = detection_step(input_models, **pars)

    step.do_detection()

    assert all(x.exists() for x in outlier_files_path)

    # clean up
    [os.remove(file_path) for file_path in outlier_files_path]

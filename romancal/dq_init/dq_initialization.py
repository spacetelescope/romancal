import copy
import functools
import logging

import numpy as np
from astropy.time import Time
from roman_datamodels.datamodels import FpsModel, RampModel, ScienceRawModel, TvacModel
from roman_datamodels.dqflags import pixel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def to_ramp_model(model):
    if isinstance(model, RampModel):
        return model

    if isinstance(model, (FpsModel | TvacModel)):
        science_raw = ScienceRawModel.create_fake_data(
            {
                **model,
                "data": model.data.value,
                "amp33": model.amp33.value,
            }
        )

        # work-around issue with roman_datamodels not copying the nested read_pattern
        science_raw.meta.exposure.read_pattern = copy.deepcopy(
            model.meta.exposure.read_pattern
        )

        # move statistics to avoid conflicts with meta.statistics
        if "statistics" in science_raw.meta:
            science_raw["extras"] = {
                "tvac": {"meta": {"statistics": science_raw.meta.pop("statistics")}}
            }

        # drop old tagged scalars
        science_raw.meta.file_date = Time(science_raw.meta.file_date)
        for key, value in science_raw.items():
            if hasattr(value, "_tag") and isinstance(value, str):
                new_value = str(value)
            else:
                continue
            *branches, leaf = key.split(".")
            functools.reduce(getattr, branches, science_raw)[leaf] = new_value
        return to_ramp_model(science_raw)

    ramp_model = RampModel.create_from_model(
        {
            **model,
            "data": None,
        }
    )

    ramp_model.meta.cal_step = {}
    for step_name in ramp_model.schema_info("required")["roman"]["meta"]["cal_step"][
        "required"
    ].info:
        ramp_model.meta.cal_step[step_name] = model.meta.get("cal_step", {}).get(
            step_name, "INCOMPLETE"
        )

    shape = model.data.shape
    ramp_model.data = model.data.astype(np.float32)
    ramp_model.pixeldq = np.zeros(shape[1:], dtype=np.uint32)

    # check if the input model has a resultantdq from SDF
    if hasattr(ramp_model, "resultantdq"):
        ramp_model.groupdq = ramp_model.pop("resultantdq")
    else:
        ramp_model.groupdq = np.zeros(shape, dtype=np.uint8)

    # check for exposure data_problem
    if isinstance(data_problem := ramp_model.meta.exposure.data_problem, bool):
        ramp_model.meta.exposure.data_problem = "True" if data_problem else None

    return ramp_model


def do_dqinit(model, mask, expand_gw_flagging=0):
    """Convert model to a RampModel and update DQ flags.

    Parameters
    ----------
    model : Roman data model, ScienceRawModel, FpsModel
        model for which to update DQ

    mask : MaskRefModel or None
        reference mask model to use to update model DQ

    expand_gw_flagging : int
        expand GW-flagged region by this many pixels

    Returns
    -------
    RampModel constructed from model, with updated DQ flags.
    """
    is_tvac = isinstance(model, (FpsModel | TvacModel))

    # Convert to RampModel
    output_model = to_ramp_model(model)

    # guide window range information
    x_start = int(output_model.meta.guide_star.window_xstart)
    x_stop = int(output_model.meta.guide_star.window_xstop)
    y_start = int(output_model.meta.guide_star.window_ystart)
    y_stop = int(output_model.meta.guide_star.window_ystop)

    npix = output_model.pixeldq.shape[0]
    if x_start >= 0 and x_start < npix and y_start >= 0 and y_start < npix:
        # set pixeldq array to GW_AFFECTED_DATA (2**4) for the given range
        output_model.pixeldq[:, x_start:x_stop] = pixel.GW_AFFECTED_DATA
        log.info(
            f"Flagging rows from: {x_start} to {x_stop} as affected by guide window read"
        )

        # expand guide window if requested
        if expand_gw_flagging > 0:
            x_start = np.clip(x_start - expand_gw_flagging, 0, npix)
            x_stop = np.clip(x_stop + expand_gw_flagging, 0, npix)
            y_start = np.clip(y_start - expand_gw_flagging, 0, npix)
            y_stop = np.clip(y_stop + expand_gw_flagging, 0, npix)

        output_model.pixeldq[y_start:y_stop, x_start:x_stop] |= (
            pixel.DO_NOT_USE | pixel.GW_AFFECTED_DATA
        )
    else:
        log.info(
            f"Invalid guide window location: {x_start}, {x_stop}, {y_start}, {y_stop}"
        )

    # the reference read has been subtracted from the science data
    # in the L1 files.  Add it back into the data.
    # the TVAC files are special and there the reference read was
    # already added back in
    reference_read = getattr(output_model, "reference_read", None)
    if reference_read is not None and not is_tvac:
        output_model.data += reference_read
        del output_model.reference_read
    reference_amp33 = getattr(output_model, "reference_amp33", None)
    if reference_amp33 is not None and not is_tvac:
        output_model.amp33 += reference_amp33
        del output_model.reference_amp33

    # If a data encoding offset was added to the data, remove it
    data_encoding_offset = getattr(model.meta.instrument, "data_encoding_offset", 0)
    output_model.data -= data_encoding_offset

    if mask is not None and output_model.pixeldq.shape == mask.dq.shape:
        output_model.pixeldq |= mask.dq
        output_model.meta.cal_step.dq_init = "COMPLETE"
    else:
        log.warning("Mask data array is None or not the same shape as the science data")
        log.warning("Mask is not updated and step is marked skipped")
        output_model.meta.cal_step.dq_init = "SKIPPED"

    output_model.border_ref_pix_right = output_model.data[:, :, -4:].copy()
    output_model.border_ref_pix_left = output_model.data[:, :, :4].copy()
    output_model.border_ref_pix_top = output_model.data[:, :4, :].copy()
    output_model.border_ref_pix_bottom = output_model.data[:, -4:, :].copy()

    output_model.dq_border_ref_pix_right = output_model.pixeldq[:, -4:].copy()
    output_model.dq_border_ref_pix_left = output_model.pixeldq[:, :4].copy()
    output_model.dq_border_ref_pix_top = output_model.pixeldq[:4, :].copy()
    output_model.dq_border_ref_pix_bottom = output_model.pixeldq[-4:, :].copy()

    return output_model

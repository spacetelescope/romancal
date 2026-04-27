import logging

import numpy as np
from roman_datamodels.datamodels import FpsModel, RampModel, ScienceRawModel, TvacModel
from roman_datamodels.dqflags import pixel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


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
    try:
        # note that this succeeds even for ScienceRawModels
        output_model = ScienceRawModel.from_tvac_raw(model)
    except ValueError:
        output_model = model

    # Convert to RampModel
    output_model = RampModel.from_science_raw(output_model)

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
    reference_read = getattr(model, "reference_read", None)
    if reference_read is not None and not is_tvac:
        output_model.data += reference_read
        del output_model.reference_read
    reference_amp33 = getattr(model, "reference_amp33", None)
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

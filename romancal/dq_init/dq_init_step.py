#! /usr/bin/env python

import roman_datamodels as rdm
from roman_datamodels.datamodels import RampModel
from roman_datamodels.dqflags import pixel

from romancal.dq_init import dq_initialization
from romancal.stpipe import RomanStep

__all__ = ["DQInitStep"]


class DQInitStep(RomanStep):
    """Initialize the Data Quality extension from the
    mask reference file.

    The dq_init step initializes the pixeldq attribute of the
    input datamodel using the MASK reference file.  For some
    Guiding and Image model types, initalize the dq attribute of
    the input model instead.  The dq attribute of the MASK model
    is bitwise OR'd with the pixeldq (or dq) attribute of the
    input model.
    """

    reference_file_types = ["mask"]

    def process(self, input):
        """Perform the dq_init calibration step

        Parameters
        ----------
        input : Roman datamodel
            input roman datamodel

        Returns
        -------
        output_model : Roman datamodel
            result roman datamodel
        """
        # Open datamodel
        input_model = rdm.open(input)

        # Convert to RampModel
        output_model = RampModel.from_science_raw(input_model)

        # guide window range information
        x_start = input_model.meta.guidestar.gw_window_xstart
        x_end = input_model.meta.guidestar.gw_window_xsize + x_start
        # set pixeldq array to GW_AFFECTED_DATA (2**4) for the given range
        output_model.pixeldq[int(x_start) : int(x_end), :] = pixel.GW_AFFECTED_DATA
        self.log.info(
            f"Flagging rows from: {x_start} to {x_end} as affected by guide window read"
        )

        # Get reference file path
        reference_file_name = self.get_reference_file(output_model, "mask")

        # Test for reference file
        if reference_file_name != "N/A" and reference_file_name is not None:
            # If there are mask files, perform dq step
            # Open the relevant reference files as datamodels
            reference_file_model = rdm.open(reference_file_name)
            self.log.debug(f"Using MASK ref file: {reference_file_name}")

            # Apply the DQ step, in place
            dq_initialization.do_dqinit(
                output_model,
                reference_file_model,
            )

            # copy original border reference file arrays (data and dq)
            # to their own attributes. they will also remain attached to
            # the science data until they are trimmed at ramp_fit
            # these arrays include the overlap regions in the corners

            output_model.border_ref_pix_right = output_model.data[:, :, -4:].copy()
            output_model.border_ref_pix_left = output_model.data[:, :, :4].copy()
            output_model.border_ref_pix_top = output_model.data[:, :4, :].copy()
            output_model.border_ref_pix_bottom = output_model.data[:, -4:, :].copy()

            output_model.dq_border_ref_pix_right = output_model.pixeldq[:, -4:].copy()
            output_model.dq_border_ref_pix_left = output_model.pixeldq[:, :4].copy()
            output_model.dq_border_ref_pix_top = output_model.pixeldq[:4, :].copy()
            output_model.dq_border_ref_pix_bottom = output_model.pixeldq[-4:, :].copy()

        else:
            # Skip DQ step if no mask files
            reference_file_model = None
            self.log.warning("No MASK reference file found.")
            self.log.warning("DQ initialization step will be skipped.")

            output_model.meta.cal_step.dq_init = "SKIPPED"

        # Close the input and reference files
        input_model.close()

        try:
            reference_file_model.close()
        except AttributeError:
            pass

        if self.save_results:
            try:
                self.suffix = "dqinit"
            except AttributeError:
                self["suffix"] = "dqinit"

        return output_model

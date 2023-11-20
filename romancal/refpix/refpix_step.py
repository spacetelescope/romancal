import logging

import roman_datamodels as rdm

from romancal.refpix import refpix
from romancal.stpipe import RomanStep

__all__ = ["RefPixStep"]

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class RefPixStep(RomanStep):
    """
    RefpixStep: Module for correcting the the science data using the reference
        pixels
    """

    spec = """
    remove_offset = boolean(default=True) # Turn on or off removing the data offset
    # prior to the reference pixel correction, then returning the offset afterwords.
    remove_trends = boolean(default=True) # Turn on or off removing the boundary
    # linear trends
    cosine_interpolate = boolean(default=True) # Turn on or off the cosine
    # interpolation of the reference pixels
    fft_interpolate = boolean(default=True) # Turn on or off the FFT interpolation
    # of the reference pixel padded values.
    """

    reference_file_types = ["refpix"]

    def process(self, input):
        """
        Perform the reference pixel correction
        """

        # open the input data model
        log.debug(f"Opening the science data: {input}")
        datamodel = rdm.open(input)

        # Get the reference file
        ref_file = self.get_reference_file(datamodel, "refpix")

        # Check for a valid reference file
        if ref_file == "N/A":
            self.log.warning("No REFPIX reference file found")
            self.log.warning("Reference pixel correction step will be skipped")
            result = datamodel
            result.meta.cal_step.refpix = "SKIPPED"
            return result

        log.debug(f"Opening the reference file: {ref_file}")
        with rdm.open(ref_file) as refs:
            # Run the correction
            log.debug("Running the reference pixel correction")
            output = refpix.run_steps(
                datamodel,
                refs,
                self.remove_offset,
                self.remove_trends,
                self.cosine_interpolate,
                self.fft_interpolate,
            )
            # Update the step status
            output.meta.cal_step["refpix"] = "COMPLETE"
            if self.save_results:
                try:
                    self.suffix = "refpix"
                except AttributeError:
                    self["suffix"] = "refpix"
        return output

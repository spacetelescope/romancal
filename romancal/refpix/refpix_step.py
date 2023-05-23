import roman_datamodels as rdm

from romancal.refpix import refpix
from romancal.stpipe import RomanStep

__all__ = ["RefpixStep"]


class RefpixStep(RomanStep):
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
    fft_interpolate = boolean(default=False) # Turn on or off the FFT interpolation
    # of the reference pixel padded values.
    """

    reference_file_types = ["refpix"]

    def process(self, input):
        """
        Perform the reference pixel correction
        """

        # open the input data model
        with rdm.open(input) as datamodel:

            # Get the reference file
            ref_file = self.get_reference_file(datamodel, "refpix")

            with rdm.open(ref_file) as refs:

                # Run the correction
                control = refpix.Control.from_step(self)
                return refpix.run_steps(datamodel, refs, control)

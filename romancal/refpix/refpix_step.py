import roman_datamodels as rdm

from romancal.refpix import refpix
from romancal.stpipe import RomanStep

__all__ = ["RefpixStep"]


class RefpixStep(RomanStep):
    """
    RefpixStep: Module for correcting the the science data using the reference
        pixels
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

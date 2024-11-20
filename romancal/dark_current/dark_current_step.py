#! /usr/bin/env python

from roman_datamodels import datamodels as rdm

from romancal.stpipe import RomanStep

__all__ = ["DarkCurrentStep"]


class DarkCurrentStep(RomanStep):
    """
    DarkCurrentStep: Performs dark current correction by subtracting
    dark current reference data from the input science data model.
    """

    class_alias = "dark"

    spec = """
        dark_output = output_file(default = None) # Dark corrected model
    """

    reference_file_types = ["dark"]

    def process(self, input):
        if isinstance(input, rdm.DataModel):
            input_model = input
        else:
            # Open the input data model
            input_model = rdm.open(input)

        # Get the name of the dark reference file to use
        self.dark_name = self.get_reference_file(input_model, "dark")
        # Check for a valid reference file
        if self.dark_name == "N/A":
            self.log.warning("No DARK reference file found")
            self.log.warning("Dark current step will be skipped")
            result = input_model
            result.meta.cal_step.dark = "SKIPPED"
            return result

        self.log.info("Using DARK reference file: %s", self.dark_name)

        # Open dark model
        with rdm.open(self.dark_name) as dark_model:
            # Temporary patch to utilize stcal dark step until MA table support
            # is fully implemented
            if "nresultants" not in dark_model.meta.exposure:
                dark_model.meta.exposure["nresultants"] = dark_model.data.shape[0]

            # Do the dark correction
            out_model = input_model
            nresultants = len(input_model.meta.exposure["read_pattern"])
            out_model.data -= dark_model.data[:nresultants].value
            out_model.pixeldq |= dark_model.dq
            out_model.meta.cal_step.dark = "COMPLETE"

            # Save dark data to file
            if self.dark_output is not None:
                dark_model.save(self.dark_output)
                # not clear to me that this makes any sense for Roman

        if self.save_results:
            try:
                self.suffix = "darkcurrent"
            except AttributeError:
                self["suffix"] = "darkcurrent"

        return out_model

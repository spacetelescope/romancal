#! /usr/bin/env python

from roman_datamodels import datamodels as rdm

from romancal.stpipe import RomanStep

__all__ = ["DarkCurrentStep"]


class DarkCurrentStep(RomanStep):
    """
    DarkCurrentStep: Performs dark current correction by subtracting
    dark current reference data from the input science data model.
    """

    spec = """
        dark_output = output_file(default = None) # Dark corrected model
    """

    reference_file_types = ["dark"]

    def process(self, input):
        # Open the input data model
        with rdm.open(input, lazy_load=False) as input_model:
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
            dark_model = rdm.open(self.dark_name)

            # Temporary patch to utilize stcal dark step until MA table support
            # is fully implemented
            if "ngroups" not in dark_model.meta.exposure:
                dark_model.meta.exposure["ngroups"] = dark_model.data.shape[0]
            if "nframes" not in dark_model.meta.exposure:
                dark_model.meta.exposure["nframes"] = input_model.meta.exposure.nframes
            if "groupgap" not in dark_model.meta.exposure:
                dark_model.meta.exposure["groupgap"] = (
                    input_model.meta.exposure.groupgap
                )

            # Do the dark correction
            out_model = input_model
            nresultants = len(input_model.meta.exposure["read_pattern"])
            out_model.data -= dark_model.data[:nresultants]
            out_model.pixeldq |= dark_model.dq
            out_model.meta.cal_step.dark = "COMPLETE"

            # Save dark data to file
            if self.dark_output is not None:
                dark_model.save(self.dark_output)
                # not clear to me that this makes any sense for Roman
            dark_model.close()

        if self.save_results:
            try:
                self.suffix = "darkcurrent"
            except AttributeError:
                self["suffix"] = "darkcurrent"
        return out_model

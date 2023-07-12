"""
Detect and flag outlier in a science image
"""

import roman_datamodels as rdm

from ..stpipe import RomanStep

__all__ = ["OutlierDetectionStep"]


class OutlierDetectionStep(RomanStep):
    """Detect and flag outliers in a science image."""

    def process(self, input):
        input_model = rdm.open(input, lazy_load=False)

        # No reference files
        # reference_file_model = {}
        # Do the outlier detection
        # output_model = outlier_detection.do_correction(
        #     input_model,
        #     reference_file_model,
        # )

        output_model = input_model

        # Close the input and reference files
        input_model.close()

        if self.save_results:
            try:
                self.suffix = "outlier_detection"
            except AttributeError:
                self["suffix"] = "outlier_detection"

        return output_model

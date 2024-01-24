"""Public common step definition for OutlierDetection processing."""

import os
from functools import partial

from romancal.datamodels import ModelContainer
from romancal.outlier_detection import outlier_detection

from ..stpipe import RomanStep

__all__ = ["OutlierDetectionStep"]


class OutlierDetectionStep(RomanStep):
    """Flag outlier bad pixels and cosmic rays in DQ array of each input image.

    Input images can be listed in an input association file or already wrapped
    with a ModelContainer. DQ arrays are modified in place.

    Parameters
    -----------
    input_data : `~romancal.datamodels.container.ModelContainer`
        A `~romancal.datamodels.container.ModelContainer` object.

    """

    class_alias = "outlier_detection"

    # The members of spec needs to be a super-set of all parameters needed
    # by the various versions of the outlier_detection algorithms, and each
    # version will pick and choose what they need while ignoring the rest.
    spec = """
        weight_type = option('ivm','exptime',default='ivm') # Weighting type to use to create the median image
        pixfrac = float(default=1.0) # Fraction by which input pixels are shrunk before being drizzled onto the output image grid
        kernel = string(default='square') # Shape of the kernel used for flux distribution onto output images
        fillval = string(default='INDEF') # Value assigned to output pixels that have zero weight or no flux during drizzling
        nlow = integer(default=0) # The number of low values in each pixel stack to ignore when computing the median value
        nhigh = integer(default=0) # The number of high values in each pixel stack to ignore when computing the median value
        maskpt = float(default=0.7) # Percentage of weight image values below which they are flagged as bad pixels
        grow = integer(default=1) # The distance beyond the rejection limit for additional pixels to be rejected in an image
        snr = string(default='5.0 4.0') # The signal-to-noise values to use for bad pixel identification
        scale = string(default='1.2 0.7') # The scaling factor applied to derivative used to identify bad pixels
        backg = float(default=0.0) # User-specified background value to subtract during final identification step
        kernel_size = string(default='7 7') # Size of kernel to be used during resampling of the data
        save_intermediate_results = boolean(default=False) # Specifies whether or not to write out intermediate products to disk
        resample_data = boolean(default=True) # Specifies whether or not to resample the input images when performing outlier detection
        good_bits = string(default="0")  # DQ bit value to be considered 'good'
        allowed_memory = float(default=None)  # Fraction of memory to use for the combined image
        in_memory = boolean(default=False) # Specifies whether or not to keep all intermediate products and datamodels in memory
    """  # noqa: E501

    def process(self, input_models):
        """Perform outlier detection processing on input data."""

        self.input_models = input_models
        self.input_container = isinstance(self.input_models, ModelContainer)
        self.skip = False

        # validation
        if self.input_container and (
            len(self.input_models) >= 2
            and all(
                model.meta.exposure.type == "WFI_IMAGE" for model in self.input_models
            )
        ):
            # Setup output path naming if associations are involved.
            asn_id = self.input_models.asn_table.get("asn_id", None)
            if asn_id is not None:
                _make_output_path = self.search_attr(
                    "_make_output_path", parent_first=True
                )
                self._make_output_path = partial(_make_output_path, asn_id=asn_id)

            detection_step = outlier_detection.OutlierDetection
            pars = {
                "weight_type": self.weight_type,
                "pixfrac": self.pixfrac,
                "kernel": self.kernel,
                "fillval": self.fillval,
                "nlow": self.nlow,
                "nhigh": self.nhigh,
                "maskpt": self.maskpt,
                "grow": self.grow,
                "snr": self.snr,
                "scale": self.scale,
                "backg": self.backg,
                "kernel_size": self.kernel_size,
                "save_intermediate_results": self.save_intermediate_results,
                "resample_data": self.resample_data,
                "good_bits": self.good_bits,
                "allowed_memory": self.allowed_memory,
                "in_memory": self.in_memory,
                "make_output_path": self.make_output_path,
                "resample_suffix": "i2d",
            }

            self.log.debug(
                f"Using {detection_step.__name__} class for outlier_detection"
            )

            # Set up outlier detection, then do detection
            step = detection_step(self.input_models, **pars)
            step.do_detection()

            state = "COMPLETE"

            if not self.save_intermediate_results:
                self.log.debug(
                    "The following files will be deleted since \
                    save_intermediate_results=False:"
                )
            for model in self.input_models:
                model.meta.cal_step["outlier_detection"] = state
                if not self.save_intermediate_results:
                    #  Remove unwanted files
                    crf_path = self.make_output_path(basepath=model.meta.filename)
                    # These lines to be used when/if outlier_i2d files follow
                    # output_dir crf_file = os.path.basename(crf_path)
                    # outlr_path = crf_path.replace(crf_file, outlr_file)
                    outlr_file = model.meta.filename.replace("cal", "outlier_i2d")
                    blot_path = crf_path.replace("crf", "blot")
                    median_path = blot_path.replace("blot", "median")

                    for fle in [outlr_file, blot_path, median_path]:
                        if os.path.isfile(fle):
                            os.remove(fle)
                            self.log.debug(f"    {fle}")

            return self.input_models

        # if input is not valid then log a warning message and
        # skip outlier detection step
        self.log.warning(
            "Input is not a ModelContainer, \
            does not contains >= 2 elements, \
            or the elements contain the wrong exposure type \
            (i.e. meta.exposure.type != 'WFI_IMAGE')."
        )
        if self.input_container:
            for model in self.input_models:
                model.meta.cal_step["outlier_detection"] = "SKIPPED"
        self.skip = True
        return self.input_models

"""Public common step definition for OutlierDetection processing."""

import contextlib
import os
from functools import partial

from romancal.datamodels import ModelContainer
from romancal.outlier_detection import outlier_detection

from ..stpipe import RomanStep

__all__ = ["OutlierDetectionStep"]


class OutlierDetectionStep(RomanStep):
    """Flag outlier bad pixels and cosmic rays in DQ array of each input image.

    Input images can be listed in an input association file or already opened
    with a ModelContainer.  DQ arrays are modified in place.

    Parameters
    -----------
    input_data : asn file or ~romancal.datamodels.ModelContainer
        Single filename association table, or a datamodels.ModelContainer.

    """

    class_alias = "outlier_detection"

    # The members of spec needs to be a super-set of all parameters needed
    # by the various versions of the outlier_detection algorithms, and each
    # version will pick and choose what they need while ignoring the rest.
    spec = """
        weight_type = option('ivm','exptime',default='ivm')
        pixfrac = float(default=1.0)
        kernel = string(default='square') # drizzle kernel
        fillval = string(default='INDEF')
        nlow = integer(default=0)
        nhigh = integer(default=0)
        maskpt = float(default=0.7)
        grow = integer(default=1)
        snr = string(default='5.0 4.0')
        scale = string(default='1.2 0.7')
        backg = float(default=0.0)
        kernel_size = string(default='7 7')
        threshold_percent = float(default=99.8)
        ifu_second_check = boolean(default=False)
        save_intermediate_results = boolean(default=False)
        resample_data = boolean(default=True)
        good_bits = string(default="~DO_NOT_USE")  # DQ flags to allow
        scale_detection = boolean(default=False)
        search_output_file = boolean(default=False)
        allowed_memory = float(default=None)  # Fraction of memory to use for the combined image
        in_memory = boolean(default=False)
    """

    def process(self, input_models):
        """Perform outlier detection processing on input data."""

        self.input_models = input_models
        self.input_container = isinstance(self.input_models, ModelContainer)
        # Setup output path naming if associations are involved.
        asn_id = None
        with contextlib.suppress(AttributeError, KeyError):
            asn_id = self.input_models.meta.asn_table.asn_id
        if asn_id is None:
            asn_id = self.search_attr("asn_id")
        if asn_id is not None:
            # TODO: revisit this once we have a better idea of the associations
            _make_output_path = self.search_attr("_make_output_path", parent_first=True)

            self._make_output_path = partial(_make_output_path, asn_id=asn_id)

        # Setup outlier detection parameters
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
            "threshold_percent": self.threshold_percent,
            "ifu_second_check": self.ifu_second_check,
            "allowed_memory": self.allowed_memory,
            "in_memory": self.in_memory,
            "save_intermediate_results": self.save_intermediate_results,
            "resample_data": self.resample_data,
            "good_bits": self.good_bits,
            "make_output_path": self.make_output_path,
        }

        # Add logic here to select which version of OutlierDetection
        # needs to be used depending on the input data
        if self.input_container:
            single_model = self.input_models[0]
        else:
            single_model = self.input_models
        exptype = single_model.meta.exposure.type
        self.check_input()

        if exptype == "WFI_IMAGE":
            # imaging with resampling
            detection_step = outlier_detection.OutlierDetection
            pars["resample_suffix"] = "i2d"
        else:
            self.log.error(
                "Outlier detection failed for unknown/unsupported ",
                f"exposure type: {exptype}",
            )
            self.valid_input = False

        if not self.valid_input:
            if self.input_container:
                for model in self.input_models:
                    model.meta.cal_step.outlier_detection = "SKIPPED"
            else:
                self.input_models.meta.cal_step.outlier_detection = "SKIPPED"
            self.skip = True
            return self.input_models

        self.log.debug(f"Using {detection_step.__name__} class for outlier_detection")

        # Set up outlier detection, then do detection
        step = detection_step(self.input_models, **pars)
        step.do_detection()

        state = "COMPLETE"
        if self.input_container:
            if not self.save_intermediate_results:
                self.log.debug(
                    "The following files will be deleted since save_intermediate_results=False:"
                )
            for model in self.input_models:
                model.meta.cal_step.outlier_detection = state
                if not self.save_intermediate_results:
                    #  Remove unwanted files
                    crf_path = self.make_output_path(basepath=model.meta.filename)
                    #  These lines to be used when/if outlier_i2d files follow output_dir
                    #  crf_file = os.path.basename(crf_path)
                    #  outlr_path = crf_path.replace(crf_file, outlr_file)
                    outlr_file = model.meta.filename.replace("cal", "outlier_i2d")
                    blot_path = crf_path.replace("crf", "blot")
                    median_path = blot_path.replace("blot", "median")

                    for fle in [outlr_file, blot_path, median_path]:
                        if os.path.isfile(fle):
                            os.remove(fle)
                            self.log.debug(f"    {fle}")
        else:
            self.input_models.meta.cal_step.outlier_detection = state
        return self.input_models

    def check_input(self):
        """Use this method to determine whether input is valid or not."""
        if self.input_container:
            ninputs = len(self.input_models)
            if not isinstance(self.input_models, ModelContainer):
                self.log.warning("Input is not a ModelContainer")
                self.log.warning("Outlier detection step will be skipped")
                self.valid_input = False
            elif ninputs < 2:
                self.log.warning(f"Input only contains {ninputs} exposure")
                self.log.warning("Outlier detection step will be skipped")
                self.valid_input = False
            else:
                self.valid_input = True
                self.log.info(f"Performing outlier detection on {ninputs} inputs")

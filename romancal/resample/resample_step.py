from __future__ import annotations

import logging
import os
from copy import deepcopy
from typing import TYPE_CHECKING

import asdf
import numpy as np
from roman_datamodels import datamodels

from ..datamodels import ModelLibrary
from ..stpipe import RomanStep
from .resample import ResampleData

if TYPE_CHECKING:
    from typing import ClassVar

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["ResampleStep"]

# conversion factor from steradian to squared arcsec
SR_TO_ARCSEC2 = 4.254517e10


class ResampleStep(RomanStep):
    """
    Resample input data onto a regular grid using the drizzle algorithm.

    .. note::
        When supplied via ``output_wcs``, a custom WCS overrides other custom
        WCS parameters such as ``output_shape`` (now computed from by
        ``output_wcs.bounding_box``), ``crpix``

    Parameters
    -----------
    input : str, `roman_datamodels.datamodels.DataModel`, or `~romancal.datamodels.ModelLibrary`
        If a string is provided, it should correspond to either a single ASDF filename
        or an association filename. Alternatively, a single DataModel instance can be
        provided instead of an ASDF filename. Multiple files can be processed via
        either an association file or wrapped by a
        `~romancal.datamodels.ModelLibrary`.

    Returns
    -------
    : `roman_datamodels.datamodels.MosaicModel`
        A mosaic datamodel with the final output frame.
    """

    class_alias = "resample"

    spec = """
        pixfrac = float(default=1.0)
        kernel = string(default='square')
        fillval = string(default='NAN' )
        weight_type = option('ivm', 'exptime', None, default='ivm')
        output_shape = int_list(min=2, max=2, default=None)  # [x, y] order
        crpix = float_list(min=2, max=2, default=None)
        crval = float_list(min=2, max=2, default=None)
        rotation = float(default=None)
        pixel_scale_ratio = float(default=1.0) # Ratio of output to input pixel scale
        pixel_scale = float(default=None) # Absolute pixel scale in arcsec
        output_wcs = string(default='')  # Custom output WCS.
        resample_on_skycell = boolean(default=True)  # if association contains skycell information use it for the wcs
        in_memory = boolean(default=True)
        good_bits = string(default='~DO_NOT_USE+NON_SCIENCE')  # The good bits to use for building the resampling mask.
    """

    reference_file_types: ClassVar = []

    def process(self, input):
        # There is no way to check for minimum values in output_shape
        # within the step spec so check them here.
        if self.output_shape is not None:
            for v in self.output_shape:
                if v < 1:
                    raise ValueError(
                        f"output shape values must be >= 1: {self.output_shape}"
                    )

        if isinstance(input, datamodels.DataModel):
            input_models = ModelLibrary([input])
            # set output filename from meta.filename found in the first datamodel
            output_filename = input.meta.filename
        elif isinstance(input, str):
            # either a single asdf filename or an association filename
            try:
                # association filename
                input_models = ModelLibrary(input)
            except Exception:
                # single ASDF filename
                input_models = ModelLibrary([input])
            output_filename = input_models.asn["products"][0]["name"]
        elif isinstance(input, ModelLibrary):
            input_models = input
            if "name" in input_models.asn["products"][0]:
                output_filename = input_models.asn["products"][0]["name"]
            else:
                # set output filename using the common prefix of all datamodels
                output_filename = f"{os.path.commonprefix([x['expname'] for x in input_models.asn['products'][0]['members']])}.asdf"
                if len(output_filename) == 0:
                    # set default filename if no common prefix can be determined
                    output_filename = "resample_output.asdf"
        else:
            raise TypeError(
                "Input must be an ASN filename, a ModelLibrary, "
                "a single ASDF filename, or a single Roman DataModel."
            )

        if not len(input_models):
            raise ValueError("At least 1 file must be provided")

        # Issue a warning about the use of exptime weighting
        if self.weight_type == "exptime":
            self.log.warning("Use of EXPTIME weighting will result in incorrect")
            self.log.warning("propagated errors in the resampled product")

        output_wcs = self._load_custom_wcs(self.output_wcs, self.output_shape)

        if output_wcs is None:
            wcs_kwargs = {
                "pscale_ratio": self.pixel_scale_ratio,
                "pscale": self.pixel_scale,
                "rotation": self.rotation,
                "shape": None if self.output_shape is None else self.output_shape[::-1],
                "crpix": self.crpix,
                "crval": self.crval,
            }
        else:
            wcs_kwargs = None

        # Call the resampling routine
        resamp = ResampleData(
            input_models,
            output_wcs,
            self.pixfrac,
            self.kernel,
            self.fillval,
            self.weight_type,
            self.good_bits,
            True,
            True,
            "from_var",
            True,
            True,
            self.resample_on_skycell,
            wcs_kwargs,
        )
        result = resamp.resample_group(range(len(input_models)))
        result.meta.filename = output_filename

        self._final_updates(result, input_models)

        return result

    def _final_updates(self, model, input_models):
        model.meta.cal_step["resample"] = "COMPLETE"

        # if pixel_scale exists, it will override pixel_scale_ratio.
        # calculate the actual value of pixel_scale_ratio based on pixel_scale
        # because source_catalog uses this value from the header.
        model.meta.resample.pixel_scale_ratio = (
            self.pixel_scale / np.sqrt(model.meta.photometry.pixel_area * SR_TO_ARCSEC2)
            if self.pixel_scale
            else self.pixel_scale_ratio
        )
        model.meta.resample.pixfrac = self.pixfrac
        if model.meta.photometry.pixel_area is not None:
            model.meta.photometry.pixel_area *= model.meta.resample.pixel_scale_ratio**2
        model.meta.resample["good_bits"] = self.good_bits

    @staticmethod
    def _load_custom_wcs(asdf_wcs_file, output_shape):
        if not asdf_wcs_file:
            return None

        with asdf.open(asdf_wcs_file) as af:
            wcs = deepcopy(af.tree["wcs"])

        if output_shape is not None:
            wcs.array_shape = output_shape[::-1]
        elif wcs.pixel_shape is not None:
            wcs.array_shape = wcs.pixel_shape[::-1]
        elif wcs.bounding_box is not None:
            wcs.array_shape = tuple(
                int(axs[1] - axs[0] + 0.5)
                for axs in wcs.bounding_box.bounding_box(order="C")
            )
        elif wcs.array_shape is None:
            raise ValueError(
                "Step argument 'output_shape' is required when custom WCS "
                "does not have neither of 'array_shape', 'pixel_shape', or "
                "'bounding_box' attributes set."
            )

        return {"wcs": wcs}

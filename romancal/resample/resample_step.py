import logging
import os
from copy import deepcopy

import asdf
import numpy as np
from roman_datamodels import datamodels
from stcal.alignment import util
from stpipe.extern.configobj.configobj import ConfigObj
from stpipe.extern.configobj.validate import Validator

from ..datamodels import ModelContainer
from ..stpipe import RomanStep
from . import resample

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["ResampleStep"]


# Force use of all DQ flagged data except for DO_NOT_USE and NON_SCIENCE
GOOD_BITS = "~DO_NOT_USE+NON_SCIENCE"


class ResampleStep(RomanStep):
    """
    Resample input data onto a regular grid using the drizzle algorithm.

    .. note::
        When supplied via ``output_wcs``, a custom WCS overrides other custom
        WCS parameters such as ``output_shape`` (now computed from by
        ``output_wcs.bounding_box``), ``crpix``

    Parameters
    -----------
    input : str, `roman_datamodels.datamodels.DataModel`, or `~romancal.datamodels.container.ModelContainer`
        If a string is provided, it should correspond to either a single ASDF filename
        or an association filename. Alternatively, a single DataModel instance can be
        provided instead of an ASDF filename. Multiple files can be processed via
        either an association file or wrapped by a
        `~romancal.datamodels.container.ModelContainer`.

    Returns
    -------
    : `roman_datamodels.datamodels.MosaicModel`
        A mosaic datamodel with the final output frame.
    """  # noqa: E501

    class_alias = "resample"

    spec = """
        pixfrac = float(default=1.0) # change back to None when drizpar reference files are updated
        kernel = string(default='square') # change back to None when drizpar reference files are updated
        fillval = string(default='INDEF' ) # change back to None when drizpar reference files are updated
        weight_type = option('ivm', 'exptime', None, default='ivm')  # change back to None when drizpar ref update
        output_shape = int_list(min=2, max=2, default=None)  # [x, y] order
        crpix = float_list(min=2, max=2, default=None)
        crval = float_list(min=2, max=2, default=None)
        rotation = float(default=None)
        pixel_scale_ratio = float(default=1.0) # Ratio of input to output pixel scale
        pixel_scale = float(default=None) # Absolute pixel scale in arcsec
        output_wcs = string(default='')  # Custom output WCS.
        single = boolean(default=False)
        blendheaders = boolean(default=True)
        allowed_memory = float(default=None)  # Fraction of memory to use for the combined image.
        in_memory = boolean(default=True)
    """  # noqa: E501

    reference_file_types = []

    def process(self, input):
        if isinstance(input, datamodels.DataModel):
            input_models = ModelContainer([input])
            # set output filename from meta.filename found in the first datamodel
            output = input_models[0].meta.filename
            self.blendheaders = False
        elif isinstance(input, str):
            # either a single asdf filename or an association filename
            try:
                # association filename
                input_models = ModelContainer(input)
            except Exception:
                # single ASDF filename
                input_models = ModelContainer([input])
            if hasattr(input_models, "asn_table") and len(input_models.asn_table):
                # set output filename from ASN table
                output = input_models.asn_table["products"][0]["name"]
            elif hasattr(input_models[0], "meta"):
                # set output filename from meta.filename found in the first datamodel
                output = input_models[0].meta.filename
        elif isinstance(input, ModelContainer):
            input_models = input
            # set output filename using the common prefix of all datamodels
            output = (
                f"{os.path.commonprefix([x.meta.filename for x in input_models])}.asdf"
            )
            if len(output) == 0:
                # set default filename if no common prefix can be determined
                output = "resample_output.asdf"
        else:
            raise TypeError(
                "Input must be an ASN filename, a ModelContainer, "
                "a single ASDF filename, or a single Roman DataModel."
            )

        # Check that input models are 2D images
        if len(input_models[0].data.shape) != 2:
            # resample can only handle 2D images, not 3D cubes, etc
            raise RuntimeError(f"Input {input_models[0]} is not a 2D image.")

        self.wht_type = self.weight_type
        self.log.info("Setting drizzle's default parameters...")
        kwargs = self.set_drizzle_defaults()
        kwargs["allowed_memory"] = self.allowed_memory

        # Issue a warning about the use of exptime weighting
        if self.wht_type == "exptime":
            self.log.warning("Use of EXPTIME weighting will result in incorrect")
            self.log.warning("propagated errors in the resampled product")

        # Custom output WCS parameters.
        # Modify get_drizpars if any of these get into reference files:
        kwargs["output_shape"] = self._check_list_pars(
            self.output_shape, "output_shape", min_vals=[1, 1]
        )
        kwargs["output_wcs"] = self._load_custom_wcs(
            self.output_wcs, kwargs["output_shape"]
        )
        kwargs["crpix"] = self._check_list_pars(self.crpix, "crpix")
        kwargs["crval"] = self._check_list_pars(self.crval, "crval")
        kwargs["rotation"] = self.rotation
        kwargs["pscale"] = self.pixel_scale
        kwargs["pscale_ratio"] = self.pixel_scale_ratio
        kwargs["in_memory"] = self.in_memory

        # Call the resampling routine
        resamp = resample.ResampleData(input_models, output=output, **kwargs)
        result = resamp.do_drizzle()

        for model in result:
            self._final_updates(model, input_models, kwargs)
        if len(result) == 1:
            result = result[0]

        input_models.close()
        return result

    def _final_updates(self, model, input_models, kwargs):
        model.meta.cal_step["resample"] = "COMPLETE"
        util.update_s_region_imaging(model)
        if (
            input_models.asn_pool_name is not None
            and input_models.asn_table_name is not None
        ):
            # update ASN attributes
            model.meta.asn.pool_name = input_models.asn_pool_name
            model.meta.asn.table_name = input_models.asn_table_name

        # if pixel_scale exists, it will override pixel_scale_ratio.
        # calculate the actual value of pixel_scale_ratio based on pixel_scale
        # because source_catalog uses this value from the header.
        model.meta.resample.pixel_scale_ratio = (
            self.pixel_scale / np.sqrt(model.meta.photometry.pixelarea_arcsecsq)
            if self.pixel_scale
            else self.pixel_scale_ratio
        )
        model.meta.resample.pixfrac = kwargs["pixfrac"]
        self.update_phot_keywords(model)

    @staticmethod
    def _check_list_pars(vals, name, min_vals=None):
        """
        Check if a specific keyword parameter is properly formatted.

        Parameters
        ----------
        vals : list or tuple
            A list or tuple containing a pair of values currently assigned to the
            keyword parameter `name`. Both values must be either `None` or not `None`.
        name : str
            The name of the keyword parameter.
        min_vals : list or tuple, optional
            A list or tuple containing a pair of minimum values to be assigned
            to `name`, by default None.

        Returns
        -------
        None or list
            If either `vals` is set to `None` (or both of its elements), the
            returned result will be `None`. Otherwise, the returned result will be
            a list containing the current values assigned to `name`.

        Raises
        ------
        ValueError
            This error will be raised if any of the following conditions are met:
            - the number of elements of `vals` is not 2;
            - the currently assigned values of `vals` are smaller than the
            minimum value provided;
            - one element is `None` and the other is not `None`.
        """
        if vals is None:
            return None
        if len(vals) != 2:
            raise ValueError(f"List '{name}' must have exactly two elements.")
        n = sum(x is None for x in vals)
        if n == 2:
            return None
        elif n == 0:
            if min_vals and sum(x >= y for x, y in zip(vals, min_vals)) != 2:
                raise ValueError(
                    f"'{name}' values must be larger or equal to {list(min_vals)}"
                )
            return list(vals)
        else:
            raise ValueError(f"Both '{name}' values must be either None or not None.")

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

        return wcs

    def update_phot_keywords(self, model):
        """Update pixel scale keywords"""
        if model.meta.photometry.pixelarea_steradians is not None:
            model.meta.photometry.pixelarea_steradians *= (
                model.meta.resample.pixel_scale_ratio**2
            )
        if model.meta.photometry.pixelarea_arcsecsq is not None:
            model.meta.photometry.pixelarea_arcsecsq *= (
                model.meta.resample.pixel_scale_ratio**2
            )

    def set_drizzle_defaults(self):
        """Set the default parameters for drizzle."""
        configspec = self.load_spec_file()
        config = ConfigObj(configspec=configspec)
        if config.validate(Validator()):
            kwargs = config.dict()

        if self.pixfrac is None:
            self.pixfrac = 1.0
        if self.kernel is None:
            self.kernel = "square"
        if self.fillval is None:
            self.fillval = "INDEF"
        # Force definition of good bits
        kwargs["good_bits"] = GOOD_BITS
        kwargs["pixfrac"] = self.pixfrac
        kwargs["kernel"] = str(self.kernel)
        kwargs["fillval"] = str(self.fillval)
        #  self.weight_type has a default value of None
        # The other instruments read this parameter from a reference file
        if self.wht_type is None:
            self.wht_type = "ivm"

        kwargs["wht_type"] = str(self.wht_type)
        kwargs["pscale_ratio"] = self.pixel_scale_ratio
        kwargs.pop("pixel_scale_ratio")

        for k, v in kwargs.items():
            if k in [
                "pixfrac",
                "kernel",
                "fillval",
                "wht_type",
                "pscale_ratio",
            ]:
                log.info("  using: %s=%s", k, repr(v))

        return kwargs

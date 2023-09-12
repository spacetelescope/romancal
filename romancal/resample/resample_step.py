import logging
import os
import re
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
    input : str, ~roman_datamodels.datamodels.DataModel, or ~romancal.datamodels.ModelContainer
        If a string is provided, it should correspond to either a single ASDF filename
        or an association filename. Alternatively, a single DataModel instance can be
        provided instead of an ASDF filename.
        Multiple files can be processed via either an association file or wrapped by a
        ModelContainer.

    Returns
    -------
    result : ~roman_datamodels.datamodels.MosaicModel
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

    # TODO: provide 'drizpars' file (then remove _set_spec_defaults?)
    reference_file_types = []

    def process(self, input):
        if isinstance(input, datamodels.DataModel):
            input_models = ModelContainer([input])
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
            if hasattr(input_models, "asn_table") and len(
                input_models.asn_table
            ):
                output = input_models.asn_table["products"][0]["name"]
            elif hasattr(input_models[0], "meta"):
                output = input_models[0].meta.filename
        elif isinstance(input, ModelContainer):
            input_models = input
            output = f"{os.path.commonprefix([x.meta.filename for x in input_models])}.asdf"
            if len(output) == 0:
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

        #  Get drizzle parameters reference file, if there is one
        self.wht_type = self.weight_type
        if "drizpars" in self.reference_file_types:
            ref_filename = self.get_reference_file(input_models[0], "drizpars")
            self.log.info(f"Using drizpars reference file: {ref_filename}")
            kwargs = self.get_drizpars(ref_filename, input_models)
        else:
            # no drizpars reference file found
            self.log.info("No drizpars reference file found.")
            kwargs = self._set_spec_defaults()

        kwargs["allowed_memory"] = self.allowed_memory

        # Issue a warning about the use of exptime weighting
        if self.wht_type == "exptime":
            self.log.warning(
                "Use of EXPTIME weighting will result in incorrect"
            )
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
        model.meta.resample["pixel_scale_ratio"] = (
            self.pixel_scale / np.sqrt(model.meta.photometry.pixelarea_arcsecsq)
            if self.pixel_scale
            else self.pixel_scale_ratio
        )
        model.meta.resample["pixfrac"] = kwargs["pixfrac"]
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
            raise ValueError(
                f"Both '{name}' values must be either None or not None."
            )

    @staticmethod
    def _load_custom_wcs(asdf_wcs_file, output_shape):
        if not asdf_wcs_file:
            return None

        with asdf.open(asdf_wcs_file) as af:
            wcs = deepcopy(af.tree["wcs"])

        if output_shape is not None or wcs is None:
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

    def get_drizpars(self, ref_filename, input_models):
        """
        Extract drizzle parameters from reference file.

        This method extracts parameters from the drizpars reference file and
        uses those to set defaults on the following ResampleStep configuration
        parameters:

        pixfrac = float(default=None)
        kernel = string(default=None)
        fillval = string(default=None)
        wht_type = option('ivm', 'exptime', None, default=None)

        Once the defaults are set from the reference file, if the user has
        used a resample.cfg file or run ResampleStep using command line args,
        then these will overwrite the defaults pulled from the reference file.
        """
        with datamodels.DrizParsModel(ref_filename) as drpt:
            drizpars_table = drpt.data

        num_groups = len(input_models.group_names)
        filtname = input_models[0].meta.instrument.filter
        row = None
        filter_match = False
        # look for row that applies to this set of input data models
        for n, filt, num in zip(
            range(0, len(drizpars_table)),
            drizpars_table["filter"],
            drizpars_table["numimages"],
        ):
            # only remember this row if no exact match has already been made for
            # the filter. This allows the wild-card row to be anywhere in the
            # table; since it may be placed at beginning or end of table.

            if str(filt) == "ANY" and not filter_match and num_groups >= num:
                row = n
            # always go for an exact match if present, though...
            if filtname == filt and num_groups >= num:
                row = n
                filter_match = True

        # With presence of wild-card rows, code should never trigger this logic
        if row is None:
            self.log.error(
                "No row found in %s matching input data.", ref_filename
            )
            raise ValueError

        # Define the keys to pull from drizpars reffile table.
        # All values should be None unless the user set them on the command
        # line or in the call to the step

        drizpars = dict(
            pixfrac=self.pixfrac,
            kernel=self.kernel,
            fillval=self.fillval,
            wht_type=self.weight_type,
        )

        # For parameters that are set in drizpars table but not set by the
        # user, use these.  Otherwise, use values set by user.
        reffile_drizpars = {k: v for k, v in drizpars.items() if v is None}
        user_drizpars = {k: v for k, v in drizpars.items() if v is not None}

        # read in values from that row for each parameter
        for k in reffile_drizpars:
            if k in drizpars_table.names:
                reffile_drizpars[k] = drizpars_table[k][row]

        # Convert the strings in the FITS binary table from np.bytes_ to str
        for k, v in reffile_drizpars.items():
            if isinstance(v, np.bytes_):
                reffile_drizpars[k] = v.decode("UTF-8")

        all_drizpars = {**reffile_drizpars, **user_drizpars}

        kwargs = dict(
            good_bits=GOOD_BITS,
            single=self.single,
            blendheaders=self.blendheaders,
        )

        kwargs.update(all_drizpars)

        for k, v in kwargs.items():
            self.log.debug(f"   {k}={v}")

        return kwargs

    def _set_spec_defaults(self):
        """NIRSpec currently has no default drizpars reference file, so default
        drizzle parameters are not set properly.  This method sets them.

        Remove this class method when a drizpars reffile is delivered.
        """
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

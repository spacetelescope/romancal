import copy
import logging

import numpy as np
from astropy import units as u
from drizzle import cdrizzle, util
from roman_datamodels import datamodels, maker_utils

from ..datamodels import ModelLibrary
from . import meta_blender, resample_utils, resampler

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["OutputTooLargeError", "ResampleData"]


class OutputTooLargeError(MemoryError):
    """Raised when the output is too large for in-memory instantiation"""


class ResampleData:
    """
    This is the controlling routine for the resampling process.

    Notes
    -----
    This routine performs the following operations:

      1. Extracts parameter settings from input model, such as pixfrac,
         weight type, exposure time (if relevant), and kernel, and merges
         them with any user-provided values.
      2. Creates output WCS based on input images and define mapping function
         between all input arrays and the output array. Alternatively, a custom,
         user-provided WCS object can be used instead.
      3. Updates output data model with output arrays from drizzle, including
         a record of metadata from all input models.
    """

    def __init__(
        self,
        input_models,
        output=None,
        single=False,
        pixfrac=1.0,
        kernel="square",
        fillval="INDEF",
        wht_type="ivm",
        good_bits="0",
        pscale_ratio=1.0,
        pscale=None,
        **kwargs,
    ):
        """
        Parameters
        ----------
        input_models : ~romancal.datamodels.ModelLibrary
            A `~romancal.datamodels.ModelLibrary` object containing the data
            to be processed.

        output : str
            filename for output

        kwargs : dict
            Other parameters.

            .. note::
                ``output_shape`` is in the ``x, y`` order.

            .. note::
                ``in_memory`` controls whether or not the resampled
                array from ``resample_many_to_many()``
                should be kept in memory or written out to disk and
                deleted from memory. Default value is `True` to keep
                all products in memory.
        """
        # outlier detection calls this with:
        # - input_models
        # - single=True
        # - **pars....
        #   - weight_type
        #   - pixfrac
        #   - kernel
        #   - fillval
        #   - nlow
        #   - nhigh
        #   - maskpt
        #   - grow
        #   - snr
        #   - scale
        #   - backg
        #   - kernel_size
        #   - save_intermediate_results
        #   - resample_data
        #   - good_bits
        #   - allowed_memory
        #   - in_memory
        #   - make_output_path
        #   - resample_suffix
        # resample calls this with:
        # - input_models
        # - output=output (is output_filename below)
        # - **kwargs...
        #   - everything in spec
        #   - overwrite good_bits, pixfrac, kernel, fillval, wht_type, pscale_ratio (from pixel_scale_ratio)
        #   - allowed_memory
        #   - output_shape? (default None)
        #   - output_wcs? (default '', becomes None)
        #   - crpix? (default None)
        #   - crval? (default None)
        #   - rotation
        #   - pscale (from pixel_scale, default None)
        #   - pscale_ratio (from pixel_scale_ratio again?, default 1.0)
        #   - in_memory
        if (input_models is None) or (len(input_models) == 0):
            raise ValueError(
                "No input has been provided. Input must be a non-empty ModelLibrary"
            )

        self.input_models = input_models
        self.output_filename = output
        self.pscale_ratio = pscale_ratio
        self.single = single
        self.pixfrac = pixfrac
        self.kernel = kernel
        self.fillval = fillval
        self.weight_type = wht_type
        self.good_bits = good_bits
        self.in_memory = kwargs.get("in_memory", True)

        log.info(f"Driz parameter kernel: {self.kernel}")
        log.info(f"Driz parameter pixfrac: {self.pixfrac}")
        log.info(f"Driz parameter fillval: {self.fillval}")
        log.info(f"Driz parameter weight_type: {self.weight_type}")

        output_wcs = kwargs.get("output_wcs", None)
        output_shape = kwargs.get("output_shape", None)
        crpix = kwargs.get("crpix", None)
        crval = kwargs.get("crval", None)
        rotation = kwargs.get("rotation", None)

        if pscale is not None:
            log.info(f"Output pixel scale: {pscale} arcsec.")
            pscale /= 3600.0
        else:
            log.info(f"Output pixel scale ratio: {pscale_ratio}")

        # build the output WCS object
        if output_wcs:
            # use the provided WCS object
            self.output_wcs = output_wcs
            if output_shape is not None:
                self.output_wcs.array_shape = output_shape[::-1]
        else:
            with self.input_models:
                models = list(self.input_models)
                # determine output WCS based on all inputs, including a reference WCS
                self.output_wcs = resample_utils.make_output_wcs(
                    models,
                    pscale_ratio=self.pscale_ratio,
                    pscale=pscale,
                    rotation=rotation,
                    shape=None if output_shape is None else output_shape[::-1],
                    crpix=crpix,
                    crval=crval,
                )
                for i, m in enumerate(models):
                    self.input_models.shelve(m, i, modify=False)

        log.debug(f"Output mosaic size: {self.output_wcs.array_shape}")

        # NOTE: should we enable memory allocation?

        # can_allocate, required_memory = datamodels.util.check_memory_allocation(
        #     self.output_wcs.array_shape,
        #     kwargs['allowed_memory'],
        #     datamodels.ImageModel
        # )
        # if not can_allocate:
        #     raise OutputTooLargeError(
        #         f'Combined ImageModel size {self.output_wcs.array_shape} '
        #         f'requires {bytes2human(required_memory)}. '
        #         f'Model cannot be instantiated.'
        #     )

        # NOTE: wait for William to fix bug in datamodels' init and then
        # use datamodels.ImageModel(shape=(nx, ny)) instead of mk_datamodel()

        # self.blank_output = maker_utils.mk_datamodel(
        #     datamodels.MosaicModel, shape=tuple(self.output_wcs.array_shape)
        # )

        # with self.input_models:
        #     models = list(self.input_models)

        #     # update meta.basic
        #     populate_mosaic_basic(self.blank_output, models)

        #     # update meta.cal_step
        #     self.blank_output.meta.cal_step = maker_utils.mk_l3_cal_step(
        #         **models[0].meta.cal_step.to_flat_dict()
        #     )

        #     # Update the output with all the component metas
        #     populate_mosaic_individual(self.blank_output, models)

        #     # update meta data and wcs
        #     l2_into_l3_meta(self.blank_output.meta, models[0].meta)
        #     self.blank_output.meta.wcs = self.output_wcs
        #     gwcs_into_l3(self.blank_output, self.output_wcs)
        #     self.blank_output.cal_logs = stnode.CalLogs()
        #     self.blank_output["individual_image_cal_logs"] = [
        #         model.cal_logs for model in models
        #     ]
        #     for i, m in enumerate(models):
        #         self.input_models.shelve(m, i, modify=False)

    def do_drizzle(self):
        """Pick the correct drizzling mode based on ``self.single``."""
        if self.single:
            return self.resample_many_to_many()
        else:
            return self.resample_many_to_one()

    def resample_many_to_many(self):
        """Resample many inputs to many outputs where outputs have a common frame.

        Coadd only different detectors of the same exposure (e.g. map SCA 1 and
        10 onto the same output image), as they image different areas of the
        sky.

        Used for outlier detection
        """
        # this requires:
        # -- computed --
        # - self.blank_output (FIXME to remove this)
        #
        # -- from args or computed --
        # - self.output_wcs
        #
        # -- from args --
        # - self.input_models
        # - self.pixfrac
        # - self.kernel
        # - self.fillval
        # - self.weight_type
        # - self.good_bits
        # - self.in_memory
        # so what does this use in the class
        output_list = []
        for group_id, indices in self.input_models.group_indices.items():
            output_model = maker_utils.mk_datamodel(
                datamodels.MosaicModel, shape=tuple(self.output_wcs.array_shape)
            )
            output_model.meta["resample"] = maker_utils.mk_resample()
            output_model.meta.wcs = copy.deepcopy(self.output_wcs)
            blender = meta_blender.MetaBlender(output_model)

            # copy over asn information
            if (asn_pool := self.input_models.asn.get("asn_pool", None)) is not None:
                output_model.meta.asn.pool_name = asn_pool
            if (
                asn_table_name := self.input_models.asn.get("table_name", None)
            ) is not None:
                output_model.meta.asn.table_name = asn_table_name

            with self.input_models:
                example_image = self.input_models.borrow(indices[0])

                # Determine output file type from input exposure filenames
                # Use this for defining the output filename
                indx = example_image.meta.filename.rfind(".")
                output_type = example_image.meta.filename[indx:]
                output_root = "_".join(
                    example_image.meta.filename.replace(output_type, "").split("_")[:-1]
                )
                output_model.meta.filename = f"{output_root}_outlier_i2d{output_type}"

                self.input_models.shelve(example_image, indices[0], modify=False)

                # Initialize the output with the wcs
                # FIXME don't use the model yet...
                driz = resampler.Resampler(
                    output_model.data.value,
                    output_model.weight,
                    output_model.context.astype(
                        "int32"
                    ),  # FIXME is this needed? it makes a copy
                    self.output_wcs,
                    pixfrac=self.pixfrac,
                    kernel=self.kernel,
                    fillval=self.fillval,
                    rollover_context=True,
                )

                log.info(f"{len(indices)} exposures to drizzle together")
                for index in indices:
                    img = self.input_models.borrow(index)
                    # TODO: should weight_type=None here?
                    inwht = resample_utils.build_driz_weight(
                        img, weight_type=self.weight_type, good_bits=self.good_bits
                    )

                    # apply sky subtraction
                    if (
                        hasattr(img.meta, "background")
                        and img.meta.background.subtracted is False
                        and img.meta.background.level is not None
                    ):
                        data = img.data - img.meta.background.level
                    else:
                        data = img.data

                    driz.add_image(
                        data.value,
                        img.meta.wcs,
                        inwht,
                    )
                    del data
                    blender.blend(img)
                    self.input_models.shelve(img, index)

                # cast context array to uint32
                # FIXME do I need to assign weight and data back?
                output_model.context = output_model.context.astype(
                    "uint32"
                )  # FIXME another copy

                blender.finalize()

                # copy over asn information
                if not self.in_memory:
                    # Write out model to disk, then return filename
                    output_name = output_model.meta.filename
                    output_model.save(output_name)
                    log.info(f"Exposure {output_name} saved to file")
                    output_list.append(output_name)
                else:
                    output_list.append(output_model)
                del blender
                del output_model

        return ModelLibrary(output_list)

    def resample_many_to_one(self):
        """Resample and coadd many inputs to a single output.
        Used for level 3 resampling
        """
        # this requires:
        # -- computed --
        # - self.blank_output (FIXME to remove this)
        #
        # -- from args or computed --
        # - self.output_wcs
        #
        # -- from args --
        # - self.output_filename
        # - self.input_models
        # - self.pixfrac
        # - self.kernel
        # - self.fillval
        # - self.weight_type
        # - self.good_bits
        #
        # - self.resample_variance_array (function call)
        #   - self.drizzle_arrays (function call)
        #   - (other stuff from above)
        # - self.resample_exposure_time (function call)
        #   - (same as resample_variance_array)
        # - self.update_exposure_times (function call)
        output_model = maker_utils.mk_datamodel(
            datamodels.MosaicModel, shape=tuple(self.output_wcs.array_shape)
        )
        output_model.meta.wcs = copy.deepcopy(self.output_wcs)
        blender = meta_blender.MetaBlender(output_model)

        # TODO should this also be in many_to_many?
        output_model.meta.filename = self.output_filename
        output_model.meta["resample"] = maker_utils.mk_resample()
        output_model.meta.resample.weight_type = self.weight_type
        output_model.meta.resample.pointings = len(self.input_models.group_names)

        # copy over asn information
        if (asn_pool := self.input_models.asn.get("asn_pool", None)) is not None:
            output_model.meta.asn.pool_name = asn_pool
        if (
            asn_table_name := self.input_models.asn.get("table_name", None)
        ) is not None:
            output_model.meta.asn.table_name = asn_table_name

        # Initialize the output with the wcs
        # FIXME don't use a model yet...
        driz = resampler.Resampler(
            output_model.data.value,
            output_model.weight,
            output_model.context.astype("int32"),  # FIXME this casts and copies
            self.output_wcs,
            pixfrac=self.pixfrac,
            kernel=self.kernel,
            fillval=self.fillval,
        )

        log.info("Resampling science data")
        with self.input_models:
            for i, img in enumerate(self.input_models):
                inwht = resample_utils.build_driz_weight(
                    img,
                    weight_type=self.weight_type,
                    good_bits=self.good_bits,
                )
                if (
                    hasattr(img.meta, "background")
                    and img.meta.background.subtracted is False
                    and img.meta.background.level is not None
                ):
                    data = img.data - img.meta.background.level
                else:
                    data = img.data

                driz.add_image(
                    data.value,
                    img.meta.wcs,
                    inwht,
                )
                blender.blend(img)
                del data, inwht
                self.input_models.shelve(img, i, modify=False)

        # record the actual filenames (the expname from the association)
        # for each file used to generate the output_model
        output_model.meta.resample["members"] = [
            m["expname"] for m in self.input_models.asn["products"][0]["members"]
        ]

        # FIXME reuse pixmap etc for the variance and exposure time resampling
        # Resample variances array in self.input_models to output_model
        self.resample_variance_array("var_rnoise", output_model)
        self.resample_variance_array("var_poisson", output_model)
        self.resample_variance_array("var_flat", output_model)

        # Make exposure time image
        exptime_tot = self.resample_exposure_time(output_model)

        # TODO: fix unit here
        output_model.err = u.Quantity(
            np.sqrt(
                np.nansum(
                    [
                        output_model.var_rnoise,
                        output_model.var_poisson,
                        output_model.var_flat,
                    ],
                    axis=0,
                )
            ),
            unit=output_model.err.unit,
        )

        self.update_exposure_times(output_model, exptime_tot)

        # TODO: fix RAD to expect a context image datatype of int32
        # TODO do I need to copy back data and weight
        output_model.context = output_model.context.astype(
            np.uint32
        )  # FIXME this cast copies
        blender.finalize()

        return ModelLibrary([output_model])

    def resample_variance_array(self, name, output_model):
        """Resample variance arrays from ``self.input_models`` to the ``output_model``.

        Resample the ``name`` variance array to the same name in ``output_model``,
        using a cumulative sum.

        This modifies ``output_model`` in-place.
        """
        output_wcs = self.output_wcs
        inverse_variance_sum = np.full_like(output_model.data.value, np.nan)

        log.info(f"Resampling {name}")
        with self.input_models:
            for i, model in enumerate(self.input_models):
                variance = getattr(model, name)
                if variance is None or variance.size == 0:
                    log.debug(
                        f"No data for '{name}' for model "
                        f"{repr(model.meta.filename)}. Skipping ..."
                    )
                    continue
                elif variance.shape != model.data.shape:
                    log.warning(
                        f"Data shape mismatch for '{name}' for model "
                        f"{repr(model.meta.filename)}. Skipping..."
                    )
                    continue

                # create a unit weight map for all the input pixels with science data
                inwht = resample_utils.build_driz_weight(
                    model, weight_type=None, good_bits=self.good_bits
                )

                resampled_variance = np.zeros_like(output_model.data)
                outwht = np.zeros_like(output_model.data)
                outcon = np.zeros_like(output_model.context, dtype="i4")

                xmin, xmax, ymin, ymax = resample_utils.resample_range(
                    variance.shape, model.meta.wcs.bounding_box
                )

                # resample the variance array (fill "unpopulated" pixels with NaNs)
                self.drizzle_arrays(
                    variance,
                    inwht,
                    model.meta.wcs,
                    output_wcs,
                    resampled_variance,
                    outwht,
                    outcon,
                    pixfrac=self.pixfrac,
                    kernel=self.kernel,
                    fillval=np.nan,
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax,
                )

                # Add the inverse of the resampled variance to a running sum.
                # Update only pixels (in the running sum) with valid new values:
                mask = resampled_variance > 0

                inverse_variance_sum[mask] = np.nansum(
                    [
                        inverse_variance_sum[mask],
                        np.reciprocal(resampled_variance[mask]),
                    ],
                    axis=0,
                )
                self.input_models.shelve(model, i, modify=False)

        # We now have a sum of the inverse resampled variances.  We need the
        # inverse of that to get back to units of variance.
        # TODO: fix unit here
        output_variance = u.Quantity(
            np.reciprocal(inverse_variance_sum), unit=u.MJy**2 / u.sr**2
        )

        setattr(output_model, name, output_variance)

    def resample_exposure_time(self, output_model):
        """Resample the exposure time from ``self.input_models`` to the
        ``output_model``.

        Create an exposure time image that is the drizzled sum of the input
        images.
        """
        output_wcs = self.output_wcs
        exptime_tot = np.zeros(output_model.data.shape, dtype="f4")

        log.info("Resampling exposure time")
        with self.input_models:
            for i, model in enumerate(self.input_models):
                exptime = np.full(
                    model.data.shape, model.meta.exposure.effective_exposure_time
                )

                # create a unit weight map for all the input pixels with science data
                inwht = resample_utils.build_driz_weight(
                    model, weight_type=None, good_bits=self.good_bits
                )

                resampled_exptime = np.zeros_like(output_model.data)
                outwht = np.zeros_like(output_model.data)
                outcon = np.zeros_like(output_model.context, dtype="i4")
                # drizzle wants an i4, but datamodels wants a u4.

                xmin, xmax, ymin, ymax = resample_utils.resample_range(
                    exptime.shape, model.meta.wcs.bounding_box
                )

                # resample the exptime array
                self.drizzle_arrays(
                    exptime * u.s,  # drizzle_arrays expects these to have units
                    inwht,
                    model.meta.wcs,
                    output_wcs,
                    resampled_exptime,
                    outwht,
                    outcon,
                    pixfrac=1,  # for exposure time images, always use pixfrac = 1
                    kernel=self.kernel,
                    fillval=0,
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax,
                )

                exptime_tot += resampled_exptime.value
                self.input_models.shelve(model, i, modify=False)

        return exptime_tot

    def update_exposure_times(self, output_model, exptime_tot):
        """Update exposure time metadata (in-place)."""
        m = exptime_tot > 0
        total_exposure_time = np.mean(exptime_tot[m]) if np.any(m) else 0
        max_exposure_time = np.max(exptime_tot)
        log.info(
            f"Mean, max exposure times: {total_exposure_time:.1f}, "
            f"{max_exposure_time:.1f}"
        )
        exposure_times = {"start": [], "end": []}
        with self.input_models:
            for group_id, indices in self.input_models.group_indices.items():
                index = indices[0]
                model = self.input_models.borrow(index)
                exposure_times["start"].append(model.meta.exposure.start_time)
                exposure_times["end"].append(model.meta.exposure.end_time)
                self.input_models.shelve(model, index, modify=False)

        # Update some basic exposure time values based on output_model
        output_model.meta.basic.mean_exposure_time = total_exposure_time
        output_model.meta.basic.time_first_mjd = min(exposure_times["start"]).mjd
        output_model.meta.basic.time_last_mjd = max(exposure_times["end"]).mjd
        output_model.meta.basic.max_exposure_time = max_exposure_time
        output_model.meta.resample.product_exposure_time = max_exposure_time

    @staticmethod
    def drizzle_arrays(
        insci,
        inwht,
        input_wcs,
        output_wcs,
        outsci,
        outwht,
        outcon,
        uniqid=1,
        xmin=None,
        xmax=None,
        ymin=None,
        ymax=None,
        pixfrac=1.0,
        kernel="square",
        fillval="INDEF",
        wtscale=1.0,
    ):
        """
        Low level routine for performing 'drizzle' operation on one image.

        The interface is compatible with STScI code. All images are Python
        `ndarrays`, instead of filenames. File handling (input and output) is
        performed by the calling routine.

        Parameters
        ----------

        insci : 2d array
            A 2d `numpy` array containing the input image to be drizzled.

        inwht : 2d array
            A 2d `numpy` array containing the pixel by pixel weighting.
            Must have the same dimensions as insci. If none is supplied,
            the weighting is set to one.

        input_wcs : `gwcs.wcs.WCS` object
            The world coordinate system of the input image.

        output_wcs : `gwcs.wcs.WCS` object
            The world coordinate system of the output image.

        outsci : 2d array
            A 2d `numpy` array containing the output image produced by
            drizzling. On the first call it should be set to zero.
            Subsequent calls it will hold the intermediate results.  This
            is modified in-place.

        outwht : 2d array
            A 2d `numpy` array containing the output counts. On the first
            call it should be set to zero. On subsequent calls it will
            hold the intermediate results.  This is modified in-place.

        outcon : 2d or 3d array, optional
            A 2d or 3d `numpy` array holding a bitmap of which image was an input
            for each output pixel. Should be integer zero on first call.
            Subsequent calls hold intermediate results.  This is modified
            in-place.

        uniqid : int, optional
            The id number of the input image. Should be one the first time
            this function is called and incremented by one on each subsequent
            call.

        xmin : float, None, optional
            on the input image. Only pixels on the input image inside this
            rectangle will have their flux added to the output image. Xmin
            sets the minimum value of the x dimension. The x dimension is the
            dimension that varies quickest on the image. If the value is zero,
            no minimum will be set in the x dimension. All four parameters are
            zero based, counting starts at zero.

        xmax : float, None, optional
            Sets the maximum value of the x dimension on the bounding box
            of the input image. If the value is zero, no maximum will
            be set in the x dimension, the full x dimension of the output
            image is the bounding box.

        ymin : float, None, optional
            Sets the minimum value in the y dimension on the bounding box. The
            y dimension varies less rapidly than the x and represents the line
            index on the input image. If the value is zero, no minimum  will be
            set in the y dimension.

        ymax : float, None, optional
            Sets the maximum value in the y dimension. If the value is zero, no
            maximum will be set in the y dimension,  the full x dimension
            of the output image is the bounding box.

        pixfrac : float, optional
            The fraction of a pixel that the pixel flux is confined to. The
            default value of 1 has the pixel flux evenly spread across the image.
            A value of 0.5 confines it to half a pixel in the linear dimension,
            so the flux is confined to a quarter of the pixel area when the square
            kernel is used.

        kernel: str, optional
            The name of the kernel used to combine the input. The choice of
            kernel controls the distribution of flux over the kernel. The kernel
            names are: `'square'`, `'gaussian'`, `'point'`, `'tophat'`, `'turbo'`,
            `'lanczos2'`, and `'lanczos3'`. The `'square'` kernel is the default.

        fillval: str, optional
            The value a pixel is set to in the output if the input image does
            not overlap it. The default value of INDEF does not set a value.

        Returns
        -------
        : tuple
            A tuple with three values: a version string, the number of pixels
            on the input image that do not overlap the output image, and the
            number of complete lines on the input image that do not overlap the
            output input image.

        """

        # Insure that the fillval parameter gets properly interpreted for use with tdriz
        fillval = "INDEF" if util.is_blank(str(fillval)) else str(fillval)
        if insci.dtype > np.float32:
            insci = insci.astype(np.float32)

        # Add input weight image if it was not passed in
        if inwht is None:
            inwht = np.ones_like(insci)

        # Compute what plane of the context image this input would
        # correspond to:
        planeid = int((uniqid - 1) / 32)

        # Check if the context image has this many planes
        if outcon.ndim == 2:
            nplanes = 1
        elif outcon.ndim == 3:
            nplanes = outcon.shape[0]
        else:
            nplanes = 0

        if nplanes <= planeid:
            raise IndexError("Not enough planes in drizzle context image")

        # Alias context image to the requested plane if 3d
        if outcon.ndim == 3:
            outcon = outcon[planeid]

        # Compute the mapping between the input and output pixel coordinates
        # for use in drizzle.cdrizzle.tdriz
        pixmap = resample_utils.calc_gwcs_pixmap(input_wcs, output_wcs, insci.shape)

        log.debug(f"Pixmap shape: {pixmap[:,:,0].shape}")
        log.debug(f"Input Sci shape: {insci.shape}")
        log.debug(f"Output Sci shape: {outsci.shape}")

        log.info(f"Drizzling {insci.shape} --> {outsci.shape}")

        _vers, _nmiss, _nskip = cdrizzle.tdriz(
            insci.astype(np.float32).value,
            inwht,
            pixmap,
            outsci.value,
            outwht.value,
            outcon,
            uniqid=uniqid,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
            pixfrac=pixfrac,
            kernel=kernel,
            in_units="cps",
            expscale=1.0,
            wtscale=wtscale,
            fillstr=fillval,
        )

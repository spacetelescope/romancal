import logging

import numpy as np
from astropy import units as u
from drizzle import cdrizzle, util
from roman_datamodels import datamodels, maker_utils

from ..datamodels import ModelContainer
from . import gwcs_drizzle, resample_utils

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
        blendheaders=True,
        pixfrac=1.0,
        kernel="square",
        fillval="INDEF",
        wht_type="ivm",
        good_bits=0,
        pscale_ratio=1.0,
        pscale=None,
        **kwargs,
    ):
        """
        Parameters
        ----------
        input_models : list of objects
            list of data models, one for each input image

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
        if (
            (input_models is None)
            or (len(input_models) == 0)
            or (not any(input_models))
        ):
            raise ValueError(
                "No input has been provided. Input should be a list of datamodel(s) or "
                "a ModelContainer."
            )

        self.input_models = input_models
        self.output_filename = output
        self.pscale_ratio = pscale_ratio
        self.single = single
        self.blendheaders = blendheaders
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
            # determine output WCS based on all inputs, including a reference WCS
            self.output_wcs = resample_utils.make_output_wcs(
                self.input_models,
                pscale_ratio=self.pscale_ratio,
                pscale=pscale,
                rotation=rotation,
                shape=None if output_shape is None else output_shape[::-1],
                crpix=crpix,
                crval=crval,
            )

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

        self.blank_output = maker_utils.mk_datamodel(
            datamodels.MosaicModel, shape=tuple(self.output_wcs.array_shape)
        )

        # update meta data and wcs
        # note we have made this input_model_0 variable so that if
        # meta includes lazily-loaded objects, that we can successfully
        # copy them into the metadata.  Directly running input_models[0].meta
        # below can lead to input_models[0] going out of scope after
        # meta is loaded but before the dictionary is constructed,
        # which can lead to seek on closed file errors if
        # meta contains lazily loaded objects.
        input_model_0 = input_models[0]
        l2_into_l3_meta(self.blank_output.meta, input_model_0.meta)
        self.blank_output.meta.wcs = self.output_wcs
        gwcs_into_l3(self.blank_output, self.output_wcs)

        self.output_models = ModelContainer()

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
        for exposure in self.input_models.models_grouped:
            output_model = self.blank_output
            output_model.meta["resample"] = maker_utils.mk_resample()
            # Determine output file type from input exposure filenames
            # Use this for defining the output filename
            indx = exposure[0].meta.filename.rfind(".")
            output_type = exposure[0].meta.filename[indx:]
            output_root = "_".join(
                exposure[0].meta.filename.replace(output_type, "").split("_")[:-1]
            )
            output_model.meta.basic.filename = f"{output_root}_outlier_i2d{output_type}"

            # Initialize the output with the wcs
            driz = gwcs_drizzle.GWCSDrizzle(
                output_model,
                pixfrac=self.pixfrac,
                kernel=self.kernel,
                fillval=self.fillval,
            )

            log.info(f"{len(exposure)} exposures to drizzle together")
            output_list = []
            for img in exposure:
                img = datamodels.open(img)
                # TODO: should weight_type=None here?
                inwht = resample_utils.build_driz_weight(
                    img, weight_type=self.weight_type, good_bits=self.good_bits
                )

                # apply sky subtraction
                if not hasattr(img.meta, "background"):
                    self._create_background_attribute(img)
                blevel = img.meta.background.level
                if not img.meta.background.subtracted and blevel is not None:
                    data = img.data - blevel
                else:
                    data = img.data

                xmin, xmax, ymin, ymax = resample_utils.resample_range(
                    data.shape, img.meta.wcs.bounding_box
                )

                driz.add_image(
                    data,
                    img.meta.wcs,
                    inwht=inwht,
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax,
                )
                del data
                img.close()

            # cast context array to uint32
            output_model.context = output_model.context.astype("uint32")
            if not self.in_memory:
                # Write out model to disk, then return filename
                output_name = output_model.meta.basic.filename
                output_model.save(output_name)
                log.info(f"Exposure {output_name} saved to file")
                output_list.append(output_name)
            else:
                output_list.append(output_model.copy())

            self.output_models = ModelContainer(output_list, return_open=self.in_memory)
            output_model.data *= 0.0
            output_model.weight *= 0.0

        return self.output_models

    def resample_many_to_one(self):
        """Resample and coadd many inputs to a single output.
        Used for level 3 resampling
        """
        output_model = self.blank_output.copy()
        output_model.meta.basic.filename = self.output_filename
        output_model.meta["resample"] = maker_utils.mk_resample()
        output_model.meta.resample["members"] = []
        output_model.meta.resample.weight_type = self.weight_type
        output_model.meta.resample.pointings = len(self.input_models.models_grouped)

        if self.blendheaders:
            log.info("Skipping blendheaders for now.")

        # Initialize the output with the wcs
        driz = gwcs_drizzle.GWCSDrizzle(
            output_model,
            outwcs=self.output_wcs,
            pixfrac=self.pixfrac,
            kernel=self.kernel,
            fillval=self.fillval,
        )

        log.info("Resampling science data")
        for img in self.input_models:
            inwht = resample_utils.build_driz_weight(
                img, weight_type=self.weight_type, good_bits=self.good_bits
            )
            if not hasattr(img.meta, "background"):
                self._create_background_attribute(img)
            blevel = img.meta.background.level
            if not img.meta.background.subtracted and blevel is not None:
                data = img.data - blevel
            else:
                data = img.data

            xmin, xmax, ymin, ymax = resample_utils.resample_range(
                data.shape, img.meta.wcs.bounding_box
            )

            driz.add_image(
                data,
                img.meta.wcs,
                inwht=inwht,
                xmin=xmin,
                xmax=xmax,
                ymin=ymin,
                ymax=ymax,
            )
            del data, inwht
            output_model.meta.resample.members.append(str(img.meta.filename))

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
        output_model.context = output_model.context.astype(np.uint32)

        self.output_models.append(output_model)

        return self.output_models

    def _create_background_attribute(self, img):
        img.meta["background"] = {}
        img.meta.background["level"] = 0
        img.meta.background["subtracted"] = True

    def resample_variance_array(self, name, output_model):
        """Resample variance arrays from ``self.input_models`` to the ``output_model``.

        Resample the ``name`` variance array to the same name in ``output_model``,
        using a cumulative sum.

        This modifies ``output_model`` in-place.
        """
        output_wcs = self.output_wcs
        inverse_variance_sum = np.full_like(output_model.data.value, np.nan)

        log.info(f"Resampling {name}")
        for model in self.input_models:
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
            outcon = np.zeros_like(output_model.context)

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
        for model in self.input_models:
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
        for exposure in self.input_models.models_grouped:
            exposure_times["start"].append(exposure[0].meta.exposure.start_time)
            exposure_times["end"].append(exposure[0].meta.exposure.end_time)

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


def l2_into_l3_meta(l3_meta, l2_meta):
    """Update the level 3 meta with info from the level 2 meta

    Parameters
    ----------
    l3_meta : dict
        The meta to update. This is updated in-place

    l2_meta : stnode
        The Level 2-like meta to pull from

    Notes
    -----
    The list of meta that is pulled from the Level 2 meta into the Level 3 meta is as follows:
    basic.visit: observation.visit
    basic.segment: observation.segment
    basic.pass: observation.pass
    basic.program: observation.program
    basic.survey: obervation.survey
    basic.optical_element: optical_element
    basic.instrument: instrument.name
    basic.telescope: telescope
    program: program
    """
    l3_meta.basic.visit = l2_meta.observation.visit
    l3_meta.basic.segment = l2_meta.observation.segment
    l3_meta.basic["pass"] = l2_meta.observation["pass"]
    l3_meta.basic.program = l2_meta.observation.program
    l3_meta.basic.survey = l2_meta.observation.survey
    l3_meta.basic.optical_element = l2_meta.instrument.optical_element
    l3_meta.basic.instrument = l2_meta.instrument.name
    l3_meta.basic.telescope = l2_meta.telescope
    l3_meta.coordinates = l2_meta.coordinates
    l3_meta.program = l2_meta.program


def gwcs_into_l3(model, wcsinfo):
    """Update the Level 3 wcsinfo block from a GWCS object

    Parameters
    ----------
    model : `DataModel`
        The model whose meta is to be updated.

    wcsinfo : `GWCS`
        GWCS info to transfer into the `meta.wcsinfo` block

    Notes
    -----
    Some models/parameters in the GWCS object have explicit names, such as
    'crpix1'. However, some do not and hence have to be accessed explicitly
    by indexing. This is fragile and will be a source of issues.
    """
    l3_wcsinfo = model.meta.wcsinfo
    transform = wcsinfo.forward_transform

    # Basic WCS info
    l3_wcsinfo.projection = "TAN"
    l3_wcsinfo.rotation_matrix = transform["pc_rotation_matrix"].matrix.value.tolist()
    l3_wcsinfo.dec_ref = transform.lat_6.value
    l3_wcsinfo.ra_ref = transform.lon_6.value
    # l3_wcsinfo.x_ref = center of mosaic?
    # l3_wcsinfo.y_ref = center of mosaic?
    # l3_wcsinfo.s_region = from gwcs.bounding_box

    # Mosaic info
    # l3_wcsinfo.dec_center = from center of mosaic
    # l3_wcsinfo.pixel_scale = ???
    # l3_wcsinfo.pixel_scale_local = ???
    # l3_wcsinfo.pixel_shape = image shape
    # l3_wcsinfo.ra_center = from center of mosaic
    # l3_wcsinfo.ra_corn1 = from bounding box
    # l3_wcsinfo.ra_corn2 = from bounding box
    # l3_wcsinfo.ra_corn3 = from bounding box
    # l3_wcsinfo.ra_corn4 = from bounding box
    # l3_wcsinfo.dec_corn1 = from bounding box
    # l3_wcsinfo.dec_corn2 = from bounding box
    # l3_wcsinfo.dec_corn3 = from bounding box
    # l3_wcsinfo.dec_corn4 = from bounding box
    # l3_wcsinfo.orientat = ???
    # l3_wcsinfo.orientat_local = ???

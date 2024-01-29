"""Primary code for performing outlier detection on Roman observations."""

import logging
import warnings
from functools import partial

import numpy as np
from astropy.stats import sigma_clip
from astropy.units import Quantity
from drizzle.cdrizzle import tblot
from roman_datamodels import datamodels as rdm
from scipy import ndimage

from romancal.datamodels import ModelContainer
from romancal.lib import dqflags
from romancal.resample import resample
from romancal.resample.resample_utils import build_driz_weight, calc_gwcs_pixmap

from ..stpipe import RomanStep

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

DO_NOT_USE = dqflags.pixel["DO_NOT_USE"]
OUTLIER = dqflags.pixel["OUTLIER"]


__all__ = ["OutlierDetection", "flag_cr", "abs_deriv"]


class OutlierDetection:
    """Main class for performing outlier detection.

    This is the controlling routine for the outlier detection process.
    It loads and sets the various input data and parameters needed by
    the various functions and then controls the operation of this process
    through all the steps used for the detection.

    Notes
    -----
    This routine performs the following operations::

      1. Extracts parameter settings from input model and merges
         them with any user-provided values
      2. Resamples all input images into grouped observation mosaics.
      3. Creates a median image from all grouped observation mosaics.
      4. Blot median image to match each original input image.
      5. Perform statistical comparison between blotted image and original
         image to identify outliers.
      6. Updates input data model DQ arrays with mask of detected outliers.

    """

    default_suffix = "i2d"

    def __init__(self, input_models, **pars):
        """
        Initialize the class with input ModelContainers.

        Parameters
        ----------
        input_models : ~romancal.datamodels.container.ModelContainer
            A `~romancal.datamodels.container.ModelContainer` object containing the data
            to be processed.

        pars : dict, optional
            Optional user-specified parameters to modify how outlier_detection
            will operate.  Valid parameters include:
            - resample_suffix

        """
        self.input_models = input_models

        self.outlierpars = dict(pars)
        self.resample_suffix = f"_outlier_{self.default_suffix if pars.get('resample_suffix') is None else pars.get('resample_suffix')}.asdf"
        log.debug(f"Defined output product suffix as: {self.resample_suffix}")

        # Define how file names are created
        self.make_output_path = pars.get(
            "make_output_path", partial(RomanStep._make_output_path, None)
        )

    def do_detection(self):
        """Flag outlier pixels in DQ of input images."""  # self._convert_inputs()
        pars = self.outlierpars

        if pars["resample_data"]:
            # Start by creating resampled/mosaic images for
            # each group of exposures
            resamp = resample.ResampleData(
                self.input_models, single=True, blendheaders=False, **pars
            )
            drizzled_models = resamp.do_drizzle()

        else:
            # for non-dithered data, the resampled image is just the original image
            drizzled_models = self.input_models
            for i in range(len(self.input_models)):
                drizzled_models[i].weight = build_driz_weight(
                    self.input_models[i],
                    weight_type="ivm",
                    good_bits=pars["good_bits"],
                )

        # Initialize intermediate products used in the outlier detection
        median_model = rdm.open(drizzled_models[0]).copy()

        # Perform median combination on set of drizzled mosaics
        median_model.data = Quantity(
            self.create_median(drizzled_models), unit=median_model.data.unit
        )
        median_model_output_path = self.make_output_path(
            basepath=median_model.meta.basic.filename.replace(
                self.resample_suffix, ".asdf"
            ),
            suffix="median",
        )
        median_model.save(median_model_output_path)
        log.info(f"Saved model in {median_model_output_path}")

        if pars["resample_data"]:
            # Blot the median image back to recreate each input image specified
            # in the original input list/ASN/ModelContainer
            blot_models = self.blot_median(median_model)

        else:
            # Median image will serve as blot image
            blot_models = ModelContainer(return_open=False)
            for i in range(len(self.input_models)):
                blot_models.append(median_model)

        # Perform outlier detection using statistical comparisons between
        # each original input image and its blotted version of the median image
        self.detect_outliers(blot_models)

        # clean-up (just to be explicit about being finished with
        # these results)
        del median_model, blot_models

    def create_median(self, resampled_models):
        """Create a median image from the singly resampled images.

        NOTES
        -----
        This version is simplified from astrodrizzle's version in the
        following ways:
        - type of combination: fixed to 'median'
        - 'minmed' not implemented as an option
        """
        maskpt = self.outlierpars.get("maskpt", 0.7)

        log.info("Computing median")

        # Compute weight means without keeping DataModel for eacn input open
        # Start by insuring that the ModelContainer does NOT open and keep each datamodel
        ropen_orig = resampled_models._return_open
        resampled_models._return_open = False  # turn off auto-opening of models
        # keep track of resulting computation for each input resampled datamodel
        weight_thresholds = []
        # For each model, compute the bad-pixel threshold from the weight arrays
        for resampled in resampled_models:
            m = rdm.open(resampled)
            weight = m.weight
            # necessary in order to assure that mask gets applied correctly
            if hasattr(weight, "_mask"):
                del weight._mask
            mask_zero_weight = np.equal(weight, 0.0)
            mask_nans = np.isnan(weight)
            # Combine the masks
            weight_masked = np.ma.array(
                weight, mask=np.logical_or(mask_zero_weight, mask_nans)
            )
            # Sigma-clip the unmasked data
            weight_masked = sigma_clip(weight_masked, sigma=3, maxiters=5)
            mean_weight = np.mean(weight_masked)
            # Mask pixels where weight falls below maskpt percent
            weight_threshold = mean_weight * maskpt
            weight_thresholds.append(weight_threshold)
            # close and delete the model, just to explicitly try to keep the memory as clean as possible
            m.close()
            del m
        # Reset ModelContainer attribute to original value
        resampled_models._return_open = ropen_orig

        # Now, set up buffered access to all input models
        resampled_models.set_buffer(1.0)  # Set buffer at 1Mb
        resampled_sections = resampled_models.get_sections()
        median_image = np.empty(
            (resampled_models.imrows, resampled_models.imcols),
            resampled_models.imtype,
        )
        median_image[:] = np.nan  # initialize with NaNs

        for resampled_sci, resampled_weight, (row1, row2) in resampled_sections:
            # Create a mask for each input image, masking out areas where there is
            # no data or the data has very low weight
            badmasks = []
            for weight, weight_threshold in zip(resampled_weight, weight_thresholds):
                badmask = np.less(weight, weight_threshold)
                log.debug(
                    "Percentage of pixels with low weight: {}".format(
                        np.sum(badmask) / len(weight.flat) * 100
                    )
                )
                badmasks.append(badmask)

            # Fill resampled_sci array with nan's where mask values are True
            for f1, f2 in zip(resampled_sci, badmasks):
                for elem1, elem2 in zip(f1, f2):
                    elem1[elem2] = np.nan

            del badmasks

            # For a stack of images with "bad" data replaced with Nan
            # use np.nanmedian to compute the median.
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    action="ignore",
                    message="All-NaN slice encountered",
                    category=RuntimeWarning,
                )
                median_image[row1:row2] = np.nanmedian(resampled_sci, axis=0)
            del resampled_sci, resampled_weight

        return median_image

    def blot_median(self, median_model):
        """Blot resampled median image back to the detector images."""
        interp = self.outlierpars.get("interp", "linear")
        sinscl = self.outlierpars.get("sinscl", 1.0)

        # Initialize container for output blot images
        blot_models = []

        log.info("Blotting median")
        for model in self.input_models:
            blotted_median = model.copy()

            # clean out extra data not related to blot result
            blotted_median.err *= 0.0  # None
            blotted_median.dq *= 0  # None

            # apply blot to re-create model.data from median image
            blotted_median.data = Quantity(
                gwcs_blot(median_model, model, interp=interp, sinscl=sinscl),
                unit=blotted_median.data.unit,
            )

            model_path = self.make_output_path(
                basepath=model.meta.filename, suffix="blot"
            )
            blotted_median.save(model_path)
            log.info(f"Saved model in {model_path}")

            # Append model name to the ModelContainer so it is not passed in memory
            blot_models.append(model_path)

        return ModelContainer(blot_models, return_open=False)

    def detect_outliers(self, blot_models):
        """Flag DQ array for cosmic rays in input images.

        The science frame in each ImageModel in input_models is compared to
        the corresponding blotted median image in blot_models.  The result is
        an updated DQ array in each ImageModel in input_models.

        Parameters
        ----------
        input_models: JWST ModelContainer object
            data model container holding science ImageModels, modified in place

        blot_models : JWST ModelContainer object
            data model container holding ImageModels of the median output frame
            blotted back to the wcs and frame of the ImageModels in
            input_models

        Returns
        -------
        None
            The dq array in each input model is modified in place

        """
        log.info("Flagging outliers")
        for i, (image, blot) in enumerate(zip(self.input_models, blot_models)):
            blot = rdm.open(blot)
            flag_cr(image, blot, **self.outlierpars)
            self.input_models[i] = image


def flag_cr(
    sci_image,
    blot_image,
    snr="5.0 4.0",
    scale="1.2 0.7",
    backg=0,
    resample_data=True,
    **kwargs,
):
    """Masks outliers in science image by updating DQ in-place

    Mask blemishes in dithered data by comparing a science image
    with a model image and the derivative of the model image.

    Parameters
    ----------
    sci_image : ~romancal.DataModel.ImageModel
        the science data

    blot_image : ~romancal.DataModel.ImageModel
        the blotted median image of the dithered science frames

    snr : str
        Signal-to-noise ratio

    scale : str
        scaling factor applied to the derivative

    backg : float
        Background value (scalar) to subtract

    resample_data : bool
        Boolean to indicate whether blot_image is created from resampled,
        dithered data or not
    """
    snr1, snr2 = (float(val) for val in snr.split())
    scale1, scale2 = (float(val) for val in scale.split())

    # Get background level of science data if it has not been subtracted, so it
    # can be added into the level of the blotted data, which has been
    # background-subtracted
    if (
        hasattr(sci_image.meta, "background")
        and sci_image.meta.background.subtracted is False
        and sci_image.meta.background.level is not None
    ):
        subtracted_background = sci_image.meta.background.level
        log.debug(f"Adding background level {subtracted_background} to blotted image")
    else:
        # No subtracted background.  Allow user-set value, which defaults to 0
        subtracted_background = backg

    sci_data = sci_image.data
    blot_data = blot_image.data
    blot_deriv = abs_deriv(blot_data.value)
    err_data = np.nan_to_num(sci_image.err)

    # create the outlier mask
    if resample_data:  # dithered outlier detection
        blot_data += subtracted_background
        diff_noise = np.abs(sci_data - blot_data)

        # Create a boolean mask based on a scaled version of
        # the derivative image (dealing with interpolating issues?)
        # and the standard n*sigma above the noise
        threshold1 = scale1 * blot_deriv + snr1 * err_data.value
        mask1 = np.greater(diff_noise.value, threshold1)

        # Smooth the boolean mask with a 3x3 boxcar kernel
        kernel = np.ones((3, 3), dtype=int)
        mask1_smoothed = ndimage.convolve(mask1, kernel, mode="nearest")

        # Create a 2nd boolean mask based on the 2nd set of
        # scale and threshold values
        threshold2 = scale2 * blot_deriv + snr2 * err_data.value
        mask2 = np.greater(diff_noise.value, threshold2)

        # Final boolean mask
        cr_mask = mask1_smoothed & mask2

    else:  # stack outlier detection
        diff_noise = np.abs(sci_data - blot_data)

        # straightforward detection of outliers for non-dithered data since
        # err_data includes all noise sources (photon, read, and flat for baseline)
        cr_mask = np.greater(diff_noise.value, snr1 * err_data.value)

    # Count existing DO_NOT_USE pixels
    count_existing = np.count_nonzero(sci_image.dq & DO_NOT_USE)

    # Update the DQ array in the input image.
    sci_image.dq = np.bitwise_or(sci_image.dq, cr_mask * (DO_NOT_USE | OUTLIER))

    # Report number (and percent) of new DO_NOT_USE pixels found
    count_outlier = np.count_nonzero(sci_image.dq & DO_NOT_USE)
    count_added = count_outlier - count_existing
    percent_cr = count_added / (sci_image.shape[0] * sci_image.shape[1]) * 100
    log.info(f"New pixels flagged as outliers: {count_added} ({percent_cr:.2f}%)")


def abs_deriv(array):
    """Take the absolute derivate of a numpy array."""
    tmp = np.zeros(array.shape, dtype=np.float64)
    out = np.zeros(array.shape, dtype=np.float64)

    tmp[1:, :] = array[:-1, :]
    tmp, out = _absolute_subtract(array, tmp, out)
    tmp[:-1, :] = array[1:, :]
    tmp, out = _absolute_subtract(array, tmp, out)

    tmp[:, 1:] = array[:, :-1]
    tmp, out = _absolute_subtract(array, tmp, out)
    tmp[:, :-1] = array[:, 1:]
    tmp, out = _absolute_subtract(array, tmp, out)

    return out


def _absolute_subtract(array, tmp, out):
    tmp = np.abs(array - tmp)
    out = np.maximum(tmp, out)
    tmp = tmp * 0.0
    return tmp, out


def gwcs_blot(median_model, blot_img, interp="poly5", sinscl=1.0):
    """
    Resample the output/resampled image to recreate an input image based on
    the input image's world coordinate system

    Parameters
    ----------
    median_model : `~roman_datamodels.datamodels.MosaicModel`

    blot_img : datamodel
        Datamodel containing header and WCS to define the 'blotted' image

    interp : str, optional
        The type of interpolation used in the resampling. The
        possible values are "nearest" (nearest neighbor interpolation),
        "linear" (bilinear interpolation), "poly3" (cubic polynomial
        interpolation), "poly5" (quintic polynomial interpolation),
        "sinc" (sinc interpolation), "lan3" (3rd order Lanczos
        interpolation), and "lan5" (5th order Lanczos interpolation).

    sincscl : float, optional
        The scaling factor for sinc interpolation.
    """
    blot_wcs = blot_img.meta.wcs

    # Compute the mapping between the input and output pixel coordinates
    pixmap = calc_gwcs_pixmap(blot_wcs, median_model.meta.wcs, blot_img.data.shape)
    log.debug(f"Pixmap shape: {pixmap[:, :, 0].shape}")
    log.debug(f"Sci shape: {blot_img.data.shape}")

    pix_ratio = 1
    log.info(f"Blotting {blot_img.data.shape} <-- {median_model.data.shape}")

    outsci = np.zeros(blot_img.shape, dtype=np.float32)

    # Currently tblot cannot handle nans in the pixmap, so we need to give some
    # other value.  -1 is not optimal and may have side effects.  But this is
    # what we've been doing up until now, so more investigation is needed
    # before a change is made.  Preferably, fix tblot in drizzle.
    pixmap[np.isnan(pixmap)] = -1
    tblot(
        median_model.data,
        pixmap,
        outsci,
        scale=pix_ratio,
        kscale=1.0,
        interp=interp,
        exptime=1.0,
        misval=0.0,
        sinscl=sinscl,
    )

    return outsci

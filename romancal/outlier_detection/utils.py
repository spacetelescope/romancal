import copy
import logging
from functools import partial

import numpy as np
from roman_datamodels.dqflags import pixel
from stcal.outlier_detection.median import MedianComputer
from stcal.outlier_detection.utils import (
    compute_weight_threshold,
    flag_crs,
    flag_resampled_crs,
    gwcs_blot,
)

from romancal.resample.resample import ResampleData
from romancal.resample.resample_utils import build_driz_weight

from . import _fileio

__all__ = ["detect_outliers"]


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def _median_with_resampling(
    input_models,
    resamp,
    maskpt,
    save_intermediate_results,
    make_output_path,
    buffer_size=None,
):
    """
    Compute median of resampled data from models in a library.

    Parameters
    ----------
    input_models : ModelLibrary
        The input datamodels.

    resamp : resample.resample.ResampleData object
        The controlling object for the resampling process.

    maskpt : float
        The weight threshold for masking out low weight pixels.

    save_intermediate_results : bool
        if True, save the drizzled models and median model to fits.

    make_output_path : function
        The functools.partial instance to pass to save_median.

    buffer_size : int
        The size of chunk in bytes that will be read into memory when computing the median.
        This parameter has no effect if the input library has its on_disk attribute
        set to False. If None or 0 the buffer size will be set to the size of one
        resampled image.
    """
    in_memory = not input_models._on_disk
    indices_by_group = list(input_models.group_indices.values())
    nresultants = len(indices_by_group)
    example_model = None
    median_wcs = resamp.output_wcs

    with input_models:
        for i, indices in enumerate(indices_by_group):
            drizzled_model = resamp.resample_group(indices)

            if save_intermediate_results:
                # write the drizzled model to file
                _fileio.save_drizzled(drizzled_model, make_output_path)

            if i == 0:
                input_shape = (nresultants, *drizzled_model.data.shape)
                dtype = drizzled_model.data.dtype
                computer = MedianComputer(input_shape, in_memory, buffer_size, dtype)
                example_model = drizzled_model

            weight_threshold = compute_weight_threshold(drizzled_model.weight, maskpt)
            drizzled_model.data[drizzled_model.weight < weight_threshold] = np.nan
            computer.append(drizzled_model.data, i)
            del drizzled_model

    # Perform median combination on set of drizzled mosaics
    median_data = computer.evaluate()

    if save_intermediate_results:
        # drizzled model already contains asn_id
        _fileio.save_median(
            example_model,
            median_data,
            median_wcs,
            partial(make_output_path, asn_id=None),
        )

    return median_data, median_wcs


def _median_without_resampling(
    input_models,
    maskpt,
    weight_type,
    good_bits,
    save_intermediate_results,
    make_output_path,
    buffer_size=None,
):
    """
    Compute median of data from models in a library.

    Parameters
    ----------
    input_models : ModelLibrary
        The input datamodels.

    maskpt : float
        The weight threshold for masking out low weight pixels.

    weight_type : str
        The type of weighting to use when combining images. Options are:
        'ivm' (inverse variance) or 'exptime' (exposure time).

    good_bits : int
        The bit values that are considered good when determining the
        data quality of the input.

    save_intermediate_results : bool
        if True, save the models and median model to fits.

    make_output_path : function
        The functools.partial instance to pass to save_median.

    buffer_size : int
        The size of chunk in bytes that will be read into memory when computing the median.
        This parameter has no effect if the input library has its on_disk attribute
        set to False. If None or 0 the buffer size will be set to the size of one
        input image.

    """
    in_memory = not input_models._on_disk
    nresultants = len(input_models)
    example_model = None

    with input_models:
        for i in range(len(input_models)):
            model = input_models.borrow(i)
            wht = build_driz_weight(
                model,
                # FIXME this was hard-coded to "ivm"
                weight_type="ivm",
                good_bits=good_bits,
            )

            if save_intermediate_results:
                # write the model to file
                _fileio.save_drizzled(model, make_output_path)

            if i == 0:
                input_shape = (nresultants, *model.data.shape)
                dtype = model.data.dtype
                computer = MedianComputer(input_shape, in_memory, buffer_size, dtype)
                example_model = model
                median_wcs = copy.deepcopy(model.meta.wcs)

            weight_threshold = compute_weight_threshold(wht, maskpt)

            data_copy = model.data.copy()
            data_copy[wht < weight_threshold] = np.nan
            computer.append(data_copy, i)

            input_models.shelve(model, i, modify=True)
            del data_copy, model

    # Perform median combination on set of drizzled mosaics
    median_data = computer.evaluate()

    if save_intermediate_results:
        _fileio.save_median(example_model, median_data, median_wcs, make_output_path)

    return median_data, median_wcs


def _flag_resampled_model_crs(
    image,
    median_data,
    median_wcs,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
    save_intermediate_results,
    make_output_path,
):
    blot = gwcs_blot(median_data, median_wcs, image.data.shape, image.meta.wcs, 1.0)

    # Get background level of science data if it has not been subtracted, so it
    # can be added into the level of the blotted data, which has been
    # background-subtracted
    if (
        hasattr(image.meta, "background")
        and image.meta.background.subtracted is False
        and image.meta.background.level is not None
    ):
        backg = image.meta.background.level
        log.debug(
            f"Adding background level {image.meta.background.level} to blotted image"
        )

    cr_mask = flag_resampled_crs(
        image.data, image.err, blot, snr1, snr2, scale1, scale2, backg
    )

    # update the dq flags in-place
    image.dq |= cr_mask * np.uint32(pixel.DO_NOT_USE | pixel.OUTLIER)
    log.info(f"{np.count_nonzero(cr_mask)} pixels marked as outliers")


def _flag_model_crs(image, median_data, snr):
    cr_mask = flag_crs(image.data, image.err, median_data, snr)

    # Update dq array in-place
    image.dq |= cr_mask * np.uint32(pixel.DO_NOT_USE | pixel.OUTLIER)

    log.info(f"{np.count_nonzero(cr_mask)} pixels marked as outliers")


def detect_outliers(
    library,
    weight_type,
    pixfrac,
    kernel,
    fillval,
    maskpt,
    snr1,
    snr2,
    scale1,
    scale2,
    backg,
    save_intermediate_results,
    resample_data,
    good_bits,
    in_memory,
    make_output_path,
):
    # setup ResampleData
    # call
    if resample_data:
        resamp = ResampleData(
            library,
            None,
            pixfrac,
            kernel,
            fillval,
            # FIXME prior code provided weight_type when only wht_type is understood
            # both default to 'ivm' but tests that set this to something else did
            # not change the resampling weight type. For now, disabling it to match
            # the previous behavior.
            "ivm",
            good_bits,
            False,
            False,
            False,
            "",
        )
        median_data, median_wcs = _median_with_resampling(
            library,
            resamp,
            maskpt,
            save_intermediate_results=save_intermediate_results,
            make_output_path=make_output_path,
        )
    else:
        median_data, median_wcs = _median_without_resampling(
            library,
            maskpt,
            weight_type,
            good_bits,
            save_intermediate_results=save_intermediate_results,
            make_output_path=make_output_path,
        )

    # Perform outlier detection using statistical comparisons between
    # each original input image and its blotted version of the median image
    with library:
        for image in library:
            if resample_data:
                _flag_resampled_model_crs(
                    image,
                    median_data,
                    median_wcs,
                    snr1,
                    snr2,
                    scale1,
                    scale2,
                    backg,
                    save_intermediate_results,
                    make_output_path,
                )
            else:
                _flag_model_crs(image, median_data, snr1)

            # mark step as complete
            image.meta.cal_step["outlier_detection"] = "COMPLETE"

            library.shelve(image, modify=True)

    return library

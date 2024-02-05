"""
Create a source catalog for tweakreg
"""

import logging

import astropy.units as u
import numpy as np
from asdf import AsdfFile
from astropy.stats import SigmaClip
from astropy.table import Table
from photutils.background import (
    Background2D,
    MeanBackground,
    MedianBackground,
    ModeEstimatorBackground,
)
from photutils.detection import DAOStarFinder
from roman_datamodels import datamodels as rdm
from roman_datamodels import maker_utils

from romancal.lib import dqflags, psf
from romancal.stpipe import RomanStep

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["SourceDetectionStep"]


class SourceDetectionStep(RomanStep):
    """
    SourceDetectionStep: Detect point-like sources in image to create a catalog
    for alignment in tweakreg.
    """

    spec = """
        kernel_fwhm = float(default=None)  # DAOStarFinder:Size of Gaussian kernel,
        # in pixels.
        sharplo = float(default=0.)  # DAOStarFinder: Lower bound for sharpness.
        # Typical values of sharpness range from 0 (flat) to 1 (delta function).
        sharphi = float(default=1.0)  # DAOStarFinder: Upper bound for sharpness.
        roundlo = float(default=-1.0)  # DAOStarFinder: Lower bound for roundness.
        # A circular source will have a zero roundness. A source extended in x or
        # y will have a negative or positive roundness, respectively.
        roundhi = float(default=1.0)  # DAOStarFinder: Upper bound for roundness.
        peakmax = float(default=1000.0)  # Upper limit on brightest pixel in sources.
        max_sources = float(default=None)  # Max number of sources, choosing brightest.
        scalar_threshold = float(default=None) # Detection threshold, to
        # be used for entire image. Assumed to be in same units as data, and is
        # an absolute threshold over background.
        calc_threshold = boolean(default=True) # Calculate a single absolute
        # detection threshold from image based on background.
        snr_threshold = float(default=3.0)  # if calc_threshold_img,
        # the SNR for the threshold image.
        bkg_estimator = string(default='median')  # if calc_threshold_img,
        # choice of mean, median, or mode.
        bkg_boxsize = integer(default=9)  # if calc_threshold_img,
        # size of box in pixels for 2D background.
        bkg_sigma = float(default=2.0) # if calc_threshold_img,
        # n sigma for sigma clipping bkgrnd.
        bkg_filter_size = integer(default=3) # if calc_threshold_img,
        # size of Gaussian kernel for background.
        save_catalogs = boolean(default=False) # Save source catalog to file?
        # Will overwrite an existing catalog of the same name.
        output_cat_filetype = option('asdf', 'ecsv', default='asdf') # Used if
        #save_catalogs=True - file type of output catalog.
        fit_psf = boolean(default=False)  # fit source PSFs for accurate astrometry
    """

    def process(self, input):
        with rdm.open(input, lazy_load=False) as input_model:
            if input_model.meta.exposure.type != "WFI_IMAGE":
                # Check to see if attempt to find sources in non-Image data
                log.info("Skipping source detection for spectral exposure.")
                input_model.meta.cal_step.source_detection = "SKIPPED"
                return input_model

            # remove units from data in this step.
            # DAOStarFinder requires unitless input
            if hasattr(input_model.data, "unit"):
                self.data = input_model.data.value
            else:
                self.data = input_model.data

            # mask DO_NOT_USE pixels

            self.coverage_mask = (
                (dqflags.pixel["DO_NOT_USE"]) & input_model.dq
            ).astype(bool)

            filt = input_model.meta.instrument["optical_element"]

            # if a pre-determined threshold value for detection for the whole
            # image is provided, use this
            if self.scalar_threshold is not None:
                threshold = float(self.scalar_threshold)
                log.info(f"Using a detection threshold of {threshold}.")

            # otherwise, if specified, calculate a scalar threshold from the
            # image by calculating a 2D background image, using this to create
            # a 2d threshold image, and using the median of the 2d threshold
            # image as the scalar detection threshold for the whole image
            elif self.calc_threshold:
                log.info("Determining detection threshold from image.")
                bkg = self._calc_2D_background()
                threshold_img = bkg.background + self.snr_threshold * bkg.background_rms
                threshold = np.median(threshold_img)
                log.info(f"Calculated a detection threshold of {threshold} from image.")
            else:
                threshold = np.median(self.data)

            if self.kernel_fwhm is None:
                # estimate the FWHM of the PSF using rough estimate
                # from the diffraction limit:
                central_wavelength = psf.filter_central_wavelengths[filt] * u.um
                diffraction_limit = float(1.22 * central_wavelength / (2.4 * u.m))
                assume_pixel_scale = 0.11 * u.arcsec / u.pix
                expect_psf_pix = np.sqrt(
                    ((diffraction_limit * u.rad).to(u.arcsec) / assume_pixel_scale) ** 2
                    + 1 * u.pix**2
                ).to_value(u.pix)
                self.kernel_fwhm = expect_psf_pix

            log.info("Detecting sources with DAOStarFinder, using entire image array.")
            daofind = DAOStarFinder(
                fwhm=self.kernel_fwhm,
                threshold=threshold,
                sharplo=self.sharplo,
                sharphi=self.sharphi,
                roundlo=self.roundlo,
                roundhi=self.roundhi,
                brightest=self.max_sources,
                peakmax=self.peakmax,
            )

            if self.scalar_threshold is not None:
                # if absolute threshold is provided
                sources = daofind(self.data, mask=self.coverage_mask)

            elif self.calc_threshold is not None:
                # subtract background from data if calculating abs. threshold
                sources = daofind(self.data - bkg.background, mask=self.coverage_mask)
            else:
                sources = daofind(self.data, mask=self.coverage_mask)

            # reduce table to minimal number of columns, just source ID,
            # positions, and fluxes
            columns = ["id", "xcentroid", "ycentroid", "flux"]

            if sources:
                catalog = sources[columns]
                log.info(f"Found {len(catalog)} sources.")
            else:
                # if no sources were detected, return an empty table
                self.log.warning("No sources detected, returning empty catalog.")
                catalog = Table(
                    names=columns,
                    dtype=(int, np.float64, np.float64, np.float64),
                )

            if self.fit_psf and len(sources):
                # refine astrometry by fitting PSF models to the detected sources
                log.info("Constructing a gridded PSF model.")
                detector = input_model.meta.instrument["detector"].replace("WFI", "SCA")
                gridded_psf_model, _ = psf.create_gridded_psf_model(
                    path_prefix="tmp",
                    filt=filt,
                    detector=detector,
                    overwrite=True,
                    logging_level="ERROR",
                )

                log.info(
                    "Fitting a PSF model to sources for improved astrometric precision."
                )
                psf_photometry_table, photometry = psf.fit_psf_to_image_model(
                    image_model=input_model,
                    psf_model=gridded_psf_model,
                    x_init=catalog["xcentroid"],
                    y_init=catalog["ycentroid"],
                    exclude_out_of_bounds=True,
                )

            catalog_as_recarray = catalog.as_array()

            # create meta.source detection section in file
            # if save_catalogs is True, this will be updated with the
            # attribute 'tweakreg_catalog_name' to point to the location
            # of the catalog on disk. If save_catalogs is false, this section
            # will be updated to contain the catalog to pass to TweakReg

            # tweakreg_catalog_name will be saved to the final output file,
            # while tweakreg_catalog is intended to be deleted by TweakRegStep
            input_model.meta["source_detection"] = maker_utils.mk_source_detection()

            # if 'save_catalogs'= True, also save the output catalog to a file
            # (format specified by output_cat_filetype) and add an attribute
            # to the file that contains the path to this file
            if self.save_catalogs:
                cat_filename = input_model.meta.filename.replace(".asdf", "")
                cat_filename += f"_tweakreg_catalog.{self.output_cat_filetype}"
                log.info(f"Saving catalog to file: {cat_filename}.")

                if self.output_cat_filetype == "asdf":
                    tree = {"tweakreg_catalog": catalog_as_recarray}
                    ff = AsdfFile(tree)
                    ff.write_to(cat_filename)
                else:
                    catalog.write(cat_filename, format="ascii.ecsv", overwrite=True)

                input_model.meta.source_detection["tweakreg_catalog_name"] = (
                    cat_filename
                )
            else:
                if self.fit_psf:
                    # PSF photometry centroid results are stored in an astropy table
                    # in columns "x_fit" and "y_fit", which we'll rename for
                    # compatibility with tweakreg:
                    catalog = psf_photometry_table.copy()

                    # rename the PSF photometry table's columns to match the
                    # expectated columns in tweakreg:
                    for old_name, new_name in zip(
                        ["x_fit", "y_fit"], ["xcentroid", "ycentroid"]
                    ):
                        catalog.rename_column(old_name, new_name)

                    input_model.meta.source_detection["psf_catalog"] = catalog
                else:
                    # only attach catalog to file if its being passed to the next step
                    # and save_catalogs is false, since it is not in the schema
                    input_model.meta.source_detection["tweakreg_catalog"] = (
                        catalog_as_recarray
                    )
            input_model.meta.cal_step["source_detection"] = "COMPLETE"

        # just pass input model to next step - catalog is stored in meta
        return input_model

    def _calc_2D_background(self):
        """Calculates a 2D background image.

        Calculates the background value for the input image in boxes specified by
        self.bkg_box_size. A mean, median, or mode estimator may be used (set
        by `bkg_estimator`). The pixels in each box will be sigma clipped,
        using a sigma specified by `bkg_sigma`."""

        filter_size = (
            self.bkg_filter_size,
            self.bkg_filter_size,
        )  # square size specified
        box_size = np.asarray(self.bkg_boxsize).astype(int)  # must be integer

        if self.bkg_estimator == "median":
            bkg_estimator = MedianBackground()
        elif self.bkg_estimator == "mean":
            bkg_estimator = MeanBackground()
        elif self.bkg_estimator == "mode":
            bkg_estimator = ModeEstimatorBackground()
        else:
            raise ValueError("bkg_estimator must be one of 'mean', 'median', or 'mode'")

        sigma_clip = SigmaClip(self.bkg_sigma)

        try:
            bkg_2D = Background2D(
                self.data,
                box_size,
                filter_size=filter_size,
                coverage_mask=self.coverage_mask,
                sigma_clip=sigma_clip,
                bkg_estimator=bkg_estimator,
            )
        except ValueError:
            # use the entire unmasked array
            log.info(
                "Background could not be estimated in meshes. "
                "Using the entire unmasked array for background "
                f"estimation: bkg_boxsize={self.data.shape}."
            )

            bkg_2D = Background2D(
                self.data,
                self.data.shape,
                filter_size=filter_size,
                coverage_mask=self.coverage_mask,
                sigma_clip=sigma_clip,
                bkg_estimator=bkg_estimator,
                exclude_percentile=100.0,
            )

        return bkg_2D

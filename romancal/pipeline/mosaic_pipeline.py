#!/usr/bin/env python
import logging
import re
from os.path import basename, isfile

import asdf
import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy.modeling import models
from gwcs import WCS, coordinate_frames

import romancal.datamodels.filetype as filetype
from romancal.datamodels import ModelLibrary

# step imports
from romancal.flux import FluxStep
from romancal.outlier_detection import OutlierDetectionStep
from romancal.patch_match import patch_match
from romancal.resample import ResampleStep
from romancal.skymatch import SkyMatchStep
from romancal.source_catalog import SourceCatalogStep

from ..stpipe import RomanPipeline

__all__ = ["MosaicPipeline"]

# Define logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)


class MosaicPipeline(RomanPipeline):
    """
    MosaicPipeline: Apply all calibration steps to the roman data
    to produce level 3 products. Included steps are:
    ``skymatch``, ``outlier_detection`` and ``resample``.
    """

    class_alias = "roman_mos"
    spec = """
        save_results = boolean(default=False)
        on_disk = boolean(default=False)
    """

    # Define aliases to steps
    step_defs = {
        "flux": FluxStep,
        "skymatch": SkyMatchStep,
        "outlier_detection": OutlierDetectionStep,
        "resample": ResampleStep,
        "sourcecatalog": SourceCatalogStep,
    }

    # start the actual processing
    def process(self, input):
        """Process the Roman WFI data from Level 2 to Level 3"""

        log.info("Starting Roman mosaic level calibration pipeline ...")
        if isinstance(input, str):
            input_filename = basename(input)
        else:
            input_filename = None

        # open the input file
        file_type = filetype.check(input)
        if file_type == "asdf":
            log.info("The level three pipeline input needs to be an association")
            exit(0)
            return

        if file_type == "asn":
            input = ModelLibrary(input, on_disk=self.on_disk)
            self.flux.suffix = "flux"
            result = self.flux(input)
            self.skymatch.suffix = "skymatch"
            result = self.skymatch(result)
            self.outlier_detection.suffix = "outlier_detection"
            result = self.outlier_detection(result)
            #
            # check to see if the product name contains a skycell name & if true get the skycell record
            product_name = input.asn["products"][0]["name"]
            try:
                skycell_name = input.asn["target"]
            except IndexError:
                skycell_name = ""
            skycell_record = []

            # if this is a valid skycell name load the database and get the skycell record
            if re.match(r"r\d{3}\w{2}\d{2}x\d{2}y\d{2}", skycell_name):
                skycell_record = patch_match.PATCH_TABLE[
                    np.where(patch_match.PATCH_TABLE["name"][:] == skycell_name)[0][0]
                ]
                log.info("Skycell record %s:", skycell_record)

                if skycell_name in skycell_record["name"]:
                    # skycell name are in the form of r270dm90x99y99
                    # example of product name "r0099101001001001001_F158_visit_r270dm90x99y99"
                    skycell_root = re.findall("^[^%\n\r]*_", product_name)[0]
                    skycell_file_name = skycell_root + skycell_name + "_i2d.asdf"

                    # check to see if there exists a skycell on disk if not create it
                    if not isfile(skycell_file_name):
                        # extract the wcs info from the record for generate_tan_wcs
                        log.info(
                            "Creating skycell image at ra: %f  dec %f",
                            float(skycell_record["ra_center"]),
                            float(skycell_record["dec_center"]),
                        )
                        skycell_wcs = generate_tan_wcs(skycell_record)
                        # skycell_wcs.bounding_box = bounding_box

                        # For resample to use an external grid we need to pass it the skycell gwcs object
                        # Currently we cannot do that directly so create an asdf file to read the skycell gwcs object
                        wcs_tree = {"wcs": skycell_wcs}
                        wcs_file = asdf.AsdfFile(wcs_tree)
                        wcs_file.write_to("skycell_wcs.asdf")
                        self.resample.output_wcs = "skycell_wcs.asdf"
                        self.resample.output_shape = (
                            int(skycell_record["nx"]),
                            int(skycell_record["ny"]),
                        )
                        log.info(
                            "Resampling using %s  and data shape %s",
                            self.resample.output_wcs,
                            self.resample.output_shape,
                        )
                        wcs_file = asdf.open(self.resample.output_wcs)
                        self.suffix = "i2d"
                        result = self.resample(result)
                        self.output_file = input.asn["products"][0]["name"]
                        # force the SourceCatalogStep to save the results
                        self.sourcecatalog.save_results = True
                        result_catalog = self.sourcecatalog(result)
                    else:
                        log.info("resampling a mosaic file is not yet supported")
                        exit(0)

            else:
                self.resample.suffix = "i2d"
                self.output_file = input.asn["products"][0]["name"]
                result = self.resample(result)
                self.sourcecatalog.save_results = True
                result_catalog = self.sourcecatalog(result)  # noqa: F841
                self.suffix = "i2d"
                if input_filename:
                    result.meta.filename = self.output_file

        return result


def generate_tan_wcs(skycell_record):
    # extract the wcs info from the record for generate_tan_wcs
    # we need the scale, ra, dec, bounding_box

    scale = float(skycell_record["pixel_scale"])
    ra_center = float(skycell_record["ra_projection_center"])
    dec_center = float(skycell_record["dec_projection_center"])
    shiftx = float(skycell_record["x0_projection"])
    shifty = float(skycell_record["y0_projection"])
    bounding_box = (
        (-0.5, -0.5 + skycell_record["nx"]),
        (-0.5, -0.5 + skycell_record["ny"]),
    )

    # components of the model
    # shift = models.Shift(shiftx) & models.Shift(shifty)

    # select a scale for the skycell image, this will come from INS and may
    # be optimized for the different survey programs
    scale_x = scale
    scale_y = scale
    # set the pixelsscale to 0.1 arcsec/pixel
    pixelscale = models.Scale(scale_x / 3600.0) & models.Scale(scale_y / 3600.0)

    pixelshift = models.Shift(-1.0 * shiftx) & models.Shift(-1.0 * shifty)
    tangent_projection = models.Pix2Sky_TAN()
    celestial_rotation = models.RotateNative2Celestial(ra_center, dec_center, 180.0)
    det2sky = pixelshift | pixelscale | tangent_projection | celestial_rotation
    detector_frame = coordinate_frames.Frame2D(
        name="detector", axes_names=("x", "y"), unit=(u.pix, u.pix)
    )
    sky_frame = coordinate_frames.CelestialFrame(
        reference_frame=coordinates.ICRS(), name="icrs", unit=(u.deg, u.deg)
    )
    wcsobj = WCS([(detector_frame, det2sky), (sky_frame, None)])
    wcsobj.bounding_box = bounding_box

    return wcsobj

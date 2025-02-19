#!/usr/bin/env python
from __future__ import annotations

import logging
import re
from os.path import basename, isfile
from typing import TYPE_CHECKING

import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy.modeling import models
from gwcs import WCS, coordinate_frames
from stcal.alignment import util as wcs_util

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

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["MosaicPipeline"]

# Define logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)


class MosaicPipeline(RomanPipeline):
    """
    MosaicPipeline: Apply all calibration steps to the roman data
    to produce level 3 products. Included steps are:
    ``flux``, ``skymatch``, ``outlier_detection``, ``resample`` and ``source catalog``.
    """

    class_alias = "roman_mos"
    spec = """
        save_results = boolean(default=False)
        on_disk = boolean(default=False)
    """

    # Define aliases to steps
    step_defs: ClassVar = {
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
            raise TypeError("The level three pipeline input needs to be an association")

        if file_type == "asn":
            input = ModelLibrary(input, on_disk=self.on_disk)
            self.flux.suffix = "flux"
            result = self.flux.run(input)
            self.skymatch.suffix = "skymatch"
            result = self.skymatch.run(result)
            self.outlier_detection.suffix = "outlier_detection"
            result = self.outlier_detection.run(result)
            #
            # check to see if the product name contains a skycell name & if true get the skycell record
            product_name = input.asn["products"][0]["name"]
            try:
                skycell_name = input.asn["target"]
            except KeyError:
                skycell_name = ""
            skycell_record = []

            # if this is a valid skycell name get the skycell record
            if re.match(r"r\d{3}\w{2}\d{2}x\d{2}y\d{2}", skycell_name):
                # check to see if the skycell coords are in the asn header if
                # so read the string and convert to a dictionary to match the patch table
                if (
                    "skycell_wcs_info" in input.asn
                    and input.asn["skycell_wcs_info"] != "none"
                ):
                    skycell_record = input.asn["skycell_wcs_info"]
                else:
                    if patch_match.PATCH_TABLE is None:
                        patch_match.load_patch_table()
                    skycell_record = patch_match.PATCH_TABLE[
                        np.where(patch_match.PATCH_TABLE["name"][:] == skycell_name)[0][
                            0
                        ]
                    ]
                log.info("Skycell record %s:", skycell_record)

                if skycell_name in skycell_record["name"]:
                    # skycell name are in the form of r270dm90x99y99
                    # example of product name "r0099101001001001001_F158_visit_r270dm90x99y99"
                    skycell_file_name = product_name + "_coadd.asdf"

                    # check to see if there exists a skycell on disk if not create it
                    if not isfile(skycell_file_name):
                        # extract the wcs info from the record for skycell_to_wcs
                        log.info(
                            "Creating skycell image at ra: %f  dec %f",
                            float(skycell_record["ra_center"]),
                            float(skycell_record["dec_center"]),
                        )
                        skycell_wcs = skycell_to_wcs(skycell_record)
                        # skycell_wcs.bounding_box = bounding_box

                        # For resample to use an external grid we need to pass it the skycell gwcs object
                        self.resample.output_wcs = skycell_wcs
                        self.resample.output_shape = (
                            int(skycell_record["nx"]),
                            int(skycell_record["ny"]),
                        )
                        log.info(
                            "Resampling using %s  and data shape %s",
                            self.resample.output_wcs,
                            self.resample.output_shape,
                        )

                        self.suffix = "coadd"
                        self.output_file = input.asn["products"][0]["name"]
                        result = self.resample.run(result)
                        self.sourcecatalog.output_file = self.output_file
                        result_catalog = self.sourcecatalog.run(result)
                    else:
                        raise NotImplementedError(
                            "Overwriting an exisiting file or resampling a mosaic file is not yet supported"
                        )

            else:
                self.resample.suffix = "coadd"
                self.output_file = input.asn["products"][0]["name"]
                result = self.resample.run(result)
                self.sourcecatalog.output_file = self.output_file
                result_catalog = self.sourcecatalog.run(result)  # noqa: F841
                self.suffix = "coadd"
                if input_filename:
                    result.meta.filename = self.output_file

        return result


def skycell_to_wcs(skycell_record):
    """From a skycell record, generate a GWCS

    Parameters
    ----------
    skycell_record : dict
        A skycell record, or row, from the skycell patches table.

    Returns
    -------
    wcsobj : wcs.GWCS
        The GWCS object from the skycell record.
    """
    wcsinfo = dict()

    # The scale is given in arcseconds per pixel. Convert to degrees.
    wcsinfo["pixel_scale"] = float(skycell_record["pixel_scale"]) / 3600.0

    # Remaining components of the wcsinfo block
    wcsinfo["ra_ref"] = float(skycell_record["ra_projection_center"])
    wcsinfo["dec_ref"] = float(skycell_record["dec_projection_center"])
    wcsinfo["x_ref"] = float(skycell_record["x0_projection"])
    wcsinfo["y_ref"] = float(skycell_record["y0_projection"])
    wcsinfo["orientat"] = float(skycell_record["orientat_projection_center"])
    wcsinfo["rotation_matrix"] = None

    # Bounding box of the skycell. Note that the center of the pixels are at (0.5, 0.5)
    bounding_box = (
        (-0.5, -0.5 + skycell_record["nx"]),
        (-0.5, -0.5 + skycell_record["ny"]),
    )

    wcsobj = wcsinfo_to_wcs(wcsinfo, bounding_box=bounding_box)
    return wcsobj


def wcsinfo_to_wcs(wcsinfo, bounding_box=None, name="wcsinfo"):
    """Create a GWCS from the L3 wcsinfo meta

    Parameters
    ----------
    wcsinfo : dict or MosaicModel.meta.wcsinfo
        The L3 wcsinfo to create a GWCS from.

    bounding_box : None or 4-tuple
        The bounding box in detector/pixel space. Form of input is:
        (x_left, x_right, y_bottom, y_top)

    name : str
        Value of the `name` attribute of the GWCS object.

    Returns
    -------
    wcs : wcs.GWCS
        The GWCS object created.
    """
    pixelshift = models.Shift(-wcsinfo["x_ref"], name="crpix1") & models.Shift(
        -wcsinfo["y_ref"], name="crpix2"
    )
    pixelscale = models.Scale(wcsinfo["pixel_scale"], name="cdelt1") & models.Scale(
        wcsinfo["pixel_scale"], name="cdelt2"
    )
    tangent_projection = models.Pix2Sky_TAN()
    celestial_rotation = models.RotateNative2Celestial(
        wcsinfo["ra_ref"], wcsinfo["dec_ref"], 180.0
    )

    matrix = wcsinfo.get("rotation_matrix", None)
    if matrix:
        matrix = np.array(matrix)
    else:
        orientat = wcsinfo.get("orientat", 0.0)
        matrix = wcs_util.calc_rotation_matrix(
            np.deg2rad(orientat), v3i_yangle=0.0, vparity=1
        )
        matrix = np.reshape(matrix, (2, 2))
    rotation = models.AffineTransformation2D(matrix, name="pc_rotation_matrix")
    det2sky = (
        pixelshift | rotation | pixelscale | tangent_projection | celestial_rotation
    )

    detector_frame = coordinate_frames.Frame2D(
        name="detector", axes_names=("x", "y"), unit=(u.pix, u.pix)
    )
    sky_frame = coordinate_frames.CelestialFrame(
        reference_frame=coordinates.ICRS(), name="icrs", unit=(u.deg, u.deg)
    )
    wcsobj = WCS([(detector_frame, det2sky), (sky_frame, None)], name=name)

    if bounding_box:
        wcsobj.bounding_box = bounding_box

    return wcsobj

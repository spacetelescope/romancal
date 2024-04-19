#!/usr/bin/env python
import logging
from os.path import basename, isfile
import numpy as np
import re
import functools
from astropy.modeling import models
from astropy import coordinates
from astropy import units as u
from gwcs import WCS, coordinate_frames
#from gwcs.wcstools import wcs_from_fiducial
from asdf import AsdfFile
import asdf
import pdb

import romancal.datamodels.filetype as filetype
from romancal.datamodels import ModelContainer
import roman_datamodels as rdm
#from roman_datamodels import maker_utils

# step imports
from romancal.flux import FluxStep
from romancal.outlier_detection import OutlierDetectionStep
from romancal.resample import ResampleStep, resample_utils
from romancal.skymatch import SkyMatchStep
#from romancal.assign_wcs import AssignWcsStep
from romancal.assign_wcs import utils as awcs_utils

from ..stpipe import RomanPipeline

__all__ = ["HighLevelPipeline"]

# Define logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)


class HighLevelPipeline(RomanPipeline):
    """
    HighLevelPipeline: Apply all calibration steps to the roman data
    to produce level 3 products. Included steps are:
    ``skymatch``, ``outlier_detection`` and ``resample``.
    """

    class_alias = "roman_hlp"
    spec = """
        save_results = boolean(default=False)
    """

    # Define aliases to steps
    step_defs = {
        "flux": FluxStep,
        "skymatch": SkyMatchStep,
        "outlier_detection": OutlierDetectionStep,
        "resample": ResampleStep,
    }

    # start the actual processing
    def process(self, input):
        """Process the Roman WFI data from Level 2 to Level 3"""

        log.info("Starting Roman high level calibration pipeline ...")
        if isinstance(input, str):
            input_filename = basename(input)
        else:
            input_filename = None

        # open the input file
        file_type = filetype.check(input)
        if file_type == "asdf":
            log.info("The level three pipeline input needs to be an association")
            return

        if file_type == "asn":
            input = ModelContainer(input)
            self.flux.suffix = "flux"
            result = self.flux(input)
            self.skymatch.suffix = "skymatch"
            result = self.skymatch(result)
            self.outlier_detection.suffix = "outlier_detection"
            result = self.outlier_detection(result)
            #
            # This will be replaced with a call to get a record based on a skycell name
            skycell_record = np.array((463181, '2742.630.31.81', 269.81166407333336, 66.09965144, 1781.5, 1781.5, 355.9788   , 3564, 3564,   67715.5, -110484.5, 269.66579575,
                                       65.99686878, 269.64830329, 66.0952398 , 269.89132874, 66.10234972, 269.90791186, 66.00394719, 0.1, 274.28571429, 63., 0.))
            
            # check to see if the product name contains a skycell name
            product_name = input.asn_table["products"][0]["name"]
            #pdb.set_trace()

            if product_name.split("_")[3] in skycell_record[1]:
                # the skycell names are in flux, the latest format for the skycell part
                # these should be moved to a function/functions
                # seems to be r270dm90x99y99
                # example of product name "r0099101001001001001_F158_visit_2742_630_30_76_-05"
                skycell_root = re.findall('^[^%\n\r]*_', product_name)[0]
                skycell_ra = "r" + str(skycell_record[2])
                # construct the dec value with m for negative dec values
                if float(skycell_record[2]) < 0.:
                    skycell_dec = "dm" + str(np.abs(skycell_record[3]))
                else:
                    skycell_dec = "d" + str(skycell_record[3])
                skycell_file_name = skycell_root + skycell_ra + skycell_dec + "x" + str(skycell_record[4]) + \
                                                "y" + str(skycell_record[5])
                
                # check to see if there exists a skycell on disk if not create it
                if not isfile(skycell_file_name): 
                # extract the wcs info from the record for generate_tan_wcs
                    log.info("Creating skycell image at ra: %f  dec %f", float(skycell_record[2]),  float(skycell_record[3]) )
                    skycell_wcs = generate_tan_wcs( skycell_record )
                    #skycell_wcs.bounding_box = bounding_box
                    
                # For resample to use an external grid we need to pass it the skycell gwcs object
                # Currently we cannot do that directly so create an asdf file to read the skycell gwcs object
                wcs_tree = {"wcs": skycell_wcs}
                wcs_file = AsdfFile(wcs_tree)
                wcs_file.write_to("skycell_wcs.asdf")
                self.resample.output_wcs = "skycell_wcs.asdf"
                self.resample.output_shape = (int(skycell_record[7]), int(skycell_record[8]))
                log.info("Resampling using %s  and data shape %s", self.resample.output_wcs, self.resample.output_shape)
                wcs_file = asdf.open( self.resample.output_wcs)
                self.suffix = "i2d"
                result = self.resample(result)
                self.output_file = input.asn_table["products"][0]["name"]
            else:
                self.resample.suffix = "i2d"
                result = self.resample(result)
                self.suffix = "i2d"
                if input_filename:
                    result.meta.filename = self.output_file
                self.output_file = input.asn_table["products"][0]["name"]

        return result

def generate_tan_wcs(skycell_record, shiftx=0, shifty=0):
    # extract the wcs info from the record for generate_tan_wcs
    # we need the scale, ra, dec, bounding_box
    # Once we have an official skycell db the indexes below should
    # be replaced with the field name to be less fragile. 

    scale = float(skycell_record[19])
    ra_center = float(skycell_record[2])
    dec_center = float(skycell_record[3])
    bounding_box = (
        (-0.5, float(skycell_record[7]) + 0.5),
        (-0.5, float(skycell_record[8]) + 0.5),
            )
    shiftx = bounding_box[0][1]/2
    shifty = bounding_box[1][1]/2
    
    # components of the model
    shift = models.Shift(shiftx) & models.Shift(shifty)

    # select a scale for the skycell image, this will come from INS and may
    # be optimized for the different survey programs
    scale_x =  scale
    scale_y = scale
    # set the pixelsscale to 0.1 arcsec/pixel
    pixelscale = models.Scale(scale_x / 3600.) & models.Scale(scale_y / 3600.)

    pixelshift = models.Shift(-1.*shiftx) & models.Shift(-1.*shifty)
    tangent_projection = models.Pix2Sky_TAN()
    celestial_rotation = models.RotateNative2Celestial(ra_center, dec_center, 180.)
    det2sky = pixelshift | pixelscale | tangent_projection | celestial_rotation
    detector_frame = coordinate_frames.Frame2D(name="detector", axes_names=("x", "y"),unit=(u.pix, u.pix))
    sky_frame = coordinate_frames.CelestialFrame(reference_frame=coordinates.ICRS(), name='icrs', unit=(u.deg, u.deg))
    wcsobj = WCS([(detector_frame, det2sky), (sky_frame, None) ])
    wcsobj.bounding_box = bounding_box    
    
    return wcsobj

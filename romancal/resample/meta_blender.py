"""
Functions, all called during __init__ with:
- populate_mosaic_basic(output_model, input_models) -> blank_output, models
- (not a func):  mk_l3_cal_step(cal_step.to_flat_dict)  with models[0]
- populate_mosaic_individual(output_model, input_models) -> blank_output, models
- l2_int_l3_meta(l3_meta, l2_meta) -> black_output.meta, models[0].meta
- (not a func): blank_output.meta.wcs = output_wcs
- gwcs_into_l3(model, wcs) -> blank_output, output_wcs
- (not a func): blank_output.cal_logs = stnode.CalLogs()
- (not a func): blank_output["individual_image_cal_logs"] = list of input model.cal_logs

So splitting this up, some things use:
    - output wcs (this is pre-computed and in memory)
    - models[0] for:
        - cal_step
        - l2_into_l3
    - all models for:
        - populate_mosaic_basic
        - populate_mosaic_individual
        - ["individual_image_cal_logs"]

Much of this (all things that use all models) is wrong for many_to_many

(for the "output" model(s))
In many_to_many:
    - copy blank_output
    - overwrite meta.resample
    - copy over asn information (meta.asn.pool_name/table_name)
In many_to_one:
    - copy blank_output
    - overwrite meta.resample
    - set meta.resample.weight_type
    - set meta.resample.pointings
    - copy over asn information (meta.asn.pool_name/table_name)
    - set meta.resample["members"]
    - (in update_exposure_times)
        - meta.basic.mean_exposure_time
        - meta.basic_time_first_mjd (overwrites populate_mosaic_basic)
        - meta.basic.time_last_mjd (overwrites populate_mosaic_basic)
        - meta.basic.max_exposure_time
        - meta.resample.product_exposure_time
(in step)
    - set meta.cal_step["resample"] to "COMPLETE"
    - update meta.wcsinfo.s_region
    - update meta.resample.pixel_scale_ratio
    - update meta.resample.pixfrac
    - update meta.photometry.pixelarea_steradians
    - update meta.photometry.pixelarea_arcsecsq
    - update meta.resample["good_bits"]

It should be possible to combine several of the above into something that:
    - sets the wcs using output_wcs
        - meta.wcs = wcs
        - gwcs_into_l3
    - updates the initial meta (for the first model)
        - mk_l3_cal_step(cal_ste.to_flat_dict)
        - l2_into_l3_meta
    - updates meta for each subsequent model
    - accumulates meta for the final update
    - does final updates after the last model
"""

import logging

import numpy as np
from astropy.coordinates import SkyCoord
from roman_datamodels import maker_utils, stnode
from stcal.alignment.util import compute_scale

from ..assign_wcs import utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def gwcs_into_l3(model, wcs):
    """Update the Level 3 wcsinfo block from a GWCS object

    Parameters
    ----------
    model : `DataModel`
        The model whose meta is to be updated.

    wcs : `GWCS`
        GWCS info to transfer into the `meta.wcsinfo` block

    Notes
    -----
    Some models/parameters in the GWCS object have explicit names, such as
    'crpix1'. However, some do not and hence have to be accessed explicitly
    by indexing. This is fragile and will be a source of issues.
    """
    l3_wcsinfo = model.meta.wcsinfo
    transform = wcs.forward_transform

    l3_wcsinfo.projection = "TAN"
    l3_wcsinfo.pixel_shape = model.shape

    pixel_center = [(v - 1) / 2.0 for v in model.shape[::-1]]
    world_center = wcs(*pixel_center)
    l3_wcsinfo.ra_center = world_center[0]
    l3_wcsinfo.dec_center = world_center[1]
    l3_wcsinfo.pixel_scale_local = compute_scale(wcs, world_center)
    l3_wcsinfo.orientat_local = calc_pa(wcs, *world_center)

    try:
        footprint = utils.create_footprint(wcs, model.shape)
    except Exception as excp:
        log.warning("Could not determine footprint due to %s", excp)
    else:
        l3_wcsinfo.ra_corn1 = footprint[0][0]
        l3_wcsinfo.ra_corn2 = footprint[1][0]
        l3_wcsinfo.ra_corn3 = footprint[2][0]
        l3_wcsinfo.ra_corn4 = footprint[3][0]
        l3_wcsinfo.dec_corn1 = footprint[0][1]
        l3_wcsinfo.dec_corn2 = footprint[1][1]
        l3_wcsinfo.dec_corn3 = footprint[2][1]
        l3_wcsinfo.dec_corn4 = footprint[3][1]
        l3_wcsinfo.s_region = utils.create_s_region(footprint)

    try:
        l3_wcsinfo.x_ref = -transform["crpix1"].offset.value
        l3_wcsinfo.y_ref = -transform["crpix2"].offset.value
    except IndexError:
        log.warning(
            "WCS has no clear reference pixel defined by crpix1/crpix2. Assuming reference pixel is center."
        )
        l3_wcsinfo.x_ref = pixel_center[0]
        l3_wcsinfo.y_ref = pixel_center[1]
    world_ref = wcs(l3_wcsinfo.x_ref, l3_wcsinfo.y_ref)
    l3_wcsinfo.ra_ref = world_ref[0]
    l3_wcsinfo.dec_ref = world_ref[1]
    l3_wcsinfo.pixel_scale = compute_scale(wcs, world_ref)
    l3_wcsinfo.orientat = calc_pa(wcs, *world_ref)

    try:
        l3_wcsinfo.rotation_matrix = transform[
            "pc_rotation_matrix"
        ].matrix.value.tolist()
    except Exception:
        log.warning(
            "WCS has no clear rotation matrix defined by pc_rotation_matrix. Calculating one."
        )
        rotation_matrix = utils.calc_rotation_matrix(l3_wcsinfo.orientat, 0.0)
        l3_wcsinfo.rotation_matrix = utils.list_1d_to_2d(rotation_matrix, 2)


def calc_pa(wcs, ra, dec):
    """Calculate position angle at given ra,dec

    Parameters
    ----------
    wcs : GWCS
        The wcs in consideration.

    ra, dec : float, float
        The ra/dec in degrees.

    Returns
    -------
    position_angle : float
        The position angle in degrees.

    """
    delta_pix = [v for v in wcs.world_to_pixel(ra, dec)]
    delta_pix[1] += 1
    delta_coord = wcs.pixel_to_world(*delta_pix)
    coord = SkyCoord(ra, dec, frame="icrs", unit="deg")

    return coord.position_angle(delta_coord).degree


class MetaBlender:
    def __init__(self, output_model):
        self.output = output_model
        self.first = True
        # assume a model has an output wcs
        gwcs_into_l3(self.output, self.output.meta.wcs)

    def first_model(self, model):
        self.output.meta.resample = maker_utils.mk_resample()
        self.output.cal_logs = stnode.CalLogs()
        self.output["individual_image_cal_logs"] = [model.cal_logs]
        self.output.meta.cal_step = maker_utils.mk_l3_cal_step(
            **model.meta.cal_step.to_flat_dict()
        )
        self.output.meta.coordinates = model.meta.coordinates  # from l2_into_l3
        self.output.meta.program = model.meta.program  # from l2_into_l3
        # -- blend --
        self.output.meta.basic.time_first_mjd = model.meta.exposure.start_time.mjd
        self.output.meta.basic.time_last_mjd = model.meta.exposure.end_time_mjd
        self.output.meta.basic.visit = model.meta.observation.visit
        self.output.meta.basic.segment = model.meta.observation.segment
        self.output.meta.basic["pass"] = model.meta.observation[
            "pass"
        ]  # FIXME why getitem?
        self.output.meta.basic.program = model.meta.observation.program
        self.output.meta.basic.survey = (
            model.meta.observation.survey
        )  # use MULTIPLE for multiple
        # -- end blend --
        self.output.meta.basic.optical_element = model.meta.instrument.optical_element
        self.output.meta.basic.instrument = model.meta.instrument.name
        self.output.meta.basic.location_name = "TBD"
        self.output.meta.basic.product_type = "TBD"
        self.output.append_individual_image_meta(
            model.meta
        )  # was in populate_mosaic_individual

        # for time_mean_mjd we add the mean, so track all the values
        self._mean_mjd = [model.meta.exposure.mid_time.mjd]

    def blend(self, model):
        if self.first:
            self.first = False
            return self.first_model(model)
        # min first_mjd
        if model.meta.exposure.start_time_mjd < self.output.meta.basic.time_first_mjd:
            self.output.meta.basic.time_first_mjd = model.meta.exposure.start_time_mjd
        # max last_mjd
        if model.meta.exposure.end_time_mjd > self.output.meta.basic.time_last_mjd:
            self.output.meta.basic.time_last_mjd = model.meta.exposure.end_time_mjd
        # check visit matches
        if model.meta.observation.visit != self.output.meta.basic.visit:
            self.output.meta.basic.visit = -1
        # check segment matches
        if model.meta.observation.segment != self.output.meta.basic.segment:
            self.output.meta.basic.segment = -1
        # check pass matches
        if model.meta.observation["pass"] != self.output.meta.basic["pass"]:
            self.output.meta.basic["pass"] = -1
        # check program matches
        if model.meta.observation.program != self.output.meta.basic.program:
            self.output.meta.basic.program = -1
        # check survey matches (if not MULTIPLE)
        if model.meta.observation.survey != self.output.meta.basic.survey:
            self.output.meta.basic.survey = "MULTIPLE"
        # append individual image meta
        self.output.append_individual_image_meta(model.meta)
        # accumulate mid_time for mean
        self._mean_mjd.append(model.meta.exposure.mid_time.mjd)

    def finalize(self):
        # mean of self._mean_mjd into ...
        self.output.meta.basic.time_mean_mjd = np.mean(self._mean_mjd)

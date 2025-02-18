import logging

import numpy as np
from astropy.coordinates import SkyCoord
from roman_datamodels import datamodels, dqflags, maker_utils
from stcal.alignment.util import (
    compute_s_region_keyword,
    compute_scale,
)
from stcal.resample import Resample

from ..assign_wcs import utils
from .resample_utils import make_output_wcs

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["ResampleData"]


class ResampleData(Resample):
    dq_flag_name_map = dqflags.pixel

    def __init__(
        self,
        input_models,
        output_wcs,
        pixfrac,
        kernel,
        fillval,
        weight_type,
        good_bits,
        enable_ctx,
        enable_var,
        compute_err,
        wcs_kwargs=None,
    ):
        self.input_models = input_models

        if output_wcs is None:
            if wcs_kwargs is None:
                wcs_kwargs = {}
            with input_models:
                models = list(input_models)
                wcs = make_output_wcs(
                    models,
                    **wcs_kwargs,
                )
                for i, m in enumerate(models):
                    input_models.shelve(m, i, modify=False)
            # FIXME use sregions here, but this causes data differences
            # sregions = []
            # ref_wcs = None
            # ref_wcsinfo = None
            # with input_models:
            #     for model in input_models:
            #         if ref_wcs is None:
            #             ref_wcs = model.meta.wcs
            #             ref_wcsinfo = model.meta.wcsinfo
            #         sregions.append(model.meta.wcsinfo.s_region)
            #         input_models.shelve(model, modify=False)
            #     output_wcs = ref_wcs
            # TODO compute pixel_scale pixel_scale_ratio (if not defined)
            # wcs = wcs_from_sregions(
            #     sregions,
            #     ref_wcs=ref_wcs,
            #     ref_wcsinfo=ref_wcsinfo,
            #     **wcs_kwargs,
            # )
            output_wcs = {
                "wcs": wcs,
                # "pixel_scale": self.pixel_scale,
                # "pixel_scale_ratio": self.pixel_scale_ratio,
            }

        super().__init__(
            output_wcs,
            n_input_models=len(input_models),
            pixfrac=pixfrac,
            kernel=kernel,
            fillval=fillval,
            weight_type=weight_type,
            good_bits=good_bits,
            enable_ctx=enable_ctx,
            enable_var=enable_var,
            compute_err=compute_err,
        )

    def _input_model_to_dict(self, model):
        pixel_area = model.meta.photometry.pixel_area
        if pixel_area == -999999:
            pixel_area = None
        exposure_time = model.meta.exposure.exposure_time
        if exposure_time == -999999:
            exposure_time = 1
        if "background" in model.meta:
            level = model.meta.background.level
            subtracted = model.meta.background.subtracted
        else:
            level = 0
            subtracted = True
        return {
            "data": model.data,
            "dq": model.dq,
            "var_rnoise": model.var_rnoise,
            "var_poisson": model.var_poisson,
            "var_flat": model.var_flat,
            "err": model.err,
            "filename": model.meta.filename,
            "wcs": model.meta.wcs,
            "pixelarea_steradians": pixel_area,
            "group_id": model.meta.group_id,
            "measurement_time": None,  # falls back to exposure_time
            "exposure_time": exposure_time,
            "start_time": model.meta.exposure.start_time,
            "end_time": model.meta.exposure.end_time,
            "duration": None,  # unused
            "level": level,
            "subtracted": subtracted,
        }

    def _get_intensity_scale(self, model):
        # FIXME we lie about this to retain the old behavior
        return 1

    def add_model(self, model):
        model_dict = self._input_model_to_dict(model)
        # FIXME we lie about a few things to retain the old behavior
        model_dict["exposure_time"] = 1
        super().add_model(model_dict)
        # TODO blend metadata

    def finalize(self):
        # TODO finish blending
        super().finalize()

        # TODO update output model and return it
        # self.output_model is some "dict" with some unknown set of attributes

        output_model = maker_utils.mk_datamodel(
            datamodels.MosaicModel, n_images=0, shape=(0, 0)
        )

        output_model.meta.resample.good_bits = self.good_bits
        output_model.meta.resample.weight_type = self.weight_type
        output_model.meta.resample.pixfrac = self.output_model["pixfrac"]
        # output_model.meta.resample.product_exposure_time = ?

        # record the actual filenames (the expname from the association)
        # for each file used to generate the output_model
        # TODO this is incorrect for resample_group
        output_model.meta.resample["members"] = [
            m["expname"] for m in self.input_models.asn["products"][0]["members"]
        ]

        pixel_scale_ratio = self.output_model["pixel_scale_ratio"]
        if pixel_scale_ratio is not None:
            output_model.meta.resample.pixel_scale_ratio = pixel_scale_ratio

        output_model.meta.resample.pointings = len(self.input_models.group_names)
        # output_model.meta.resample.pointings = self.output_model["pointings"]

        output_model.meta.basic.location_name = self.input_models.asn.get(
            "target", "None"
        )

        # copy over asn information
        if (asn_pool := self.input_models.asn.get("asn_pool", None)) is not None:
            output_model.meta.asn.pool_name = asn_pool
        if (
            asn_table_name := self.input_models.asn.get("table_name", None)
        ) is not None:
            output_model.meta.asn.table_name = asn_table_name

        # every resampling will generate these
        output_model.data = self.output_model["data"]

        output_model.weight = self.output_model["wht"]

        # some things are conditional
        if self._enable_ctx:
            output_model.context = self.output_model["con"].astype(np.uint32)
        for arr_name in ("err", "var_rnoise", "var_poisson", "var_flat"):
            if arr_name in self.output_model:
                new_array = self.output_model[arr_name]
                if new_array is not None:
                    setattr(output_model, arr_name, new_array)

        # assign wcs to output model
        output_model.meta.wcs = self.output_wcs
        gwcs_into_l3(output_model, output_model.meta.wcs)
        return output_model

    def reset_arrays(self, n_input_models=None):
        super().reset_arrays(n_input_models)

        # FIXME even though we tell drizzle INDEF it sets
        # this to NaN if we don't provide an output array
        # so we hard code it here. This is bad since fillval
        # could be something else...
        self._driz._fillval = "INDEF"

        # TODO make model blender, could be done in __init__?

    def resample_group(self, indices):
        if self.is_finalized():
            self.reset_arrays(len(indices))

        with self.input_models:
            for index in indices:
                model = self.input_models.borrow(index)
                self.add_model(model)
                # TODO modify False?
                self.input_models.shelve(model, index)

        return self.finalize()


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

    # Fill out image-local information
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
        l3_wcsinfo.s_region = compute_s_region_keyword(footprint)

    # Fill out wcs-general information
    try:
        l3_wcsinfo.x_ref = -transform["crpix1"].offset.value
        l3_wcsinfo.y_ref = -transform["crpix2"].offset.value
    except IndexError:
        log.warning(
            "WCS has no clear reference pixel defined by crpix1/crpix2. Assuming reference pixel is center."
        )
        l3_wcsinfo.x_ref = pixel_center[0]
        l3_wcsinfo.y_ref = pixel_center[1]

    world_ref = wcs(l3_wcsinfo.x_ref, l3_wcsinfo.y_ref, with_bounding_box=False)
    l3_wcsinfo.ra_ref = world_ref[0]
    l3_wcsinfo.dec_ref = world_ref[1]

    try:
        cdelt1 = transform["cdelt1"].factor.value
        cdelt2 = transform["cdelt2"].factor.value
        l3_wcsinfo.pixel_scale = (cdelt1 + cdelt2) / 2.0
    except IndexError:
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
    delta_pix = [v for v in wcs.invert(ra, dec, with_bounding_box=False)]
    delta_pix[1] += 1
    delta_coord = SkyCoord(
        *wcs(*delta_pix, with_bounding_box=False), frame="icrs", unit="deg"
    )
    coord = SkyCoord(ra, dec, frame="icrs", unit="deg")

    return coord.position_angle(delta_coord).degree

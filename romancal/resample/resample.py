import logging

import numpy as np
from astropy.coordinates import SkyCoord
from roman_datamodels import datamodels, dqflags, maker_utils
from stcal.alignment.util import compute_s_region_keyword, compute_scale
from stcal.resample import Resample

from ..assign_wcs import utils

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
        output_filename,  # to allow meta.filename setting
    ):
        self.output_filename = output_filename

        self.input_models = input_models

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
        model_dict = {
            "data": model.data,
            "dq": model.dq,
            "filename": model.meta.filename,
            "wcs": model.meta.wcs,
            "pixelarea_steradians": pixel_area,
            "group_id": model.meta.group_id,
            "measurement_time": None,  # falls back to exposure_time
            "exposure_time": model.meta.exposure.exposure_time,
            "start_time": model.meta.exposure.start_time,
            "end_time": model.meta.exposure.end_time,
            "duration": None,  # unused
            "level": model.meta.background.level,
            "subtracted": model.meta.background.subtracted,
        }

        if self.enable_var:
            model_dict["var_rnoise"] = model.var_rnoise
            model_dict["var_poisson"] = model.var_poisson
            model_dict["var_flat"] = model.var_flat

        # TODO this is only needed when compute_err=driz_err?
        model_dict["err"] = model.err
        return model_dict

    def add_model(self, model):
        super().add_model(self._input_model_to_dict(model))
        # TODO blend metadata

    def finalize(self):
        # TODO finish blending
        super().finalize()

        # TODO update output model and return it
        # self.output_model is some "dict" with some unknown set of attributes

        output_model = maker_utils.mk_datamodel(
            datamodels.MosaicModel, n_images=0, shape=(0, 0)
        )
        output_model.meta.filename = self.output_filename

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

        output_model.meta.resample.pixel_scale_ratio = self.output_model[
            "pixel_scale_ratio"
        ]
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

        # some things are conditional
        if ctx := self.output_model.get("ctx"):
            output_model.context = ctx.astype(np.uint32)
        for arr_name in ("err", "var_rnoise", "var_poisson", "var_flat"):
            if arr_name in self.output_model:
                setattr(output_model, arr_name, self.output_model[arr_name])
        return output_model

    def reset_arrays(self, n_input_models=None):
        # TODO why is this overridden?
        super().reset_arrays(n_input_models)
        # make model blender, could be done in __init__?

    def resample_group(self, indices):
        # required by outlier detection ONLY
        # call reset_arrays (this is ALSO called by __init__ so not sure if needed)
        #   maybe only call it if finalized is True need to call with reset_output=True
        # ... add_model
        # finalize
        #
        # if we already resampled a group this instance will be "finalized" and
        # require resetting before we can add new models
        if self.finalized:
            self.reset_arrays(reset_output=True, n_output_models=len(indices))

        with self.input_models:
            for index in indices:
                model = self.input_models.borrow(index)
                self.add_model(model)
                # TODO modify False?
                self.input_models.shelve(model, index)

        return self.finalize()

    def do_drizzle(self):
        # old API used by step, jwst PR removes this and replaces it with
        # either a many_to_many or many_to_one call
        #
        # no need to call reset_arrays since it's called by __init__
        # just call add_model... then finalize
        # TODO resample_group with range(len(self.input_models))?
        with self.input_models:
            for index, model in enumerate(self.input_models):
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

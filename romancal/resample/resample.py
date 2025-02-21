import logging

import numpy as np
from roman_datamodels import datamodels, dqflags, maker_utils
from stcal.resample import Resample

from .exptime_resampler import ExptimeResampler
from .l3_wcs import assign_l3_wcs
from .meta_blender import MetaBlender
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
        compute_exptime,
        blend_meta,
        wcs_kwargs=None,
    ):
        # fillval indef doesn't define a starting value
        # since we're not resampling onto an existing array
        # we overwrite this to 0 (to match the old behavior)
        if fillval.lower() == "indef":
            fillval = 0
        self.input_models = input_models

        self._compute_exptime = compute_exptime
        self._blend_meta = blend_meta

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

            # FIXME use s_regions here, but this causes data differences
            # since s_regions are not using center=False
            # sregions = []
            # ref_wcs = None
            # ref_wcsinfo = None
            # shape = None
            # with input_models:
            #     for model in input_models:
            #         if ref_wcs is None:
            #             ref_wcs = model.meta.wcs
            #             ref_wcsinfo = model.meta.wcsinfo
            #             shape = model.data.shape
            #         sregions.append(model.meta.wcsinfo.s_region)
            #         input_models.shelve(model, modify=False)
            #     output_wcs = ref_wcs

            # if (pscale := wcs_kwargs.get("pixel_scale")) is not None:
            #     pscale /= 3600.0

            # pixel_scale_ratio = wcs_kwargs.get("pixel_scale_ratio", 1.0)

            # if pscale is None:
            #     pscale_in0 = compute_scale(
            #         ref_wcs,
            #         fiducial=np.array([ref_wcsinfo["ra_ref"], ref_wcsinfo["dec_ref"]])
            #     )
            #     pixel_scale = pscale_in0 * pixel_scale_ratio
            # else:
            #     pscale_in0 = np.rad2deg(
            #         math.sqrt(compute_mean_pixel_area(ref_wcs, shape=shape))
            #     )

            #     pixel_scale_ratio = pixel_scale / pscale_in0

            # wcs = wcs_from_sregions(
            #     sregions,
            #     ref_wcs=ref_wcs,
            #     ref_wcsinfo=ref_wcsinfo,
            #     pscale_ratio=pixel_scale_ratio,
            #     pscale=pixel_scale,
            #     shape=wcs_kwargs.get("shape"),
            #     rotation=wcs_kwargs.get("rotation"),
            #     crpix=wcs_kwargs.get("crpix"),
            #     crval=wcs_kwargs.get("crval"),
            # )

            # ps = pixel_scale
            # ps_ratio = pixel_scale_ratio

            output_wcs = {
                "wcs": wcs,
                # "pixel_scale": 3600.0 * ps,
                # "pixel_scale_ratio": ps_ratio,
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

    @property
    def compute_exptime(self):
        return self._compute_exptime

    @property
    def blend_meta(self):
        return self._blend_meta

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
        super().add_model(model_dict)
        if self.blend_meta:
            self._meta_blender.blend(model)
        if self.compute_exptime:
            self._resample_exptime(model)

    def _resample_exptime(self, model):
        if not hasattr(self, "_exptime_resampler"):
            self._exptime_resampler = ExptimeResampler(
                self.output_wcs, self.output_array_shape, self.good_bits, self.kernel
            )
        self._exptime_resampler.add_image(model)

    def finalize(self):
        super().finalize()

        if self.blend_meta:
            output_model = self._meta_blender.finalize()
        else:
            output_model = maker_utils.mk_datamodel(
                datamodels.MosaicModel, n_images=0, shape=(0, 0)
            )

        output_model.meta.resample.good_bits = self.good_bits
        output_model.meta.resample.weight_type = self.weight_type
        output_model.meta.resample.pixfrac = self.output_model["pixfrac"]
        output_model.meta.basic.product_type = "TBD"

        pixel_scale_ratio = self.output_model["pixel_scale_ratio"]
        if pixel_scale_ratio is not None:
            output_model.meta.resample.pixel_scale_ratio = pixel_scale_ratio

        output_model.meta.resample.pointings = self.output_model["pointings"]

        # copy over asn information
        output_model.meta.basic.location_name = self.input_models.asn.get(
            "target", "None"
        )
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
        if self.compute_exptime and hasattr(self, "_exptime_resampler"):
            exptime_total = self._exptime_resampler.finalize()
            m = exptime_total > 0
            total_exposure_time = np.mean(exptime_total[m]) if np.any(m) else 0
            max_exposure_time = np.max(exptime_total)
            log.info(
                f"Mean, max exposure times: {total_exposure_time:.1f}, "
                f"{max_exposure_time:.1f}"
            )
            output_model.meta.basic.mean_exposure_time = total_exposure_time
            output_model.meta.basic.max_exposure_time = max_exposure_time
            output_model.meta.resample.product_exposure_time = max_exposure_time

        if self._enable_ctx:
            output_model.context = self.output_model["con"].astype(np.uint32)

        for arr_name in ("err", "var_rnoise", "var_poisson", "var_flat"):
            if arr_name in self.output_model:
                new_array = self.output_model[arr_name]
                if new_array is not None:
                    setattr(output_model, arr_name, new_array)

        # assign wcs to output model
        assign_l3_wcs(output_model, self.output_wcs)

        return output_model

    def reset_arrays(self, n_input_models=None):
        super().reset_arrays(n_input_models)
        if self.blend_meta:
            self._meta_blender = MetaBlender()

    def resample_group(self, indices):
        if self.is_finalized():
            self.reset_arrays(len(indices))

        with self.input_models:
            for index in indices:
                model = self.input_models.borrow(index)
                self.add_model(model)
                self.input_models.shelve(model, index, modify=False)

        return self.finalize()

import logging

import numpy as np
from roman_datamodels import datamodels, dqflags
from stcal.resample import Resample

import romancal.skycell.skymap as sc
from romancal.lib.basic_utils import compute_var_rnoise

from .exptime_resampler import ExptimeResampler
from .l3_wcs import assign_l3_wcs
from .meta_blender import MetaBlender
from .resample_utils import compute_var_sky, make_output_wcs

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
        resample_on_skycell,
        wcs_kwargs=None,
        variance_array_names=None,
        propagate_dq=False,
    ):
        """
        See the base class `stcal.resample.resample.Resample` for more details.

        This is primarily used by `romancal.resample.ResampleStep`. See the step
        spec for parameter defaults and more information.

        Parameters
        ----------
        input_models : ModelLibrary
            The library of models to resample.

        output_wcs : dict, None
            Output wcs to resample to. If None one will be computed from
            input_models.

        pixfrac : float
            The fraction of a pixel that the pixel flux is confined to.

        kernel : str
            The name of the resampling kernel.

        fillval : float, None, str
            The value of output pixels that did not have contributions from
            input images' pixels.

        weight_type : str
            The weighting type for adding data.

        good_bits : int, str, None
            The bit mask of pixels to use for resampling.

        enable_ctx : bool
            If `True` include a context array in the resampled output.

        enable_var : bool
            If `True` resample variance arrays.

        compute_err : {"from_var", "driz_err"}, None
            If `None` do not compute an output model error array. For
            other options see `stcal.resample.resample.Resample`.

        compute_exptime : bool
            If `True` resample input exposure times to compute output
            exposure time information.

        blend_meta : bool
            If `True` include blended input model metadata in the output
            model.

        resample_on_skycell : bool
            If `True` and the association contains skycell information
            use the skycell wcs.

        wcs_kwargs : dict, None
            A dictionary of custom WCS parameters used to define an
            output WCS from input models' outlines. This argument is ignored
            when ``output_wcs`` is specified. See
            `romancal.resample.resample_utils.make_output_wcs`
            for supported options.

        propagate_dq : bool
            If `True`, propagate DQ during resampling. DQ flags are propagated
            by bitwise OR of all input DQ flags that contribute
            to a given output pixel.

        variance_array_names : list, None
            Variance arrays to resample.  If None, use stcal default.
        """
        # fillval indef doesn't define a starting value
        # since we're not resampling onto an existing array
        # we overwrite this to 0 (to match the old behavior)
        if fillval.lower() == "indef":
            fillval = 0
        self.input_models = input_models

        self._compute_exptime = compute_exptime
        self._blend_meta = blend_meta

        if output_wcs is None and resample_on_skycell:
            # first try to retrieve a sky cell name from the association
            try:
                skycell = sc.SkyCells.from_asns([self.input_models.asn])

                log.info(f"Skycell record: {skycell.data}")

                log.info(
                    f"Creating skycell image at ra: {skycell.radec_centers[0, 0]}  dec {skycell.radec_centers[0, 1]}",
                )
                log.info("Resampling to skycell wcs")
                output_wcs = {"wcs": skycell.wcs[0]}
            except ValueError as err:
                log.warning(f"Unable to compute skycell from input association: {err}")
                log.warning("Computing output wcs from all input wcses")

        if output_wcs is None:
            if wcs_kwargs is None:
                wcs_kwargs = {}
            wcs, ps, ps_ratio = make_output_wcs(input_models, **wcs_kwargs)
            output_wcs = {
                "wcs": wcs,
                "pixel_scale": ps,
                "pixel_scale_ratio": ps_ratio,
            }

        if variance_array_names:
            self.variance_array_names = variance_array_names

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
            propagate_dq=propagate_dq,
        )

    @property
    def compute_exptime(self):
        """Indicates if exposure time resampling is enabled."""
        return self._compute_exptime

    @property
    def blend_meta(self):
        """Indicates if metadata blending is enabled."""
        return self._blend_meta

    def _input_model_to_dict(self, model):
        """Convert an input datamodel to a dictionary suitable for the base class"""
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
            "var_rnoise": compute_var_rnoise(model),
            "var_poisson": model.var_poisson,
            "var_sky": compute_var_sky(model)
            if self.weight_type == "ivm-sky"
            else None,
            "err": model.err,
            "filename": model.meta.filename,
            "wcs": model.meta.wcs,
            "pixelarea_steradians": pixel_area,
            "group_id": self.input_models._model_to_group_id(model),
            "measurement_time": None,  # falls back to exposure_time
            "exposure_time": exposure_time,
            "start_time": model.meta.exposure.start_time,
            "end_time": model.meta.exposure.end_time,
            "duration": None,  # unused
            "level": level,
            "subtracted": subtracted,
            "effective_exposure_time": model.meta.exposure.effective_exposure_time,
            "var_flat": getattr(model, "var_flat", None),
        }

    def _get_intensity_scale(self, model):
        # always provide 1: https://github.com/spacetelescope/romancal/issues/1637
        return 1

    def add_model_hook(
        self,
        model,
        pixmap,
        pixel_scale_ratio,
        iscale,
        weight_map,
        xmin,
        xmax,
        ymin,
        ymax,
    ):
        if self.compute_exptime:
            if not hasattr(self, "_exptime_resampler"):
                self._exptime_resampler = ExptimeResampler(
                    self.output_wcs,
                    self.output_array_shape,
                    self.good_bits,
                    self.kernel,
                )
            self._exptime_resampler.add_image(
                model=model,
                pixmap=pixmap,
                pixel_scale_ratio=pixel_scale_ratio,
                xmin=xmin,
                xmax=xmax,
                ymin=ymin,
                ymax=ymax,
            )

    def add_model(self, model):
        model_dict = self._input_model_to_dict(model)
        super().add_model(model_dict)
        if self.blend_meta:
            self._meta_blender.blend(model)

    def finalize(self):
        super().finalize()

        if self.blend_meta:
            output_model = self._meta_blender.finalize()
        else:
            output_model = datamodels.MosaicModel.create_minimal()

        # copy over asn information
        output_model.meta.association.name = self.input_models.asn.get(
            "table_name", "?"
        )

        # resample parameters
        output_model.meta.resample.good_bits = self.good_bits
        output_model.meta.resample.pixel_scale_ratio = self.output_model[
            "pixel_scale_ratio"
        ]
        output_model.meta.resample.pixfrac = self.output_model["pixfrac"]
        output_model.meta.resample.pointings = self.output_model["pointings"]
        output_model.meta.resample.weight_type = self.weight_type

        # every resampling will generate these
        output_model.data = self.output_model["data"]
        output_model.weight = self.output_model["wht"]

        # some things are conditional
        if self.compute_exptime and hasattr(self, "_exptime_resampler"):
            exptime_total = self._exptime_resampler.finalize()
            m = exptime_total > 0
            total_exposure_time = (
                np.mean(exptime_total[m], dtype="f8") if np.any(m) else 0
            )
            max_exposure_time = np.max(exptime_total)
            log.info(
                f"Mean, max exposure times: {total_exposure_time:.1f}, "
                f"{max_exposure_time:.1f}"
            )
            output_model.meta.coadd_info.max_exposure_time = max_exposure_time
            output_model.meta.coadd_info.exposure_time = total_exposure_time

        if self.enable_ctx:
            output_model.context = self.output_model["con"].astype(np.uint32)

        if self.propagate_dq:
            output_model.dq = self.output_model["dq"].astype(np.uint32)

        for arr_name in ["err", *self.variance_array_names]:
            if arr_name in self.output_model:
                new_array = self.output_model[arr_name]
                if new_array is not None:
                    setattr(output_model, arr_name, new_array)

        # assign wcs to output model
        assign_l3_wcs(output_model, self.output_wcs)
        output_model.meta.wcsinfo.skycell_name = self.input_models.asn.get(
            "target", "None"
        )
        # get data release ID
        output_model.meta.data_release_id = self.input_models.asn.get(
            "data_release_id", "p"
        )

        return output_model

    def reset_arrays(self, n_input_models=None):
        super().reset_arrays(n_input_models)
        if self.blend_meta:
            self._meta_blender = MetaBlender()

    def resample_group(self, indices):
        """
        Resample the models at indices.

        Parameters
        ----------
        indices : iterable
            An iterable that returns the indices of models to resample.

        Returns
        -------
        MosaicModel
            The resampled datamodel
        """
        if self.is_finalized():
            self.reset_arrays(len(indices))

        with self.input_models:
            for index in indices:
                model = self.input_models.borrow(index)
                self.add_model(model)
                self.input_models.shelve(model, index, modify=False)

        return self.finalize()

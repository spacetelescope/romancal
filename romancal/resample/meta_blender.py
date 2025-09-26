from typing import ClassVar

from astropy.time import Time
from roman_datamodels import datamodels, stnode


class MetaBlender:
    _meta_blend_paths: ClassVar[dict[str, str | None]] = {
        "observation.execution_plan": None,
        "observation.exposure": None,
        "observation.observation": None,
        "observation.pass": None,
        "observation.program": None,
        "observation.segment": None,
        "observation.visit": None,
        "program.title": None,
        "program.investigator_name": None,
        "program.category": None,
        "program.subcategory": None,
        "program.science_category": None,
        "instrument.optical_element": None,
        "cal_step.flux": "INCOMPLETE",
        "cal_step.outlier_detection": "INCOMPLETE",
        "cal_step.skymatch": "INCOMPLETE",
    }

    def _blend_first(self, model):
        # make a blank mosic metdata node
        self._model = datamodels.MosaicModel.create_minimal()

        # FIXME assuming everything is a prompt coadd
        self._model.meta.product_type = stnode.ProductType("p_visit_coadd")

        self._meta = self._model.meta

        self._model["individual_image_cal_logs"] = []

        # for computing mean start time
        self._start_times = []

        # grab start/end times, blended with all other models below
        self._meta.coadd_info.time_first = model.meta.exposure.start_time
        self._meta.coadd_info.time_last = model.meta.exposure.end_time

        # copy over metadata from first model
        for path, default in self._meta_blend_paths.items():
            dst = self._meta
            src = model.meta
            *sub_path, key = path.split(".")
            for k in sub_path:
                src = src.get(k, {})
                dst = dst[k]
            dst[key] = src.get(key, default)

    def blend(self, model):
        if not hasattr(self, "_model"):
            self._blend_first(model)
        else:
            # for non-first only blending
            self._meta.coadd_info.time_first = min(
                self._meta.coadd_info.time_first, model.meta.exposure.start_time
            )
            self._meta.coadd_info.time_last = max(
                self._meta.coadd_info.time_last, model.meta.exposure.end_time
            )

            # blend
            for path, default in self._meta_blend_paths.items():
                dst = self._meta
                src = model.meta
                *sub_path, key = path.split(".")
                for k in sub_path:
                    src = src.get(k, {})
                    dst = dst[k]
                if key not in src:
                    continue
                if dst[key] != src[key]:
                    dst[key] = default

        self._start_times.append(model.meta.exposure.start_time)

        if cal_logs := model.meta.get("cal_logs"):
            self._model["individual_image_cal_logs"].append(cal_logs)

        self._meta.resample.members.append(model.meta.filename)
        self._model.add_image(model)

    def finalize(self):
        self._meta.coadd_info.time_mean = Time(self._start_times).mean()
        if all(
            getattr(self._meta.observation, key) is not None
            for key in ("execution_plan", "pass", "segment", "observation", "visit")
        ):
            self._meta.observation.exposure_grouping = "v{execution_plan:02d}{pass:03d}{segment:03d}{observation:03d}{visit:03d}".format(
                **self._meta.observation
            )
        else:
            self._meta.observation.exposure_grouping = None

        self._model.populate_individual_image_meta()

        return self._model

import numpy as np
from roman_datamodels import datamodels, maker_utils


class MetaBlender:
    def _blend_first(self, model):
        # make a blank mosic metdata node
        self._model = maker_utils.mk_datamodel(
            datamodels.MosaicModel, n_images=0, shape=(0, 0)
        )
        self._meta = self._model.meta

        self._model["individual_image_cal_logs"] = []

        self._mid_mjds = []

        self._meta.basic.time_first_mjd = model.meta.exposure.start_time.mjd
        self._meta.basic.time_last_mjd = model.meta.exposure.end_time.mjd

        # copy some metadata from the first input model
        # FIXME this is incorrect and these values should be
        # "blended" from all models
        self._meta.basic.visit = model.meta.observation.visit
        self._meta.basic.segment = model.meta.observation.segment
        self._meta.basic["pass"] = model.meta.observation["pass"]
        self._meta.basic.program = model.meta.observation.program
        self._meta.basic.optical_element = model.meta.instrument.optical_element
        self._meta.basic.instrument = model.meta.instrument.name
        self._meta.coordinates = model.meta.coordinates
        self._meta.program = model.meta.program
        for step_name in self._meta.cal_step:
            if hasattr(model.meta.cal_step, step_name):
                setattr(
                    self._meta.cal_step,
                    step_name,
                    getattr(model.meta.cal_step, step_name),
                )

    def blend(self, model):
        if not hasattr(self, "_model"):
            self._blend_first(model)
        else:
            # for non-first only blending
            self._meta.basic.time_first_mjd = min(
                self._meta.basic.time_first_mjd, model.meta.exposure.start_time.mjd
            )
            self._meta.basic.time_last_mjd = max(
                self._meta.basic.time_last_mjd, model.meta.exposure.end_time.mjd
            )

        self._model["individual_image_cal_logs"].append(model.meta.cal_logs)
        self._meta.resample.members.append(model.meta.filename)

        self._mid_mjds.append(model.meta.exposure.mid_time.mjd)

        # FIXME this is an existing hack that we have to reproduce here
        # since roman_datamodels was never updated
        cal_logs = model.meta.cal_logs
        del model.meta["cal_logs"]
        try:
            self._model.append_individual_image_meta(model.meta)
        finally:
            model.meta.cal_logs = cal_logs

    def finalize(self):
        self._meta.basic.time_mean_mjd = np.mean(self._mid_mjds)
        return self._model

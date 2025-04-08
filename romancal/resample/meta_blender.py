from collections import defaultdict

import numpy as np
from asdf.lazy_nodes import AsdfDictNode, AsdfListNode
from astropy.table import QTable
from roman_datamodels import datamodels, maker_utils, stnode


class MissingCellType:
    pass


MISSING_CELL = MissingCellType()


class TableBuilder:
    def __init__(self):
        self._columns = {}
        self._nrows = 0

    def _get_column(self, name):
        if name not in self._columns:
            self._columns[name] = [MISSING_CELL] * self._nrows
        return self._columns[name]

    def _santize_value(self, value):
        if isinstance(value, list | dict | AsdfDictNode | AsdfListNode):
            return str(value)
        return value

    def add_row(self, data_dict):
        updated = set()
        for key, value in data_dict.items():
            self._get_column(key).append(self._santize_value(value))
            updated.add(key)
        for key in self._columns.keys() - updated:
            self._columns[key].append(MISSING_CELL)

    def to_table(self):
        # TODO fix None and MISSING_CELL values
        return QTable(self._columns)


class MetaBlender:
    def __init__(self):
        self._tables = defaultdict(lambda: TableBuilder())

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

    def _update_tables(self, meta):
        basic_data = {}
        for key, value in meta.to_flat_dict().items():
            if key in {"wcs", "cal_logs"}:
                continue

            if isinstance(value, stnode.DNode):
                self._tables[key].add_row(value)
            else:
                basic_data[key] = value
        if basic_data:
            self._tables["basic"].add_row(basic_data)

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

        mid_time = (
            model.meta.exposure.start_time.mjd
            + model.meta.exposure.exposure_time / 60 / 60 / 24
        )
        self._mid_mjds.append(mid_time)

        self._update_tables(model.meta)

    def finalize(self):
        self._meta.basic.time_mean_mjd = np.mean(self._mid_mjds)
        self._meta.individual_image_meta = stnode.IndividualImageMeta()
        for table_name, builder in self._tables.items():
            self._meta.individual_image_meta[table_name] = builder.to_table()
        return self._model

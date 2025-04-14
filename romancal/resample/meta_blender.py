from collections import defaultdict

import numpy as np
from asdf.lazy_nodes import AsdfDictNode, AsdfListNode
from asdf.tags.core.ndarray import NDArrayType
from astropy.table import Table
from roman_datamodels import datamodels, maker_utils, stnode


class MissingCellType:
    def __str__(self):
        return "MISSING_CELL"


MISSING_CELL = MissingCellType()


class TableBuilder:
    def __init__(self):
        self._columns = {}
        self._flagged_columns = set()

    def _get_column(self, name, row_index):
        if name not in self._columns:
            if row_index:
                self._flagged_columns.add(name)
            self._columns[name] = [MISSING_CELL] * row_index
        return self._columns[name]

    def _sanitize_value(self, value):
        if isinstance(value, list | dict | AsdfDictNode | AsdfListNode):
            return str(value)
        return value

    def add_row(self, data_dict, row_index):
        updated = set()
        for key, value in data_dict.items():
            svalue = self._sanitize_value(value)
            if svalue is None:
                self._flagged_columns.add(key)
            self._get_column(key, row_index).append(svalue)
            updated.add(key)
        for key in self._columns.keys() - updated:
            self._flagged_columns.add(key)
            self._columns[key].append(MISSING_CELL)

    def to_table(self):
        array_columns = {}
        for name, values in self._columns.items():
            if name not in self._flagged_columns:
                array_columns[name] = values
                continue
            # this column has either a None or MISSING_CELL which
            # has to be filled in
            arr = np.array([v for v in values if v not in (None, MISSING_CELL)])

            # if we have no valid values, return all "None"
            if not arr.size:
                array_columns[name] = ["None"] * len(values)
                continue

            # if we have a float or int use 'nan' for missing/None values
            # this will convert int columns to floats
            if arr.dtype.kind in ("f", "i"):
                array_columns[name] = [
                    np.nan if v in (MISSING_CELL, None) else v for v in values
                ]
                continue

            # if all else fails, convert everything to a string
            array_columns[name] = [str(v) for v in values]

        return Table(array_columns)


class MetaBlender:
    def __init__(self):
        self._tables = defaultdict(lambda: TableBuilder())
        self._n_rows = 0

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
                self._tables[key].add_row(value, self._n_rows)
            elif not isinstance(value, dict | NDArrayType | Table | AsdfDictNode):
                basic_data[key] = value
        if basic_data:
            self._tables["basic"].add_row(basic_data, self._n_rows)
        self._n_rows += 1

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

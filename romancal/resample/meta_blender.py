from collections import defaultdict
from typing import ClassVar

import numpy as np
from asdf.lazy_nodes import AsdfDictNode, AsdfListNode
from asdf.tags.core.ndarray import NDArrayType
from astropy.table import Table
from astropy.time import Time
from roman_datamodels import datamodels, stnode


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
    _meta_blend_paths: ClassVar = {
        "observation.program": None,
        "observation.execution_plan": None,
        "observation.pass": None,
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

    def __init__(self):
        self._tables = defaultdict(TableBuilder)
        self._n_rows = 0

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

        self._update_tables(model.meta)

    def finalize(self):
        self._meta.coadd_info.time_mean = Time(self._start_times).mean()
        # TODO observation.observation? what if some of these values are None?
        self._meta.observation.exposure_grouping = None
        # (
        #     "v{execution_plan:02d}{pass:03d}{segment:03d}001{visit:03d}".format(
        #         **self._meta.observation
        #     )
        # )
        self._meta.individual_image_meta = {}
        for table_name, builder in self._tables.items():
            self._meta.individual_image_meta[table_name] = builder.to_table()
        return self._model

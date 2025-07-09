from collections import defaultdict

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
    def __init__(self):
        self._tables = defaultdict(TableBuilder)
        self._n_rows = 0

    def _blend_first(self, model):
        # make a blank mosic metdata node
        self._model = datamodels.MosaicModel.create_minimal()
        #     {
        #         "meta": {
        #             "file_date": Time(
        #                 "2020-01-01T00:00:00.0", format="isot", scale="utc"
        #             ),
        #         },
        #     }
        # )

        # FIXME assuming everything is a prompt coadd
        self._model.meta.product_type = stnode.ProductType("p_visit_coadd")

        self._meta = self._model.meta

        self._model["individual_image_cal_logs"] = []

        # for computing mean start time
        self._start_times = []

        # for computing combined optical element
        self._optical_elements = set()

        # grab start/end times, blended with all other models below
        self._meta.coadd_info.time_first = model.meta.exposure.start_time
        self._meta.coadd_info.time_last = model.meta.exposure.end_time

        # copy some metadata from the first input model
        # FIXME this is incorrect and these values should be
        # "blended" from all models
        for key in ("program", "execution_plan", "pass", "segment", "visit"):
            self._meta.observation[key] = model.meta.observation[key]
        for key in ("title", "investigator_name", "category", "subcategory", "science_category"):
            self._meta.program[key] = model.meta.program[key]
        # TODO we can't propagate L2 categories and subcategories
        self._meta.program.category = "CAL"
        self._meta.program.subcategory = "None"
        for step_name in ("flux", "outlier_detection", "skymatch"):
            self._meta.cal_step[step_name] = model.meta.get("cal_step", {}).get(step_name, "INCOMPLETE")

    def _update_tables(self, meta):
        # TODO unpass
        pass
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
                self._meta.coadd_info.time_first,
                model.meta.exposure.start_time
            )
            self._meta.coadd_info.time_last = max(
                self._meta.coadd_info.time_last,
                model.meta.exposure.end_time
            )

        self._start_times.append(model.meta.exposure.start_time)
        self._optical_elements.add(model.meta.instrument.optical_element)

        if cal_logs := model.meta.get("cal_logs"):
            self._model["individual_image_cal_logs"].append(cal_logs)

        self._meta.resample.members.append(model.meta.filename)

        self._update_tables(model.meta)

    def finalize(self):
        self._meta.coadd_info.time_mean = Time(self._start_times).mean()
        self._meta.instrument.optical_element = ", ".join(sorted(self._optical_elements))
        # TODO observation.observation?
        self._meta.observation.exposure_grouping = "v{execution_plan:02d}{pass:03d}{segment:03d}001{visit:03d}".format(**self._meta.observation)
        self._meta.individual_image_meta = {}
        # TODO unpass
        # for table_name, builder in self._tables.items():
        #     self._meta.individual_image_meta[table_name] = builder.to_table()
        return self._model

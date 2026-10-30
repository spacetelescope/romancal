import copy
import numbers
import pprint
from collections.abc import Mapping, Sequence, Set
from dataclasses import dataclass
from textwrap import dedent

import asdf
import numpy as np
import roman_datamodels as rdm
from astropy.table import Table
from astropy.time import Time
from gwcs import WCS
from gwcs.wcstools import grid_from_bounding_box


def _wcs_to_ra_dec(wcs):
    x, y = grid_from_bounding_box(wcs.bounding_box)
    return wcs(x, y)


def _get_model_info(tree):
    node = tree.get("roman", {})
    if "roman" not in tree:
        model_type = "Unknown"
    else:
        # look up datamodel from node
        model_type = str(type(rdm.open(asdf.AsdfFile({"roman": node}))))
    crds_context = (
        node.get("meta", {})
        .get("ref_file", {})
        .get("crds", {})
        .get("context", "Unknown context")
    )
    filename = node.get("meta", {}).get("filename", "Unknown filename")
    return {
        "model type": model_type,
        "crds_context": crds_context,
        "filename": filename,
    }


@dataclass
class DiffSettings:
    rtol: float
    atol: float
    equal_nan: bool
    ignore_paths: set[str]


class Diff:
    def __init__(self, truth, result, settings):
        # capture file information: model type, crds context
        self._truth_info = _get_model_info(truth)
        self._result_info = _get_model_info(result)
        self.settings = settings
        # keys=diff types: arrays_differ, etc
        # values=dict of {path: {result: ..., truth: ...}}
        self.diff = {}
        self._seen = set()
        self._diff(truth, result, ())

    @property
    def identical(self):
        return not bool(self.diff)

    def _skip_path(self, path):
        return ".".join(path) in self.settings.ignore_paths

    def _record(self, diff_type, path, details):
        if diff_type not in self.diff:
            self.diff[diff_type] = {}
        self.diff[diff_type][".".join(path)] = details

    def _diff_types(self, truth, result, path):
        if type(truth) is type(result):
            return
        if not {type(truth), type(result)} - {asdf.tags.core.NDArrayType, np.ndarray}:
            # allow NDArrayType == ndarray
            return
        return {"truth": type(truth), "result": type(result)}

    def _diff_wcses(self, truth, result, path):
        # for comparing wcs instances this function evaluates
        # each wcs and compares the resulting ra and dec outputs
        ra_truth, dec_truth = _wcs_to_ra_dec(truth)
        ra_result, dec_result = _wcs_to_ra_dec(result)
        details = {}
        if ra_diff := self._diff_arrays(ra_truth, ra_result, path):
            details["ra"] = ra_diff
        if dec_diff := self._diff_arrays(dec_truth, dec_result, path):
            details["dec"] = dec_diff
        return details

    def _diff_tables(self, truth, result, path):
        # compare as arrays
        if array_difference := self._diff_arrays(
            truth.as_array(), result.as_array(), path
        ):
            details = array_difference
        else:
            details = {}

        # compare metas (ignoring some specific paths)
        meta_settings = copy.copy(self.settings)
        meta_settings.ignore_paths = {"date", "version", "versions"}
        meta_diff = Diff(truth.meta, result.meta, meta_settings)
        if not meta_diff.identical:
            details["metas_differ"] = meta_diff.diff
        return details

    def _diff_times(self, truth, result, path):
        if truth == result:
            return
        return {"truth": truth, "result": result}

    def _diff_arrays(self, truth, result, path):
        details = {}
        for attr in ("shape", "dtype", "unit"):
            truth_attr = getattr(truth, attr, None)
            result_attr = getattr(result, attr, None)
            if truth_attr != result_attr:
                details[attr] = {"truth": truth_attr, "result": result_attr}
        if truth.dtype.names or result.dtype.names:
            if truth.dtype.names != result.dtype.names:
                details["names"] = {
                    "truth": truth.dtype.names,
                    "result": result.dtype.names,
                }
            if not details:
                per_column_details = {}
                # only compare values if no other difference exists
                for name in truth.dtype.names:
                    child = (*path, name)
                    if self._skip_path(child):
                        continue
                    if column_details := self._diff_arrays(
                        truth[name], result[name], child
                    ):
                        per_column_details[name] = column_details
                if per_column_details:
                    details["column_values"] = per_column_details
            return details

        if details:
            # don't do more comparisons if we have differences
            return details

        if not np.issubdtype(truth.dtype, np.number):
            # non-numeric (e.g. string, bool) arrays must match exactly
            if n_diffs := np.count_nonzero(truth != result):
                details["n_diffs"] = n_diffs
            return details

        n_diffs = np.count_nonzero(
            np.logical_not(
                np.isclose(
                    truth,
                    result,
                    rtol=self.settings.rtol,
                    atol=self.settings.atol,
                    equal_nan=self.settings.equal_nan,
                )
            )
        )

        if n_diffs == 0:
            return

        details["n_diffs"] = n_diffs

        abs_diff = np.abs(truth - result)
        # check for all nan difference to avoid errors in other comparisons
        if not np.any(np.isfinite(abs_diff)):
            details["all_nonfinite_diffs"] = True
            return details

        index = np.unravel_index(np.nanargmax(abs_diff), truth.shape)
        details["worst_abs_diff"] = {
            "index": index,
            "value": abs_diff[index],
        }
        with np.errstate(invalid="ignore", divide="ignore"):
            # 0 / 0 == nan and produces an 'invalid' error
            # 1 / 0 == inf and produces a 'divide' error
            # ignore these here for computing the fractional diff
            fractional_diff = np.abs(truth / result)

        if np.any(np.isfinite(abs_diff)):
            index = np.unravel_index(np.nanargmax(fractional_diff), truth.shape)
            details["worst_fractional_diff"] = {
                "index": index,
                "value": fractional_diff[index],
            }
        # np.quantiles(new - old, [0.0001, 0.001, 0.01, 0.1, 0.16, 0.5, 0.84, 0.9, 0.99, 0.999, 0.9999])
        return details

    def _diff_mappings(self, truth, result, path):
        for key in truth:
            child = (*path, key)
            if self._skip_path(child):
                continue
            if key in result:
                self._diff(truth[key], result[key], child)
            else:
                self._record("dictionary_item_removed", child, str(truth[key]))
        for key in result:
            child = (*path, key)
            if self._skip_path(child):
                continue
            if key not in truth:
                self._record("dictionary_item_added", child, str(result[key]))

    def _diff_sequences(self, truth, result, path):
        for index in range(max(len(truth), len(result))):
            # handle integer path elements
            child = (*path, str(index))
            if self._skip_path(child):
                continue
            try:
                truth_value = truth[index]
            except LookupError:
                self._record("iterable_item_added", child, str(result[index]))
                continue
            try:
                result_value = result[index]
            except LookupError:
                self._record("iterable_item_removed", child, str(truth_value))
                continue
            self._diff(truth_value, result_value, child)

    def _diff_sets(self, truth, result, path):
        removed = truth - result
        added = result - truth
        if removed:
            self._record("set_item_removed", path, removed)
        if added:
            self._record("set_item_added", path, added)

    def _diff_numbers(self, truth, result, path):
        if np.isclose(
            truth,
            result,
            rtol=self.settings.rtol,
            atol=self.settings.atol,
            equal_nan=self.settings.equal_nan,
        ):
            return
        return {"truth": truth, "result": result}

    def _diff_default(self, truth, result, path):
        if truth == result:
            return
        # convert values to strings to avoid carrying object references through the default path
        return {"truth": str(truth), "result": str(result)}

    def _diff(self, truth, result, path):
        if self._skip_path(path):
            return

        # avoid repeats
        # TODO fix to avoid marking scalars as seen
        # TODO fix issue where an object appearing in >1 path (1 ignored) can result in false equality (is this an issue?)
        seen_key = (id(truth), id(result))
        if seen_key in self._seen:
            return
        self._seen.add(seen_key)

        # compare types
        if details := self._diff_types(truth, result, path):
            self._record("types_differ", path, details)
            return

        # at this point, can assume types are the same
        if isinstance(truth, WCS):
            if details := self._diff_wcses(truth, result, path):
                self._record("wcs_differ", path, details)
        elif isinstance(truth, Table):
            if details := self._diff_tables(truth, result, path):
                self._record("tables_differ", path, details)
        elif isinstance(truth, Time):
            if details := self._diff_times(truth, result, path):
                self._record("times_differ", path, details)
        elif isinstance(truth, (np.ndarray, asdf.tags.core.NDArrayType)):
            if details := self._diff_arrays(truth, result, path):
                self._record("arrays_differ", path, details)
        elif isinstance(truth, Mapping):
            self._diff_mappings(truth, result, path)
        elif isinstance(truth, Sequence) and not isinstance(truth, (str, bytes)):
            self._diff_sequences(truth, result, path)
        elif isinstance(truth, Set):
            self._diff_sets(truth, result, path)
        elif isinstance(truth, numbers.Number):
            if details := self._diff_numbers(truth, result, path):
                self._record("values_differ", path, details)
        else:
            if details := self._diff_default(truth, result, path):
                self._record("values_differ", path, details)

    def report(self):
        title = dedent(
            f"""
            Diff report for:
                result file: {self._result_info["filename"]}
                    model type: {self._result_info["model type"]}
                    crds_context: {self._result_info["crds_context"]}
                truth file: {self._truth_info["filename"]}
                    model type: {self._truth_info["model type"]}
                    crds_context: {self._truth_info["crds_context"]}
            """
        )
        report = pprint.pformat(self.diff)
        return f"{title}\n{report}"

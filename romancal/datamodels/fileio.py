import warnings

import roman_datamodels.datamodels as rdm

from . import filetype
from .library import ModelLibrary
from .migration import _error_on_downgrade, update_model_version

__all__ = ["open_dataset"]


def open_dataset(
    dataset,
    *,
    update_version=False,
    downgrade_unknown_version=False,
    return_type=False,
    as_library=False,
    open_kwargs=None,
):
    """
    If needed, open the provided dataset for processing.

    This is used for handling the varied Step input types
    converting the provided input to the format expected by
    the Step as defined by the provided arguments.

    Parameters
    ----------
    dataset : str, Path, DataModel, ModelLibrary or other
        Varied input format. For str/Path open as DataModel
        for ASDF files and ModelLibrary for JSON. For
        DataModel/ModelLibrary return either unmodified
        or converted (based on the other parameters). For
        other formats (for example a list of filenames)
        the dataset is opened as a ModelLibrary.

    update_version : bool, optional
        Update the dataset to the newest DataModel (tag) version.

    downgrade_unknown_version : bool, optional
        Process newer (unknown tag) DataModels as the newest known version.

    return_type : bool, optional
        Also return the input dataset type (as returned by filetype.check).

    as_library : bool, optional
        Convert the provided dataset to a ModelLibrary instance.

    open_kwargs : dict, optional
        Keyword arguments to pass to DataModel or ModelLibrary opening.

    Returns
    -------
    opened_dataset : DataModel or ModelLibrary
        The opened dataset. If the provided dataset was already
        open and the correct version it will be returned unchanged.

    input_dataset_type : str, optional
        The type of the input dataset as determined by filetype.check
    """
    open_kwargs = open_kwargs or {}

    try:
        dataset_type = filetype.check(dataset)
    except ValueError:
        # to allow ModelLibrary to handle lists of filenames and lists of datamodels
        dataset_type = "unknown"

    match dataset_type:
        case "DataModel":
            if update_version:
                result = update_model_version(dataset)
            else:
                result = dataset

        case "ModelLibrary":
            dataset._datamodels_open_kwargs["update_version"] = update_version
            dataset._datamodels_open_kwargs["downgrade_unknown_version"] = (
                downgrade_unknown_version
            )
            result = dataset

        case "asn":
            result = ModelLibrary(
                dataset,
                update_version=update_version,
                downgrade_unknown_version=downgrade_unknown_version,
                **open_kwargs,
            )

        case "asdf":
            if open_kwargs.pop("on_disk", None) is not None:
                warnings.warn(
                    "on_disk is only supported for associations, ignoring on_disk",
                    stacklevel=2,
                )
            with _error_on_downgrade(not downgrade_unknown_version):
                model = rdm.open(dataset, **open_kwargs)
            if update_version:
                result = update_model_version(model, close_on_update=True)
            else:
                result = model

        case _:
            # on disk is only supported for associations
            if open_kwargs.get("on_disk", False):
                kwargs = open_kwargs.copy()
                kwargs.pop("on_disk")
                warnings.warn(
                    f"ModelLibrary on_disk is not supported for input {dataset}. "
                    "Disabling on_disk mode.",
                    stacklevel=2,
                )
            else:
                kwargs = open_kwargs
            result = ModelLibrary(
                dataset,
                update_version=update_version,
                downgrade_unknown_version=downgrade_unknown_version,
                **kwargs,
            )

    if as_library and isinstance(result, rdm.DataModel):
        result = ModelLibrary(
            [result],
            update_version=update_version,
            downgrade_unknown_version=downgrade_unknown_version,
            **open_kwargs,
        )

    if return_type:
        return result, dataset_type

    return result

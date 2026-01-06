import roman_datamodels.datamodels as rdm

from . import filetype
from .library import ModelLibrary
from .migration import update_model_version


def open_dataset(
    dataset,
    *,
    update_version=False,
    return_type=False,
    as_library=False,
    open_kwargs=None,
):
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
            result = dataset
        case "asn":
            result = ModelLibrary(dataset, update_version=update_version, **open_kwargs)
        case "asdf":
            model = rdm.open(dataset, **open_kwargs)
            if update_version:
                result = update_model_version(model, close_on_update=True)
            else:
                result = model
        case _:
            result = ModelLibrary(dataset, update_version=update_version, **open_kwargs)

    if as_library and isinstance(result, rdm.DataModel):
        result = ModelLibrary([result], update_version=update_version, **open_kwargs)

    if return_type:
        return result, dataset_type
    return result

"""
Implementation of the datamodels.open function.
"""
from pathlib import PurePath, Path

import asdf

from .core import RomanDataModel
from .reference_files.referencefile import ReferenceFileModel
from .reference_files.flat import FlatModel


def open(init, memmap=False, **model_kwargs):
    """
    Open an ASDF file and wrap it in a data model class.

    Parameters
    ----------
    init : str or pathlib.PurePath or asdf.AsdfFile or romancal.datamodels.RomanDataModel
        If str or `~pathlib.PurePath`, filesystem path to an ASDF file.
        If `~asdf.AsdfFile`, an open ASDF file.
        If `~romancal.datamodels.RomanDataModel`, an already open model.

    memmap : bool, optional
        Set to True to enable memory-mapping of arrays.  Ignored if the
        ``init`` argument is an AsdfFile.

    model_kwargs : dict, optional
        Additional arguments to pass when initializing the model.

    Returns
    -------
    model : RomanDataModel
        Instance of `RomanDataModel` or subclass.
    """
    if isinstance(init, (str, PurePath)):
        if isinstance(init, str):
            init = Path(init)

        if not init.is_file():
            raise FileNotFoundError(f"File at path '{init}' does not exist")

        asdf_kwargs = {
            # The copy_arrays argument instructs the asdf library to
            # make in-memory copies of arrays, so we need to send it
            # the opposite of our memmap value.
            "copy_arrays": not memmap
        }

        asdf_file = asdf.open(init, **asdf_kwargs)
    elif isinstance(init, asdf.AsdfFile):
        asdf_file = init
    elif isinstance(init, RomanDataModel):
        # Make a copy of the model so that the original instance can
        # be closed without impacting this one.
        return init.__class__(init)
    else:
        raise TypeError("init must be a path to an ASDF file, an open AsdfFile instance, or a RomanDataModel")

    model_class = _select_model_class(asdf_file)

    model = model_class(asdf_file, **model_kwargs)

    return model


def _select_model_class(asdf_file):
    """
    Select a model class based on ASDF file metadata.

    Parameters
    ----------
    asdf_file : asdf.AsdfFile

    Returns
    -------
    type
        `RomanDataModel` or subclass.
    """
    # TODO: Still need to decide how this is going to work.
    # Options include:
    # - Choose according to model_type metadata field (similar to jwst)
    # - Infer from other relevant metadata (reference file type, exposure type, etc)
    # - Use ASDF tags to automatically instantiate the correct model
    # - Use RomanDataModel for all data but select the schema according to one of the above
    # - ???

    # Check for the existence of reftype to indicate a reference file model
    reftype = asdf_file["meta"].get("reftype")
    if reftype is not None:
        # Check if flat file model
        if reftype == "FLAT":
            return FlatModel
        # Return base reference file model
        else:
            return ReferenceFileModel

    # Not a reference file model
    else:
        return RomanDataModel

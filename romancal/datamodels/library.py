import copy
import os.path
import tempfile
from collections.abc import Iterable, MutableMapping, Sequence
from pathlib import Path
from types import MappingProxyType

import asdf
from roman_datamodels import open as datamodels_open

from romancal.associations import AssociationNotValidError, load_asn

__all__ = ["LibraryError", "BorrowError", "ClosedLibraryError", "ModelLibrary"]


class LibraryError(Exception):
    """
    Generic ModelLibrary related exception
    """

    pass


class BorrowError(LibraryError):
    """
    Exception indicating an issue with model borrowing
    """

    pass


class ClosedLibraryError(LibraryError):
    """
    Exception indicating a library method was used outside of a
    ``with`` context (that "opens" the library).
    """

    pass


class _Ledger(MutableMapping):
    def __init__(self):
        self._id_to_index = {}
        self._index_to_model = {}

    def __getitem__(self, model_or_index):
        if not isinstance(model_or_index, int):
            index = self._id_to_index[id(model_or_index)]
        else:
            index = model_or_index
        return self._index_to_model[index]

    def __setitem__(self, index, model):
        self._index_to_model[index] = model
        self._id_to_index[id(model)] = index

    def __delitem__(self, model_or_index):
        if isinstance(model_or_index, int):
            index = model_or_index
            model = self._index_to_model[index]
        else:
            model = model_or_index
            index = self._id_to_index[id(model)]
        del self._id_to_index[id(model)]
        del self._index_to_model[index]

    def __iter__(self):
        # only return indexes
        return iter(self._index_to_model)

    def __len__(self):
        return len(self._id_to_index)


class _OnDiskModelStore(MutableMapping):
    def __init__(self, memmap=False, directory=None):
        self._memmap = memmap
        if directory is None:
            self._tempdir = tempfile.TemporaryDirectory(dir="")
            self._path = Path(self._tempdir.name)
        else:
            self._path = Path(directory)
        self._filenames = {}

    def __getitem__(self, key):
        if key not in self._filenames:
            raise KeyError(f"{key} is not in {self}")
        return datamodels_open(self._filenames[key], memmap=self._memmap)

    def __setitem__(self, key, value):
        if key in self._filenames:
            fn = self._filenames[key]
        else:
            model_filename = value.meta.filename
            if model_filename is None:
                model_filename = "model.asdf"
            subpath = self._path / f"{key}"
            os.makedirs(subpath)
            fn = subpath / model_filename
            self._filenames[key] = fn

        # save the model to the temporary location
        value.save(fn)

    def __del__(self):
        if hasattr(self, "_tempdir"):
            self._tempdir.cleanup()

    def __delitem__(self, key):
        del self._filenames[key]

    def __iter__(self):
        return iter(self._filenames)

    def __len__(self):
        return len(self._filenames)


class ModelLibrary(Sequence):
    """
    A "library" of models (loaded from an association file).

    Do not anger the librarian!

    The library owns all models from the association and it will handle
    opening and closing files.

    Models can be "borrowed" from the library (by iterating through the
    library or indexing a specific model). However the library must be
    "open" (used in a ``with`` context)  to borrow a model and the model
    must be "returned" before the library "closes" (the ``with`` context exits).

    >>> with library:   # doctest: +SKIP
            model = library[0]  # borrow the first model
            # do stuff with the model
            library[0] = model  # return the model

    Failing to "open" the library will result in a ClosedLibraryError.

    Failing to "return" a borrowed model will result in a BorrowError.
    """

    def __init__(
        self,
        init,
        asn_exptypes=None,
        asn_n_members=None,
        on_disk=False,
        memmap=False,
        temp_directory=None,
    ):
        self._on_disk = on_disk
        self._open = False
        self._ledger = _Ledger()

        # FIXME is there a cleaner way to pass these along to datamodels.open?
        self._memmap = memmap

        if self._on_disk:
            self._model_store = _OnDiskModelStore(memmap, temp_directory)
        else:
            self._model_store = {}

        if isinstance(init, MutableMapping):
            asn_data = init
            self._asn_dir = os.path.abspath(".")
            self._asn = init

            if asn_exptypes is not None:
                raise NotImplementedError()

            if asn_n_members is not None:
                raise NotImplementedError()

            self._members = self._asn["products"][0]["members"]

            for member in self._members:
                if "group_id" not in member:
                    filename = os.path.join(self._asn_dir, member["expname"])
                    member["group_id"] = _file_to_group_id(filename)
        elif isinstance(init, (str, Path)):
            asn_path = os.path.abspath(os.path.expanduser(os.path.expandvars(init)))
            self._asn_dir = os.path.dirname(asn_path)

            # TODO asn_table_name is there another way to handle this
            self.asn_table_name = os.path.basename(asn_path)

            # load association
            try:
                with open(asn_path) as asn_file:
                    asn_data = load_asn(asn_file)
            except AssociationNotValidError as e:
                raise OSError("Cannot read ASN file.") from e

            if asn_exptypes is not None:
                asn_data["products"][0]["members"] = [
                    m
                    for m in asn_data["products"][0]["members"]
                    if m["exptype"] in asn_exptypes
                ]

            if asn_n_members is not None:
                asn_data["products"][0]["members"] = asn_data["products"][0]["members"][
                    :asn_n_members
                ]

            # make members easier to access
            self._asn = asn_data
            self._members = self._asn["products"][0]["members"]

            # check that all members have a group_id
            # TODO base this off of the model
            for member in self._members:
                if "group_id" not in member:
                    filename = os.path.join(self._asn_dir, member["expname"])
                    member["group_id"] = _file_to_group_id(filename)
        elif isinstance(init, Iterable):  # assume a list of models
            # make a fake asn from the models
            filenames = set()
            members = []
            for index, model_or_filename in enumerate(init):
                if isinstance(model_or_filename, str):
                    # TODO supporting a list of filenames by opening them as models
                    # has issues, if this is a widely supported mode (vs providing
                    # an association) it might make the most sense to make a fake
                    # association with the filenames at load time.
                    model = datamodels_open(model_or_filename)
                else:
                    model = model_or_filename
                filename = model.meta.filename
                if filename in filenames:
                    raise ValueError(
                        f"Models in library cannot use the same filename: {filename}"
                    )
                self._model_store[index] = model
                # FIXME: output models created during resample (during outlier detection
                # an possibly others) do not have meta.observation which breaks the group_id
                # code
                try:
                    group_id = _model_to_group_id(model)
                except AttributeError:
                    group_id = str(index)
                # FIXME: assign the group id here as it may have been computed above
                # this is necessary for some tweakreg tests that pass in a list of models that
                # don't have group_ids. If this is something we want to support there may
                # be a cleaner way to do this.
                model.meta["group_id"] = group_id
                members.append(
                    {
                        "expname": filename,
                        "exptype": getattr(model.meta, "exptype", "SCIENCE"),
                        "group_id": group_id,
                    }
                )

            if asn_exptypes is not None:
                raise NotImplementedError()

            if asn_n_members is not None:
                raise NotImplementedError()

            # make a fake association
            self._asn = {
                # TODO other asn data?
                "products": [
                    {
                        "members": members,
                    }
                ],
            }
            self._members = self._asn["products"][0]["members"]

        elif isinstance(init, self.__class__):
            # TODO clone/copy?
            raise NotImplementedError()
        else:
            raise NotImplementedError()

        # make sure first model is loaded in memory (as expected by stpipe)
        if asn_n_members == 1:
            # FIXME stpipe also reaches into _models (instead of _model_store)
            self._models = [self._load_member(0)]

    def __del__(self):
        # FIXME when stpipe no longer uses '_models'
        if hasattr(self, "_models"):
            self._models[0].close()

    @property
    def asn(self):
        # return a "read only" association
        def _to_read_only(obj):
            if isinstance(obj, dict):
                return MappingProxyType(obj)
            if isinstance(obj, list):
                return tuple(obj)
            return obj

        return asdf.treeutil.walk_and_modify(self._asn, _to_read_only)

    @property
    def group_names(self):
        names = set()
        for member in self._members:
            names.add(member["group_id"])
        return names

    @property
    def group_indices(self):
        group_dict = {}
        for i, member in enumerate(self._members):
            group_id = member["group_id"]
            if group_id not in group_dict:
                group_dict[group_id] = []
            group_dict[group_id].append(i)
        return group_dict

    def __len__(self):
        return len(self._members)

    def borrow(self, index):
        if not self._open:
            raise ClosedLibraryError("ModelLibrary is not open")

        # if model was already borrowed, raise
        if index in self._ledger:
            raise BorrowError("Attempt to double-borrow model")

        if index in self._model_store:
            model = self._model_store[index]
        else:
            model = self._load_member(index)
            if not self._on_disk:
                # it's ok to keep this in memory since _on_disk is False
                self._model_store[index] = model

        # track the model is "in use"
        self._ledger[index] = model
        return model

    def __getitem__(self, index):
        # FIXME: this is here to allow the library to pass the Sequence
        # check. Removing this will require more extensive stpipe changes
        raise Exception()

    def shelve(self, model, index=None, modify=True):
        if not self._open:
            raise ClosedLibraryError("ModelLibrary is not open")

        if index is None:
            index = self._ledger[model]

        if index not in self._ledger:
            raise BorrowError("Attempt to shelve non-borrowed model")

        if modify:
            self._model_store[index] = model

        del self._ledger[index]

        # TODO should we allow this to change group_id for the member?

    def __iter__(self):
        for i in range(len(self)):
            yield self.borrow(i)

    def _load_member(self, index):
        member = self._members[index]
        filename = os.path.join(self._asn_dir, member["expname"])

        model = datamodels_open(filename, memmap=self._memmap)

        # patch model metadata with asn member info
        # TODO asn.table_name asn.pool_name here?
        for attr in ("group_id", "tweakreg_catalog", "exptype"):
            if attr in member:
                # FIXME model.meta.group_id throws an error
                # setattr(model.meta, attr, member[attr])
                model.meta[attr] = member[attr]
                if attr == "exptype":
                    # FIXME why does tweakreg expect meta.asn.exptype instead of meta.exptype?
                    model.meta["asn"] = {"exptype": member["exptype"]}
        # FIXME tweakreg also expects table_name and pool_name
        model.meta.asn["table_name"] = self.asn_table_name
        model.meta.asn["pool_name"] = self.asn["asn_pool"]
        return model

    def __copy__(self):
        # TODO make copy and deepcopy distinct and not require loading
        # all models into memory
        assert not self._on_disk
        with self:
            model_copies = []
            for i, model in enumerate(self):
                model_copies.append(model.copy())
                self.shelve(model, i, modify=False)
        return self.__class__(model_copies)

    def __deepcopy__(self, memo):
        return self.__copy__()

    def copy(self, memo=None):
        return copy.deepcopy(self, memo=memo)

    # TODO save, required by stpipe
    def save(self, path=None, dir_path=None, save_model_func=None, overwrite=True):
        # FIXME: the signature for this function can lead to many possible outcomes
        # stpipe may call this with save_model_func and path defined
        # skymatch tests call with just dir_path
        # stpipe sometimes provides overwrite=True

        if path is None:

            def path(file_path, index):
                return file_path

        elif not callable(path):

            def path(file_path, index):
                path_head, path_tail = os.path.split(file_path)
                base, ext = os.path.splitext(path_tail)
                if index is not None:
                    base = base + str(index)
                return os.path.join(path_head, base + ext)

        # FIXME: since path is the first argument this means that calling
        # ModelLibrary.save("my_directory") will result in saving all models
        # to the current directory, ignoring "my_directory" this matches
        # what was done for ModelContainer
        dir_path = dir_path if dir_path is not None else os.getcwd()

        # output_suffix = kwargs.pop("output_suffix", None)  # FIXME this was unused

        output_paths = []
        with self:
            for i, model in enumerate(self):
                if len(self) == 1:
                    index = None
                else:
                    index = i
                if save_model_func is None:
                    filename = model.meta.filename
                    output_path, output_filename = os.path.split(path(filename, index))

                    # use dir_path when provided
                    output_path = output_path if dir_path is None else dir_path

                    # create final destination (path + filename)
                    save_path = os.path.join(output_path, output_filename)

                    model.to_asdf(save_path)  # TODO save args?

                    output_paths.append(save_path)
                else:
                    output_paths.append(save_model_func(model, idx=index))

                self.shelve(model, i, modify=False)

        return output_paths

    # TODO crds_observatory, get_crds_parameters, when stpipe uses these...
    def crds_observatory(self):
        return "roman"

    def get_crds_parameters(self):
        raise NotImplementedError()

    def finalize_result(self, step, reference_files_used):
        with self:
            for i, model in enumerate(self):
                step.finalize_result(model, reference_files_used)
                self.shelve(model, i)

    def __enter__(self):
        self._open = True
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._open = False
        if exc_value:
            # if there is already an exception, don't worry about checking the ledger
            # instead allowing the calling code to raise the original error to provide
            # a more useful feedback without any chained ledger exception about
            # un-returned models
            return
        if self._ledger:
            raise BorrowError(
                f"ModelLibrary has {len(self._ledger)} un-returned models"
            )

    def index(self, attribute, copy=False):
        """
        Access a single attribute from all models
        """
        # TODO we could here implement efficient accessors for
        # certain attributes (like `meta.wcs` or `meta.wcs_info.s_region`)
        if copy:
            copy_func = lambda value: value.copy()  # noqa: E731
        else:
            copy_func = lambda value: value  # noqa: E731
        with self:
            for i, model in enumerate(self):
                attr = model[attribute]
                self.shelve(model, i, modify=False)
                yield copy_func(attr)

    def map_function(self, function, modify=False):
        with self:
            for i, model in enumerate(self):
                try:
                    yield function(model)
                finally:
                    # this is in a finally to allow cleanup if the generator is
                    # deleted after it finishes (when it's not fully consumed)
                    self.shelve(model, i, modify)


def _mapping_to_group_id(mapping):
    """
    Combine a number of file metadata values into a ``group_id`` string
    """
    return (
        "roman{program}{observation}{visit}"
        "_{visit_file_group}{visit_file_sequence}{visit_file_activity}"
        "_{exposure}"
    ).format_map(mapping)


def _file_to_group_id(filename):
    """
    Compute a "group_id" without loading the file as a DataModel

    This function will return the meta.group_id stored in the ASDF
    extension (if it exists) or a group_id calculated from the
    FITS headers.
    """
    asdf_yaml = asdf.util.load_yaml(filename)
    if group_id := asdf_yaml["roman"]["meta"].get("group_id"):
        return group_id
    return _mapping_to_group_id(asdf_yaml["roman"]["meta"]["observation"])


def _model_to_group_id(model):
    """
    Compute a "group_id" from a model using the DataModel interface
    """
    if (group_id := getattr(model.meta, "group_id", None)) is not None:
        return group_id
    return _mapping_to_group_id(model.meta.observation)

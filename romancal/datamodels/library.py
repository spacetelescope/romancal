import copy
import os.path
import tempfile
from collections.abc import Iterable, MutableMapping, Sequence
from pathlib import Path
from types import MappingProxyType

import asdf
from roman_datamodels import open as datamodels_open

from .container import ModelContainer


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


class _OnDiskModelStore(MutableMapping):
    def __init__(self, memmap=False, directory=None):
        self._memmap = memmap
        if directory is None:
            # when tem
            self._tempdir = tempfile.TemporaryDirectory(dir="")
            # TODO should I make this a path?
            self._dir = self._tempdir.name
        else:
            self._dir = directory
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
            subdir = os.path.join(self._dir, f"{key}")
            os.makedirs(subdir)
            fn = os.path.join(subdir, model_filename)
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
        self._asn_exptypes = asn_exptypes
        self._asn_n_members = asn_n_members
        self._on_disk = on_disk

        self._open = False
        self._ledger = {}

        # FIXME is there a cleaner way to pass these along to datamodels.open?
        self._memmap = memmap

        if self._on_disk:
            self._model_store = _OnDiskModelStore(memmap, temp_directory)
        else:
            self._model_store = {}

        # TODO path support
        # TODO model list support
        if isinstance(init, (str, Path)):
            self._asn_path = os.path.abspath(
                os.path.expanduser(os.path.expandvars(init))
            )
            self._asn_dir = os.path.dirname(self._asn_path)
            # load association
            # TODO why did ModelContainer make this local?
            from ..associations import AssociationNotValidError, load_asn

            try:
                with open(self._asn_path) as asn_file:
                    asn_data = load_asn(asn_file)
            except AssociationNotValidError as e:
                raise OSError("Cannot read ASN file.") from e

            if self._asn_exptypes is not None:
                asn_data["products"][0]["members"] = [
                    m
                    for m in asn_data["products"][0]["members"]
                    if m["exptype"] in self._asn_exptypes
                ]

            if self._asn_n_members is not None:
                asn_data["products"][0]["members"] = asn_data["products"][0]["members"][
                    : self._asn_n_members
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
                members.append(
                    {
                        "expname": filename,
                        "exptype": getattr(model.meta, "exptype", "SCIENCE"),
                        "group_id": _model_to_group_id(model),
                    }
                )

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

            # _asn_dir?
            # _asn_path?

        elif isinstance(init, self.__class__):
            # TODO clone/copy?
            raise NotImplementedError()
        else:
            raise NotImplementedError()

        # make sure first model is loaded in memory (as expected by stpipe)
        if self._asn_n_members == 1:
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

    # TODO we may want to not expose this as it could go out-of-sync
    # pretty easily with the actual models.
    # @property
    # def members(self):
    #     return self.asn['products'][0]['members']

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

    def __getitem__(self, index):
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

    def __setitem__(self, index, model):
        if index not in self._ledger:
            raise BorrowError("Attempt to return non-borrowed model")

        # un-track this model
        del self._ledger[index]

        # and store it
        self._model_store[index] = model

        # TODO should we allow this to change group_id for the member?

    def discard(self, index, model):
        # TODO it might be worth allowing `discard(model)` by adding
        # an index of {id(model): index} to the ledger to look up the index
        if index not in self._ledger:
            raise BorrowError("Attempt to discard non-borrowed model")

        # un-track this model
        del self._ledger[index]
        # but do not store it

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

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
        # this returns an OPEN model, it's up to calling code to close this
        return model

    def __copy__(self):
        # TODO make copy and deepcopy distinct and not require loading
        # all models into memory
        assert not self._on_disk
        with self:
            model_copies = []
            for i, model in enumerate(self):
                model_copies.append(model.copy())
                self.discard(i, model)
        return self.__class__(model_copies)

    def __deepcopy__(self, memo):
        return self.__copy__()

    def copy(self, memo=None):
        return copy.deepcopy(self, memo=memo)

    # TODO save, required by stpipe
    def save(self, dir_path=None):
        # dir_path: required by SkyMatch tests
        if dir_path is None:
            raise NotImplementedError()
        # save all models
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        with self:
            for i, model in enumerate(self):
                model.save(os.path.join(dir_path, model.meta.filename))
                self.discard(i, model)

    # TODO crds_observatory, get_crds_parameters, when stpipe uses these...

    def _to_container(self):
        # create a temporary directory
        tmpdir = tempfile.TemporaryDirectory(dir="")

        # write out all models (with filenames from member list)
        fns = []
        with self:
            for i, model in enumerate(self):
                fn = os.path.join(tmpdir.name, model.meta.filename)
                model.save(fn)
                fns.append(fn)
                self[i] = model

        # use the new filenames for the container
        # copy over "in-memory" options
        # init with no "models"
        container = ModelContainer(
            fns, save_open=not self._on_disk, return_open=not self._on_disk
        )
        # give the model container a reference to the temporary directory so it's not deleted
        container._tmpdir = tmpdir
        # FIXME container with filenames already skip finalize_result
        return container

    def finalize_result(self, step, reference_files_used):
        with self:
            for i, model in enumerate(self):
                step.finalize_result(model, reference_files_used)
                self[i] = model

    def __enter__(self):
        self._open = True
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._open = False
        # if exc_value:
        #     # if there is already an exception, don't worry about checking the ledger
        #     # instead allowing the calling code to raise the original error to provide
        #     # a more useful feedback without any chained ledger exception about
        #     # un-returned models
        #     return
        # TODO we may want to change this chain to make tracebacks and pytest output
        # easier to read.
        if self._ledger:
            raise BorrowError(
                f"ModelLibrary has {len(self._ledger)} un-returned models"
            ) from exc_value

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
            for i, model in range(len(self)):
                attr = model[attribute]
                self.discard(i, model)
                yield copy_func(attr)


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
    return _mapping_to_group_id(model.meta.observation)

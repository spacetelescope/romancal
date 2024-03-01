import contextlib
import copy
import logging
import os
import os.path as op
import re
from collections import OrderedDict
from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np
from roman_datamodels import datamodels as rdm

from romancal.lib.basic_utils import is_association

from ..associations import AssociationNotValidError, load_asn

__all__ = [
    "ModelContainer",
]

_ONE_MB = 1 << 20
RECOGNIZED_MEMBER_FIELDS = [
    "tweakreg_catalog",
]

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class ModelContainer(Sequence):
    """
    A container for holding DataModels.

    This functions like a list for holding DataModel objects.  It can be
    iterated through like a list and the datamodels within the container can be
    addressed by index. Additionally, the datamodels can be grouped by exposure.

    Parameters
    ----------
    init : path to ASN file, list of either datamodels or path to ASDF files, or `None`
        If `None`, then an empty `ModelContainer` instance is initialized, to which
        datamodels can later be added via the ``insert()``, ``append()``,
        or ``extend()`` method.

    iscopy : bool
        Presume this model is a copy. Members will not be closed
        when the model is closed/garbage-collected.

    memmap : bool
        Open ASDF file binary data using memmap (default: False)

    return_open : bool
        (optional) See notes below on usage.

    save_open : bool
        (optional) See notes below on usage.

    Examples
    --------
    To load a list of ASDF files into a `ModelContainer`:

    .. code-block:: python

        container = ModelContainer(
            [
                "/path/to/file1.asdf",
                "/path/to/file2.asdf",
                ...,
                "/path/to/fileN.asdf"
            ]
        )

    To load a list of open Roman DataModels into a `ModelContainer`:

    .. code-block:: python

        import roman_datamodels.datamodels as rdm
        data_list = [
                "/path/to/file1.asdf",
                "/path/to/file2.asdf",
                ...,
                "/path/to/fileN.asdf"
            ]
        datamodels_list = [rdm.open(x) for x in data_list]
        container = ModelContainer(datamodels_list)

    To load an ASN file into a `ModelContainer`:

    .. code-block:: python

        asn_file = "/path/to/asn_file.json"
        container = ModelContainer(asn_file)

    In any of the cases above, the content of each file in a `ModelContainer` can
    be accessed by iterating over its elements. For example, to print out the filename
    of each file, we can run:

    .. code-block:: python

        for model in container:
            print(model.meta.filename)


    Notes
    -----
    The optional parameters ``save_open`` and ``return_open`` can be
    provided to control how the `DataModel` are used by the
    :py:class:`ModelContainer`. If ``save_open`` is set to `False`, each input
    `DataModel` instance in ``init`` will be written out to disk and
    closed, then only the filename for the `DataModel` will be used to
    initialize the :py:class:`ModelContainer` object.
    Subsequent access of each member will then open the `DataModel` file to
    work with it. If ``return_open`` is also `False`, then the `DataModel`
    will be closed when access to the `DataModel` is completed. The use of
    these parameters can minimize the amount of memory used by this object
    during processing.

    .. warning:: Input files will be updated in-place with new ``meta`` attribute
        values when ASN table's members contain additional attributes.

    """

    def __init__(
        self,
        init=None,
        asn_exptypes=None,
        asn_n_members=None,
        iscopy=False,
        memmap=False,
        # always return an open datamodel
        return_open=True,
        save_open=True,
    ):
        self._models = []
        self._iscopy = iscopy
        self._memmap = memmap
        self._return_open = return_open
        self._save_open = save_open

        self.asn_exptypes = asn_exptypes
        self.asn_n_members = asn_n_members
        self.asn_table = {}
        self.asn_table_name = None
        self.asn_pool_name = None

        try:
            init = Path(init)
        except TypeError:
            if init is None:
                # don't populate container
                pass
            elif isinstance(init, Sequence):
                # only append list items to self._models if all items are either
                # not-null strings (i.e. path to an ASDF file) or instances of DataModel
                is_all_string = all(isinstance(x, str) and len(x) for x in init)
                is_all_roman_datamodels = all(
                    isinstance(x, rdm.DataModel) for x in init
                )
                is_all_path = all(isinstance(x, Path) for x in init)

                if len(init) and (is_all_string or is_all_roman_datamodels):
                    self._models.extend(init)
                elif len(init) and is_all_path:
                    # parse Path object to string
                    self._models.extend([str(x) for x in init])
                else:
                    raise TypeError(
                        "Input must be an ASN file or a list of either strings "
                        "(full path to ASDF files) or Roman datamodels."
                    )
        else:
            if is_association(init):
                self.from_asn(init)
            elif isinstance(init, Path) and init.name != "":
                try:
                    init_from_asn = self.read_asn(init)
                    self.from_asn(init_from_asn, asn_file_path=init)
                except Exception as e:
                    raise TypeError(
                        "Input must be an ASN file or a list of either strings "
                        "(full path to ASDF files) or Roman datamodels."
                    ) from e
            else:
                raise TypeError(
                    "Input must be an ASN file or a list of either strings "
                    "(full path to ASDF files) or Roman datamodels."
                )

    def __len__(self):
        return len(self._models)

    def __getitem__(self, index):
        if isinstance(index, slice):
            start = index.start
            stop = index.stop
            step = index.step
            m = self._models[start:stop:step]
            m = [
                (
                    rdm.open(item, memmap=self._memmap)
                    if (not isinstance(item, rdm.DataModel) and self._return_open)
                    else item
                )
                for item in m
            ]
        else:
            m = self._models[index]
            if not isinstance(m, rdm.DataModel) and self._return_open:
                m = rdm.open(m, memmap=self._memmap)
        return m

    def __setitem__(self, index, model):
        if isinstance(model, rdm.DataModel):
            self._models[index] = model
        else:
            raise ValueError("Only datamodels can be used.")

    def __delitem__(self, index):
        del self._models[index]

    def __iter__(self):
        for model in self._models:
            if not isinstance(model, rdm.DataModel) and self._return_open:
                model = rdm.open(model, memmap=self._memmap)
            yield model

    def insert(self, index, model):
        if isinstance(model, rdm.DataModel):
            self._models.insert(index, model)
        else:
            raise ValueError("Only datamodels can be used.")

    def append(self, model):
        if isinstance(model, rdm.DataModel):
            self._models.append(model)
        else:
            raise ValueError("Only datamodels can be used.")

    def extend(self, input_object):
        if not isinstance(input_object, (Iterable, rdm.DataModel)) or isinstance(
            input_object, str
        ):
            raise ValueError("Not a valid input object.")
        elif all(isinstance(x, rdm.DataModel) for x in input_object):
            self._models.extend(input_object)
        else:
            raise ValueError("Not a valid input object.")

    def pop(self, index=-1):
        self._models.pop(index)

    def copy(self, memo=None):
        """
        Returns a deep copy of the models in this model container.
        """
        return copy.deepcopy(self, memo=memo)

    def close(self):
        """Close all datamodels."""
        if not self._iscopy:
            for model in self._models:
                if isinstance(model, rdm.DataModel):
                    model.close()

    @staticmethod
    def read_asn(filepath):
        """
        Load ASDF files from a Roman association file.

        Parameters
        ----------
        filepath : str
            The path to an association file.
        """
        filepath = op.abspath(op.expanduser(op.expandvars(filepath)))
        try:
            with open(filepath) as asn_file:
                asn_data = load_asn(asn_file)
        except AssociationNotValidError as e:
            raise OSError("Cannot read ASN file.") from e
        return asn_data

    def from_asn(self, asn_data, asn_file_path=None):
        """
        Load ASDF files from a Roman association file.

        Parameters
        ----------
        asn_data : `~roman_datamodels.associations.Association`
            Association dictionary.

        asn_file_path : str
            Filepath of the association, if known.
        """
        # match the asn_exptypes to the exptype in the association and retain
        # only those file that match, as a list, if asn_exptypes is set to none
        # grab all the files
        if self.asn_exptypes:
            infiles = []
            logger.debug(
                f"Filtering datasets based on allowed exptypes {self.asn_exptypes}:"
            )
            for member in asn_data["products"][0]["members"]:
                if any(
                    x
                    for x in self.asn_exptypes
                    if re.match(member["exptype"], x, re.IGNORECASE)
                ):
                    infiles.append(member)
                    logger.debug(f'Files accepted for processing {member["expname"]}:')
        else:
            infiles = list(asn_data["products"][0]["members"])

        asn_dir = op.dirname(asn_file_path) if asn_file_path else ""
        # Only handle the specified number of members.
        sublist = infiles[: self.asn_n_members] if self.asn_n_members else infiles
        try:
            for member in sublist:
                filepath = op.join(asn_dir, member["expname"])
                update_model = any(attr in member for attr in RECOGNIZED_MEMBER_FIELDS)
                if update_model or self._save_open:
                    m = rdm.open(filepath, memmap=self._memmap)
                    m.meta["asn"] = {"exptype": member["exptype"]}
                    for attr, val in member.items():
                        if attr in RECOGNIZED_MEMBER_FIELDS:
                            if attr == "tweakreg_catalog":
                                val = op.join(asn_dir, val) if val.strip() else None
                            m.meta[attr] = val

                    if not self._save_open:
                        m.save(filepath)
                        m.close()
                else:
                    m = filepath

                self._models.append(m)

        except OSError:
            self.close()
            raise

        # Pull the whole association table into asn_table
        self.merge_tree(self.asn_table, asn_data)

        if asn_file_path is not None:
            self.asn_table_name = op.basename(asn_file_path)
            self.asn_pool_name = asn_data["asn_pool"]
            for model in self:
                with contextlib.suppress(AttributeError):
                    model.meta.asn["table_name"] = self.asn_table_name
                    model.meta.asn["pool_name"] = self.asn_pool_name

    def save(self, path=None, dir_path=None, save_model_func=None, **kwargs):
        """
        Write out models in container to ASDF.

        Parameters
        ----------
        path : str or func or None
            - If None, the `meta.filename` is used for each model.
            - If a string, the string is used as a root and an index is
              appended.
            - If a function, the function takes the two arguments:
              the value of model.meta.filename and the
              `idx` index, returning constructed file name.

        dir_path : str
            Directory to write out files.  Defaults to current working dir.
            If directory does not exist, it creates it.  Filenames are pulled
            from `.meta.filename` of each datamodel in the container.

        save_model_func: func or None
            Alternate function to save each model instead of
            the models `save` method. Takes one argument, the model,
            and keyword argument `idx` for an index.

        Note
        ----
        Additional parameters provided via `**kwargs` are passed on to
        `roman_datamodels.datamodels.DataModel.to_asdf`

        Returns
        -------
        output_paths: [str[, ...]]
            List of output file paths of where the models were saved.
        """
        output_paths = []
        if path is None:

            def path(filename, idx=None):
                return filename

        elif not callable(path):
            path = make_file_with_index

        # use current path if dir_path is not provided
        dir_path = dir_path if dir_path is not None else os.getcwd()
        # output filename suffix
        output_suffix = kwargs.pop("output_suffix", None)
        for idx, model in enumerate(self._models):
            if len(self) <= 1:
                idx = None
            if save_model_func is None:
                filename = model.meta.filename
                output_path, output_filename = op.split(path(filename, idx=idx))

                # use dir_path when provided
                output_path = output_path if dir_path is None else dir_path

                # handle optional modifications to filename
                base, ext = op.splitext(output_filename)
                if output_suffix is not None:
                    # add suffix to filename
                    base = "".join([base, output_suffix])
                output_filename = "".join([base, ext])

                # create final destination (path + filename)
                save_path = op.join(output_path, output_filename)

                if ext == ".asdf":
                    output_paths.append(save_path)
                    model.to_asdf(save_path, **kwargs)
                else:
                    raise ValueError(f"Unknown filetype {ext}")
            else:
                output_paths.append(save_model_func(model, idx=idx))

        return output_paths

    @property
    def models_grouped(self):
        """
        Returns a list of a list of datamodels grouped by exposure.
        Assign an ID grouping by exposure.

        Data from different detectors of the same exposure will have the
        same group id, which allows grouping by exposure.  The following
        metadata is used for grouping:

        meta.observation.program
        meta.observation.observation
        meta.observation.visit
        meta.observation.visit_file_group
        meta.observation.visit_file_sequence
        meta.observation.visit_file_activity
        meta.observation.exposure
        """
        unique_exposure_parameters = [
            "program",
            "observation",
            "visit",
            "visit_file_group",
            "visit_file_sequence",
            "visit_file_activity",
            "exposure",
        ]

        group_dict = OrderedDict()
        for i, model in enumerate(self._models):
            model = model if isinstance(model, rdm.DataModel) else rdm.open(model)

            if not self._save_open:
                model = rdm.open(model, memmap=self._memmap)

            params = [
                str(getattr(model.meta.observation, param))
                for param in unique_exposure_parameters
            ]
            try:
                group_id = "roman" + "_".join(
                    ["".join(params[:3]), "".join(params[3:6]), params[6]]
                )
                model.meta["group_id"] = group_id
            except TypeError:
                model.meta["group_id"] = f"exposure{i + 1:04d}"

            group_id = model.meta.group_id
            if not self._save_open and not self._return_open:
                model.close()
                model = self._models[i]

            if group_id in group_dict:
                group_dict[group_id].append(model)
            else:
                group_dict[group_id] = [model]

        return group_dict.values()

    def merge_tree(self, a, b):
        """
        Merge elements from tree ``b`` into tree ``a``.
        """

        def recurse(a, b):
            if isinstance(b, dict):
                if not isinstance(a, dict):
                    return copy.deepcopy(b)
                for key, val in b.items():
                    a[key] = recurse(a.get(key), val)
                return a
            return copy.deepcopy(b)

        recurse(a, b)
        return a

    @property
    def crds_observatory(self):
        """
        Get the CRDS observatory for this container.  Used when selecting
        step/pipeline parameter files when the container is a pipeline input.

        Returns
        -------
        str
        """
        return "roman"

    def get_crds_parameters(self):
        """
        Get parameters used by CRDS to select references for this model.

        Returns
        -------
        dict
        """
        crds_header = {}
        if len(self._models):
            model = self._models[0]
            model = model if isinstance(model, rdm.DataModel) else rdm.open(model)
            crds_header |= model.get_crds_parameters()

        return crds_header

    def set_buffer(self, buffer_size, overlap=None):
        """Set buffer size for scrolling section-by-section access.

        Parameters
        ----------
        buffer_size : float, None
            Define size of buffer in MB for each section.
            If `None`, a default buffer size of 1MB will be used.

        overlap : int, optional
            Define the number of rows of overlaps between sections.
            If `None`, no overlap will be used.
        """
        self.overlap = 0 if overlap is None else overlap
        self.grow = 0

        with rdm.open(self._models[0]) as model:
            imrows, imcols = model.data.shape
            data_item_size = model.data.itemsize
            data_item_type = model.data.dtype
            del model

        min_buffer_size = imcols * data_item_size

        self.buffer_size = (
            min_buffer_size if buffer_size is None else (buffer_size * _ONE_MB)
        )

        section_nrows = min(imrows, int(self.buffer_size // min_buffer_size))

        if section_nrows == 0:
            self.buffer_size = min_buffer_size
            logger.warning(
                "WARNING: Buffer size is too small to hold a single row."
                f"Increasing buffer size to {self.buffer_size / _ONE_MB}MB"
            )
            section_nrows = 1

        nbr = section_nrows - self.overlap
        nsec = (imrows - self.overlap) // nbr
        if (imrows - self.overlap) % nbr > 0:
            nsec += 1

        self.n_sections = nsec
        self.nbr = nbr
        self.section_nrows = section_nrows
        self.imrows = imrows
        self.imcols = imcols
        self.imtype = data_item_type

    def get_sections(self):
        """Iterator to return the sections from all members of the container."""

        for k in range(self.n_sections):
            e1 = k * self.nbr
            e2 = e1 + self.section_nrows

            if k == self.n_sections - 1:  # last section
                e2 = min(e2, self.imrows)
                e1 = min(e1, e2 - self.overlap - 1)

            data_list = np.empty(
                (len(self._models), e2 - e1, self.imcols), dtype=self.imtype
            )
            wht_list = np.empty(
                (len(self._models), e2 - e1, self.imcols), dtype=self.imtype
            )
            for i, model in enumerate(self._models):
                model = rdm.open(model, memmap=self._memmap)

                data_list[i, :, :] = model.data[e1:e2].copy()
                wht_list[i, :, :] = model.weight[e1:e2].copy()
                del model

            yield (data_list, wht_list, (e1, e2))


def make_file_with_index(file_path, idx):
    """Append an index to a filename

    Parameters
    ----------
    file_path: str
        The file to append the index to.
    idx: int
        An index to append


    Returns
    -------
    file_path: str
        Path with index appended
    """
    # Decompose path
    path_head, path_tail = op.split(file_path)
    base, ext = op.splitext(path_tail)
    if idx is not None:
        base = base + str(idx)
    return op.join(path_head, base + ext)

import os
import os.path
import warnings
from collections import OrderedDict
from collections.abc import Iterable

from roman_datamodels.datamodels import DataModel

import asdf
import packaging.version

# .dev is included in the version comparison to allow for correct version
# comparisons with development versions of asdf 3.0
if packaging.version.Version(asdf.__version__) < packaging.version.Version(
    "3.dev"
):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=asdf.exceptions.AsdfDeprecationWarning,
            message=r"AsdfInFits has been deprecated.*",
        )
        from asdf.fits_embed import AsdfInFits
else:
    AsdfInFits = None

__all__ = [
    "ModelContainer",
]


class ModelContainer(Iterable):
    """
    A container for holding DataModels.

    This functions like a list for holding DataModel objects.  It can be
    iterated through like a list and the datamodels within the container can be
    addressed by index. Additionally, the datamodels can be grouped by exposure.

    Parameters
    ----------
    init : file path, list of DataModels or path to ASDF files, or None

        - file path: initialize from an association table

        - list: a list of either datamodels or full path to ASDF files (as string)

        - None: initializes an empty `ModelContainer` instance, to which
          DataModels can be added via the ``append()`` method.

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
    >>> container = ModelContainer(['file1.asdf', 'file2.asdf', ..., 'fileN.asdf'])
    >>> for model in container:
    ...     print(model.meta.filename)


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
        iscopy=False,
        memmap=False,
        return_open=True,
        save_open=True,
    ):
        self._models = []
        self._iscopy = iscopy
        self._memmap = memmap
        self._return_open = return_open
        self._save_open = save_open

        if init is None:
            # don't populate container
            pass
        elif isinstance(init, Iterable):
            # only append list items to self._models if all items are either
            # strings (i.e. path to an ASDF file) or instances of DataModel
            is_all_string = all(isinstance(x, str) for x in init)
            is_all_roman_datamodels = all(
                isinstance(x, DataModel) for x in init
            )

            if is_all_string or is_all_roman_datamodels:
                self._models.extend(init)
            else:
                raise TypeError(
                    "Input must be a list of strings (full path to ASDF files) or Roman datamodels."
                )
        else:
            raise TypeError(
                "Input must be a list of either strings (full path to ASDF files) or Roman datamodels."
            )

    def __len__(self):
        return len(self._models)

    def __getitem__(self, index):
        m = self._models[index]
        if not isinstance(m, DataModel) and self._return_open:
            m = open(m, memmap=self._memmap)
        return m

    def __setitem__(self, index, model):
        self._models[index] = model

    def __iter__(self):
        for model in self._models:
            if not isinstance(model, DataModel) and self._return_open:
                model = open(model, memmap=self._memmap)
            yield model

    def save(self, dir_path=None, *args, **kwargs):
        # use current path if dir_path is not provided
        dir_path = dir_path if dir_path is not None else os.getcwd()
        # output filename suffix
        output_suffix = "cal_twkreg"
        for model in self._models:
            filename = model.meta.filename
            base, ext = os.path.splitext(filename)
            base = base.replace("cal", output_suffix)
            output_filename = "".join([base, ext])
            output_path = os.path.join(dir_path, output_filename)
            # TODO: Support gzip-compressed fits
            if ext == ".asdf":
                model.to_asdf(output_path, *args, **kwargs)
            else:
                raise ValueError(f"unknown filetype {ext}")

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
            model = model if isinstance(model, DataModel) else open(model)

            if not self._save_open:
                model = open(model, memmap=self._memmap)

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

    @property
    def to_association(self):
        pass

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

    @property
    def get_crds_parameters(self):
        """
        Get parameters used by CRDS to select references for this model.

        Returns
        -------
        dict
        """
        return {
            key: val
            for key, val in self.to_flat_dict(include_arrays=False).items()
            if isinstance(val, (str, int, float, complex, bool))
        }

import asdf
from roman_datamodels import open as datamodels_open
from stpipe.library import AbstractModelLibrary, NoGroupID

from romancal.associations import AssociationNotValidError, load_asn

__all__ = ["ModelLibrary"]


class ModelLibrary(AbstractModelLibrary):
    @property
    def crds_observatory(self):
        return "roman"

    def _model_to_filename(self, model):
        model_filename = model.meta.filename
        if model_filename is None:
            model_filename = "model.asdf"
        return model_filename

    def _datamodels_open(self, filename, **kwargs):
        return datamodels_open(filename, **kwargs)

    @classmethod
    def _load_asn(cls, asn_path):
        try:
            with open(asn_path) as asn_file:
                asn_data = load_asn(asn_file)
        except AssociationNotValidError as e:
            raise OSError("Cannot read ASN file.") from e
        return asn_data

    def _filename_to_group_id(self, filename):
        """
        Compute a "group_id" without loading the file as a DataModel

        This function will return the meta.group_id stored in the ASDF
        extension (if it exists) or a group_id calculated from the
        ASDF headers.
        """
        meta = asdf.util.load_yaml(filename)["roman"]["meta"]
        if group_id := meta.get("group_id"):
            return group_id
        if "observation" in meta:
            return _mapping_to_group_id(meta["observation"])
        raise NoGroupID(f"{filename} missing group_id")

    def _model_to_group_id(self, model):
        """
        Compute a "group_id" from a model using the DataModel interface
        """
        if (group_id := getattr(model.meta, "group_id", None)) is not None:
            return group_id
        if hasattr(model.meta, "observation"):
            return _mapping_to_group_id(model.meta.observation)
        raise NoGroupID(f"{model} missing group_id")

    def _assign_member_to_model(self, model, member):
        # roman_datamodels doesn't allow assignment of meta.group_id
        # (since it's not in the schema). To work around this use
        # __setitem__ calls here instead of setattr
        for attr in ("group_id", "tweakreg_catalog", "exptype"):
            if attr in member:
                model.meta[attr] = member[attr]
        if not hasattr(model.meta, "asn"):
            model.meta["asn"] = {}

        if "table_name" in self.asn:
            model.meta.asn["table_name"] = self.asn["table_name"]
        if "asn_pool" in self.asn:
            model.meta.asn["pool_name"] = self.asn["asn_pool"]


def _mapping_to_group_id(mapping):
    """
    Combine a number of file metadata values into a ``group_id`` string
    """
    return (
        "roman{program}{observation}{visit}"
        "_{visit_file_group}{visit_file_sequence}{visit_file_activity}"
        "_{exposure}"
    ).format_map(mapping)

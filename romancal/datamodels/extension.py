from pathlib import Path

from asdf.extension import AsdfExtension
from asdf import util


SCHEMAS_ROOT = (Path(__file__).parent/"schemas").resolve()
URI_PREFIX = "http://stsci.edu/schemas/roman_datamodel/"


class RomanDataModelExtension(AsdfExtension):
    """
    Extension that maps datamodel schema URIs to their corresponding
    locations on the disk.  This allows the asdf package to locate
    the schema content when we request validation against one of these
    URIs.
    """

    # The AsdfExtension interface requires types and tag_mapping attributes,
    # but we won't actually be using them here.
    types = []
    tag_mapping = []

    url_mapping = [
        (URI_PREFIX, util.filepath_to_url(str(SCHEMAS_ROOT)) + "/{url_suffix}.yaml"),
    ]

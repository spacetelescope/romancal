import yaml
import asdf

from romancal.datamodels.extension import SCHEMAS_ROOT


def test_schema_uri_mapping():
    """
    Confirm that each schema's URI is mapped by asdf
    to the correct schema content.
    """
    for schema_path in SCHEMAS_ROOT.glob("**/*.yaml"):
        schema = yaml.safe_load(schema_path.read_bytes())
        schema_uri = schema["id"]

        assert asdf.schema.load_schema(schema_uri) == schema

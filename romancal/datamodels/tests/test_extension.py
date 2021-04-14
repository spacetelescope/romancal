import yaml
import asdf

from rad.integration import get_resource_mappings


def test_schema_uri_mapping():
    """
    Confirm that each schema's URI is mapped by asdf
    to the correct schema content.
    """
    schema_mappings = get_resource_mappings()[0]

    schemas = [schema_mappings[item] for item in schema_mappings]
    for schema_text in schemas:
        schema = yaml.safe_load(schema_text)
        schema_uri = schema["id"]

        assert asdf.schema.load_schema(schema_uri) == schema

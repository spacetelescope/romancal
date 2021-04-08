import yaml
import asdf

from rad.extension import get_resource_mappings


def test_schema_uri_mapping():
    """
    Confirm that each schema's URI is mapped by asdf
    to the correct schema content.
    """
    schema_mappings = get_resource_mappings()[0]
    schema_paths = [item + '.yaml' for item in schema_mappings]

    for schema_path in schema_paths:
    	print(schema_path)
        schema = yaml.safe_load(schema_path)
        schema_uri = schema["id"]

        assert asdf.schema.load_schema(schema_uri) == schema

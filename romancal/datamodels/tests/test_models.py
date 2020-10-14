import asdf

from romancal import datamodels


def test_model_schemas():
    model_classes = [
        m for m in datamodels.__dict__.values()
        if isinstance(m, type) and issubclass(m, datamodels.RomanDataModel)
    ]
    for model_class in model_classes:
        try:
            asdf.schema.load_schema(model_class.schema_url)
        except Exception:
            assert False, f"Unable to load schema URI '{model_class.schema_url}' for model class {model_class.__name__}"

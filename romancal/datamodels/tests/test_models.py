import asdf

from romancal.datamodels import RomanDataModel


def test_model_schemas():
    def iter_subclasses(model_class):
        yield model_class
        for sub_class in model_class.__subclasses__():
            yield from iter_subclasses(sub_class)

    for model_class in iter_subclasses(RomanDataModel):
        try:
            asdf.schema.load_schema(model_class.schema_url)
        except Exception:
            assert False, f"Unable to load schema URI '{model_class.schema_url}' for model class {model_class.__name__}"

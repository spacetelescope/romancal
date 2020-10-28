import pytest

import asdf
import numpy as np

from romancal.datamodels import RomanDataModel
import romancal.datamodels as DataModels
from stdatamodels.validate import ValidationWarning


def test_model_schemas(tmp_path):
    def iter_subclasses(model_class):
        yield model_class
        for sub_class in model_class.__subclasses__():
            yield from iter_subclasses(sub_class)

    for model_class in iter_subclasses(RomanDataModel):
        try:
            asdf.schema.load_schema(model_class.schema_url)
        except Exception:
            assert False, f"Unable to load schema URI '{model_class.schema_url}' for model " \
                          f"class {model_class.__name__}"

def test_referencefile_model_schemas(tmp_path):
    print("XXX test_referencefile_model_schemas")
    print("XXX tmp_path = " + str(tmp_path))

    def iter_subclasses(model_class):
        yield model_class
        for sub_class in model_class.__subclasses__():
            yield from iter_subclasses(sub_class)

    file_path = tmp_path / "test.asdf"

    for model_class in iter_subclasses(DataModels.referencefile.ReferenceFileModel):
        print("XXX model_class = "+str(model_class))
        with asdf.AsdfFile() as af:
            af["meta"] = {"author":"Space Telescope Science Institute",
                          "date": "Today",
                          "description":"Test file",
                          "instrument":{
                              "detector":"WFI01",
                              "filter":"F158",
                              "name":"WFI"
                          },
                          "pedigree":"Test",
                          "reftype":"REFERENCEFILE",
                          "telescope":"ROMAN",
                          "useafter":"Today"
            }


        # Test referencefile model
        if (model_class == DataModels.referencefile.ReferenceFileModel):
            print("XXX Referencefile model found")
            af.write_to(file_path)

            assert DataModels.open(file_path)

            # with DataModels.open(file_path) as model:
            #     assert isinstance(model, DataModels.referencefile.ReferenceFileModel)

        # Test Flat file model
        elif (model_class == DataModels.flat.FlatModel):
            print("XXX Flat model found")
            af["meta"]["reftype"] = "FLAT"
            af["data"] = np.zeros((4096, 4096))
            af["dq"] = np.zeros((4096, 4096))
            af["derr"] = np.zeros((4096, 4096))
            af.write_to(file_path)

            with DataModels.open(file_path) as model:
                with pytest.warns(None):
                    print("XXX model.validate = "+str(model.validate()))
                    model.validate()
                    assert not model.validate()
                assert isinstance(model, DataModels.flat.FlatModel)

            af["meta"]["telescope"] = "NOTROMAN"
            af.write_to(file_path)
            with DataModels.open(file_path) as model:
                with pytest.warns(ValidationWarning):
                    print("XXX Bad model.validate = " + str(model.validate()))
                    print("XXX model['meta']['telescope'] = " + str(model["meta"]["telescope"]))
                    model.validate()
                    assert model["meta"]["telescope"] == "NOTROMAN"
                    assert model.validate() != None

        else:
            raise Exception("XXX Unknown referencefile model found:" + str(model_class))


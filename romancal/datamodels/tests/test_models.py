import pytest

import asdf
import numpy as np
from astropy.time import Time
from copy import deepcopy

from romancal.datamodels import RomanDataModel
import romancal.datamodels as DataModels
from stdatamodels.validate import ValidationWarning

# Helper class to iterate over model subclasses
def iter_subclasses(model_class):
    yield model_class
    for sub_class in model_class.__subclasses__():
        yield from iter_subclasses(sub_class)

def test_model_schemas(tmp_path):
    for model_class in iter_subclasses(RomanDataModel):
        try:
            asdf.schema.load_schema(model_class.schema_url)
        except Exception:
            assert False, f"Unable to load schema URI '{model_class.schema_url}' for model " \
                          f"class {model_class.__name__}"



# Helper function to test required keywords
# meta = dictionary of required key/value pairs to ensure each are required
def confirm_required_keywords(tmp_path, meta):
    file_path = tmp_path / "test.asdf"

    # Returns the "paths" for all the keys in a dictionary (e.g. keys within keys.)
    # This return is in the form of a list of lists, allowing the paths to be combined in whatever
    # fashion is appropriate
    def getkeypaths(meta_dict, path=None):
        if path is None:
            path = []
        for keyname, value in meta_dict.items():
            newpath = path + [keyname]
            if isinstance(value, dict):
                for subkeyname in getkeypaths(value, newpath):
                    yield subkeyname
            else:
                yield newpath

    # Deletes the key/value pair for a "pathed key" (e.g. a key within keys.)
    # Returns the dictionary minus that key/value pair.
    def delkeypath(meta_dict, keypath):
        keylist = keypath.split('.')
        sub_dict = meta_dict
        for subkey in keylist[:-1]:
            sub_dict = sub_dict[subkey]
        del sub_dict[keylist[-1]]
        return sub_dict

    # Cycle through each key/value pair in meta dictionary
    for keypathlist in getkeypaths(meta):
        keypathname = '.'.join(keypathlist)

        # Create temporary copy of meta dictionary and remove one key/value pair
        badmeta = deepcopy(meta)
        delkeypath(badmeta,keypathname)

        # Save asdf file with badmeta dictionary as metadata tree
        with asdf.AsdfFile(badmeta) as af:
            af.write_to(file_path)

        # Test for warning that required keyword is missing
        with DataModels.open(file_path) as model:
            with pytest.warns(ValidationWarning) as record:
                model.validate()
        assert (len(record) > 0)

# Testing core schema
def test_core_schemas(tmp_path):

    # Define meta dictionary of required key / value pairs
    meta = {"meta":
                {"instrument": {
                     "detector": "WFI01",
                     "optical_element": "F158",
                     "name": "WFI"
                 },
                 "telescope": "ROMAN"
                 }
            }

    # Test schema for required key/value pairs
    confirm_required_keywords(tmp_path, meta)

# Testing all reference file schemas
def test_referencefile_model_schemas(tmp_path):

    # Define base meta dictionary of required key / value pairs for reference files
    # NOTE 1: date is commented out because it is not properly failing its test
    #         Ticket: https://github.com/spacetelescope/stdatamodels/issues/23
    # NOTE 2: core.schema key/value pairs commented out due to not properly failing their tests
    #         Ticket: https://github.com/spacetelescope/stdatamodels/issues/7

    meta = {"meta":
                {
                    "author": "Space Telescope Science Institute",
                    "description": "Test file",
                    "pedigree": "Test",
                    "reftype": "REFERENCEFILE",
                    "useafter": Time('1999-01-01T00:00:00.123456789',format='isot', scale='utc'),
#                    "date": Time('1999-01-01T00:00:00.123456789',format='isot', scale='utc'),
#                     # Key/value pairs from core.schema
#                     "instrument": {
#                         "detector": "WFI01",
#                         "optical_element": "F158",
#                         "name": "WFI"
#                     },
#                     "telescope": "ROMAN",
                }
    }

    # Set temporary asdf file
    file_path = tmp_path / "test.asdf"

    # Cycle through reference file models to test each
    for model_class in iter_subclasses(DataModels.referencefile.ReferenceFileModel):
        # Test base referencefile model
        if (model_class == DataModels.referencefile.ReferenceFileModel):
            # Test required key/values pairs in meta dictionary
            confirm_required_keywords(tmp_path, meta)

            # Testing referencefile asdf file
            with asdf.AsdfFile(meta) as af:
                af.write_to(file_path)

                # Load base referencefile schema
                schema = asdf.schema.load_schema(
                    "http://stsci.edu/schemas/roman_datamodel/reference_files/referencefile.schema")

                # Ensure that base referencefile schema contains core schema
                found = False
                for schemaref in schema.get("allOf", [{}]):
                    if schemaref.get("$ref","") == "../core.schema":
                        found = True
                        break
                assert found

                # Test that asdf file opens properly
                assert DataModels.open(file_path)

                # Confirm that asdf file is opened as base referencefile model
                with DataModels.open(file_path) as model:
                    assert isinstance(model, DataModels.referencefile.ReferenceFileModel)

        # Test Flat file model
        elif (model_class == DataModels.flat.FlatModel):

            # Testing flat file asdf file
            with asdf.AsdfFile(meta) as af:
                # Add required flat file elements
                af["meta"]["reftype"] = "FLAT"
                af["data"] = np.zeros((4096, 4096))
                af["dq"] = np.zeros((4096, 4096))
                af["derr"] = np.zeros((4096, 4096))
                af.write_to(file_path)

                # Test that asdf file opens properly
                with DataModels.open(file_path) as model:
                    with pytest.warns(None) as record:
                        model.validate()

                    # Confirm that asdf file is opened as flat file model
                    assert isinstance(model, DataModels.flat.FlatModel)
                assert (len(record) == 0)

                # Test for invalid flat file
                af["meta"]["telescope"] = "NOTROMAN"
                af.write_to(file_path)

                with DataModels.open(file_path) as model:
                    model.validate()
                    with pytest.warns(ValidationWarning) as record:
                        model.validate()
                    assert model["meta"]["telescope"] == "NOTROMAN"
                assert (len(record) > 0)

        # Raise error if unknown reference file encountered
        else:
            raise Exception("XXX Unknown referencefile model found:" + str(model_class))

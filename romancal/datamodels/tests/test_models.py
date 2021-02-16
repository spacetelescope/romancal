import pytest

import asdf
import numpy as np
from astropy.time import Time

from romancal import datamodels
from stdatamodels.validate import ValidationWarning

# Helper class to iterate over model subclasses
def iter_subclasses(model_class, include_base_model=True):
    if include_base_model:
        yield model_class
    for sub_class in model_class.__subclasses__():
        yield from iter_subclasses(sub_class)

def test_model_schemas(tmp_path):
    for model_class in iter_subclasses(datamodels.RomanDataModel):
        asdf.schema.load_schema(model_class.schema_url)

# Testing core schema
def test_core_schema(tmp_path):
    # Set temporary asdf file
    file_path = tmp_path / "test.asdf"

    # Define meta dictionary of required key / value pairs
    meta = {"meta":
                {"instrument": {
                     "detector": "WFI01",
                     "optical_element": "F158",
                     "name": "WFI"
                 },
                 "telescope": "ROMAN",
                 "model_type": "RomanDataModel",
                 "date": "2020-12-20"
                 }
            }

    # Testing bad core asdf file
    with asdf.AsdfFile(meta) as af:
        # Test for invalid entry
        af["meta"]["telescope"] = "NOTROMAN"
        af.write_to(file_path)

        with datamodels.open(file_path) as model:
            with pytest.warns(ValidationWarning):
                model.validate()
            assert model["meta"]["telescope"] == "NOTROMAN"

# Helper function to confirm that a reference exists within a schema
def assert_referenced_schema(schema_uri, ref_uri):
    # Hide function from pytest traceback
    __tracebackhide__ = True

    # Load schema
    schema = asdf.schema.load_schema(schema_uri)

    # Search for reference in schema
    for subschema in schema.get("allOf", []):
        if asdf.generic_io.resolve_uri(schema_uri, subschema.get("$ref")) == ref_uri:
            # Reference found
            return
    # Reference not found
    assert False, f"Schema '{schema_uri}' does not reference '{ref_uri}' in a top-level allOf " \
                  f"combiner"


# Define base meta dictionary of required key / value pairs for reference files
#
# NOTE 2: core.schema key/value pairs commented out due to not properly failing their tests
#         Ticket: https://github.com/spacetelescope/stdatamodels/issues/7
REFERENCEFILE_SCHEMA_DICT = {
    "meta": {
                    "author": "Space Telescope Science Institute",
                    "description": "Test file",
                    "pedigree": "Test",
                    "reftype": "BASE",
                    "useafter": Time('1999-01-01T00:00:00.123456789',format='isot', scale='utc'),
                    "date": Time('1999-01-01T00:00:00.123456789',format='isot', scale='utc'),
#                     # Key/value pairs from core.schema
#                     "instrument": {
#                         "detector": "WFI01",
#                         "optical_element": "F158",
#                         "name": "WFI"
#                     },
#                     "telescope": "ROMAN",
#                     "model_type": "ReferenceFileModel",
        }
    }

# Testing all reference file schemas
def test_reference_file_model_base(tmp_path):
    # Set temporary asdf file
    file_path = tmp_path / "test.asdf"

    # Get reference dictionary base
    meta = REFERENCEFILE_SCHEMA_DICT

    # Testing referencefile asdf file
    with asdf.AsdfFile(meta) as af:
        af.write_to(file_path)

        # Ensure that base referencefile schema contains core schema
        assert_referenced_schema(datamodels.reference_files.referencefile.ReferenceFileModel.schema_url,
                                 datamodels.core.RomanDataModel.schema_url)

        # Test that asdf file opens properly
        assert datamodels.open(file_path)

        # Confirm that asdf file is opened as base referencefile model
        with datamodels.open(file_path) as model:
            assert isinstance(model, datamodels.reference_files.referencefile.ReferenceFileModel)


# Common tests for all ReferenceFileModel subclasses
@pytest.mark.parametrize("model_class", list(iter_subclasses(
    datamodels.reference_files.referencefile.ReferenceFileModel, include_base_model=False)))
def test_reference_file_model(tmp_path, model_class):
    # Ensure that specific reference file schema contains the base module schema
    assert_referenced_schema(model_class.schema_url,
                             datamodels.reference_files.referencefile.ReferenceFileModel.schema_url)

# FlatModel tests
def test_flat_model(tmp_path):
    # Set temporary asdf file
    file_path = tmp_path / "test.asdf"

    # Get reference dictionary base
    meta = REFERENCEFILE_SCHEMA_DICT

    # Testing flat file asdf file
    with asdf.AsdfFile(meta) as af:
        # Add required flat file elements
        af["meta"]["reftype"] = "FLAT"
        af["meta"]["model_type"] = "FlatModel"
        af["data"] = np.zeros((4096, 4096))
        af["dq"] = np.zeros((4096, 4096))
        af["derr"] = np.zeros((4096, 4096))
        af.write_to(file_path)

        # Test that asdf file opens properly
        with datamodels.open(file_path) as model:
            with pytest.warns(None):
                model.validate()

            # Confirm that asdf file is opened as flat file model
            assert isinstance(model, datamodels.reference_files.flat.FlatModel)


def test_meta_date_management(tmp_path):
    model = datamodels.RomanDataModel({
        "meta": {
            "date": Time("2000-01-01T00:00:00.000"),
            "instrument": {"name": "WFI", "detector": "WFI01", "optical_element": "F062"},
            "telescope": "ROMAN",
        }
    })
    assert model.meta.date == Time("2000-01-01T00:00:00.000")
    model.save(str(tmp_path/"test.asdf"))
    assert abs((Time.now() - model.meta.date).value) < 1.0

    model = datamodels.RomanDataModel()
    assert abs((Time.now() - model.meta.date).value) < 1.0

import pytest

from jsonschema import ValidationError
import asdf
import numpy as np
from astropy.time import Time

from roman_datamodels import datamodels
from roman_datamodels.tests import util
from roman_datamodels.validate import ValidationWarning
import roman_datamodels.rconverters as converters
from roman_datamodels.extensions import DATAMODEL_EXTENSIONS
from roman_datamodels.util import get_schema_uri_from_converter
from rad.integration import get_resource_mappings

EXPECTED_COMMON_REFERENCE = \
    {'$ref': 'http://stsci.edu/schemas/datamodels/roman/reference_files/ref_common-1.0.0'}

# Helper class to iterate over model subclasses
def iter_subclasses(model_class, include_base_model=True):
    if include_base_model:
        yield model_class
    for sub_class in model_class.__subclasses__():
        yield from iter_subclasses(sub_class)

def test_model_schemas():
    dmodels = datamodels.model_registry.keys()
    for model in dmodels:
        schema_uri = next(t for t in DATAMODEL_EXTENSIONS[0].tags 
            if t._tag_uri==model._tag)._schema_uri
        asdf.schema.load_schema(schema_uri)

# Testing core schema
def test_core_schema(tmp_path):
    # Set temporary asdf file
    file_path = tmp_path / "test.asdf"


    wfi_image = util.mk_level2_image(arrays=(10, 10))
    with asdf.AsdfFile() as af:
        af.tree = {'roman': wfi_image}
        with pytest.raises(ValidationError):
            af.tree['roman'].meta.telescope = 'NOTROMAN'
        af.tree['roman'].meta['telescope'] = 'NOTROMAN'
        with pytest.raises(ValidationError):
            af.write_to(file_path)
        af.tree['roman'].meta.telescope = 'ROMAN'
        af.write_to(file_path)
    # Now mangle the file 
    with open(file_path, 'rb') as fp:
        fcontents = fp.read()
    romanloc = fcontents.find(bytes('ROMAN', 'utf-8'))
    newcontents = fcontents[:romanloc] + bytes('X', 'utf-8') + fcontents[romanloc+1:]
    with open(file_path, 'wb') as fp:
        fp.write(newcontents)
    with pytest.raises(ValidationError): 
        with datamodels.open(file_path) as model:
            pass
    asdf.get_config().validate_on_read = False
    with datamodels.open(file_path) as model:
        assert model.meta.telescope == 'XOMAN'
    asdf.get_config().validate_on_read = True


# Testing all reference file schemas
def test_reference_file_model_base(tmp_path):
    # Set temporary asdf file
    file_path = tmp_path / "test.asdf"

    # Get all reference file classes
    tags = [t for t in DATAMODEL_EXTENSIONS[0].tags if "/reference_files/" in t._tag_uri]
    schemas = []
    for tag in tags:
        schema = asdf.schema.load_schema(tag._schema_uri)
        # Check that schema references common reference schema
        allofs = schema['properties']['meta']['allOf']
        found_common = False
        for item in allofs:
            if item == EXPECTED_COMMON_REFERENCE:
                found_common = True
        if not found_common:
            raise ValueError("Reference schema does not include ref_common")

def test_opening_flat_ref(tmp_path):
    # First make test reference file
    file_path = tmp_path / 'testflat.asdf'
    util.mk_flat(file_path)
    flat = datamodels.open(file_path)
    assert flat.meta.instrument.optical_element == 'F158'
    assert isinstance(flat, datamodels.FlatRefModel)


# FlatModel tests
def test_flat_model(tmp_path):
    # Set temporary asdf file
    file_path = tmp_path / "test.asdf"

    meta = {}
    util.add_ref_common(meta)
    meta['reftype'] = "FLAT"
    flatref = converters.FlatRef()
    flatref['meta'] = meta
    flatref.meta.instrument['optical_element'] = 'F062'
    shape = (4096, 4096)
    flatref['data'] = np.zeros(shape, dtype=np.float32)
    flatref['dq'] = np.zeros(shape, dtype=np.uint32)
    flatref['err'] = np.zeros(shape, dtype=np.float32)

    # Testing flat file asdf file
    with asdf.AsdfFile(meta) as af:
        af.tree = {'roman': flatref}
        af.write_to(file_path)

        # Test that asdf file opens properly
        with datamodels.open(file_path) as model:
            with pytest.warns(None):
                model.validate()

            # Confirm that asdf file is opened as flat file model
            assert isinstance(model, datamodels.FlatRefModel)

## not sure what the following is supposed to ensure PG

# def test_meta_date_management(tmp_path):
#     model = datamodels.RomanDataModel({
#         "meta": {
#             "date": Time("2000-01-01T00:00:00.000"),
#             "instrument": {"name": "WFI", "detector": "WFI01", "optical_element": "F062"},
#             "telescope": "ROMAN",
#         }
#     })
#     assert model.meta.date == Time("2000-01-01T00:00:00.000")
#     model.save(str(tmp_path/"test.asdf"))
#     assert abs((Time.now() - model.meta.date).value) < 1.0

#     model = datamodels.RomanDataModel()
#     assert abs((Time.now() - model.meta.date).value) < 1.0

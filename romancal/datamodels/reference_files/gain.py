from .referencefile import ReferenceFileModel

class GainModel(ReferenceFileModel):
    """
    A data model for 2D gain.

    Parameters
    __________
    data : numpy float32 array
         The gain
    """
    schema_url = "http://stsci.edu/schemas/roman_datamodel/reference_files/gain.schema"

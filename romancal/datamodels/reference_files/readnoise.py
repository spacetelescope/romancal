from .referencefile import ReferenceFileModel


class ReadNoiseModel(ReferenceFileModel):
    """
    A data model for 2D read noise.

    Parameters
    __________
    data : numpy float32 array
         The read noise
    """
    schema_url = "http://stsci.edu/schemas/roman_datamodel/reference_files/readnoise.schema"

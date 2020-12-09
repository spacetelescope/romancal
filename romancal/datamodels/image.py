from .core import RomanDataModel

class ImageModel(RomanDataModel):
    """
    A data model for 2D images.

    Parameters
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array

    zeroframe : numpy float32 array
         Zeroframe array

    var_poisson : numpy float32 array
         variance due to poisson noise

    var_rnoise : numpy float32 array
         variance due to read noise

    area : numpy float32 array
         Pixel area map array
    """
    schema_url = "http://stsci.edu/schemas/roman_datamodel/image_files/image.schema"

    def __init__(self, init=None, **kwargs):
        super().__init__(init=init, **kwargs)


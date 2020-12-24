from ..core import RomanDataModel

class ImageModel(RomanDataModel):
    """
    A data model for 2D images.

    """
    schema_url = "http://stsci.edu/schemas/roman_datamodel/image_files/image.schema"

    def __init__(self, init=None, **kwargs):
        super().__init__(init=init, **kwargs)

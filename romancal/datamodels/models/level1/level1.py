from ..core import RomanDataModel

class Level1FileModel(RomanDataModel):
    """
    Data model for level 1 data files.
    """
    schema_url = "http://stsci.edu/schemas/roman_datamodel/level1.schema"

    def __init__(self, init=None, **kwargs):
        super().__init__(init=init, **kwargs)


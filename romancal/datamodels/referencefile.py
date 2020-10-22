from stdatamodels import DataModel

class ReferenceFileModel(DataModel):
    """
    Data model for reference files.
    """
    schema_url = "http://stsci.edu/schemas/roman_datamodel/reference_files/referencefile.schema"

    def __init__(self, init=None, **kwargs):
        super().__init__(init=init, **kwargs)

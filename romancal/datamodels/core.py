from stdatamodels import DataModel


class RomanDataModel(DataModel):
    """
    Base class for Roman data models.  Also serves as a generic
    model class for data that does not identify its model type.
    """
    schema_url = "http://stsci.edu/schemas/roman_datamodel/core.schema"


from .referencefile import ReferenceFileModel
from ..dynamicdq import dynamic_mask

class DarkModel(ReferenceFileModel):
    """
    A data model for dark reference files.

    Parameters
    __________
    data : numpy float32 array
         Dark current array

    dq : numpy uint16 array
         2-D data quality array for all planes

    err : numpy float32 array
         Error array
    """
    schema_url = "http://stsci.edu/schemas/roman_datamodel/reference_files/dark.schema"

    def __init__(self, init=None, **kwargs):
        super().__init__(init=init, **kwargs)

        # Saved for DQ code implementation
        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err

from .referencefile import ReferenceFileModel
from ..dynamicdq import dynamic_mask

class FlatModel(ReferenceFileModel):
    """
    A data model for 2D flat-field images.

    Parameters
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint16 array
         Data quality array

    err : numpy float32 array
         Error array
    """
    schema_url = "http://stsci.edu/schemas/roman_datamodel/reference_files/flat.schema"

    def __init__(self, init=None, **kwargs):
        super().__init__(init=init, **kwargs)

        # Saved for DQ code implementation
        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.err = self.err

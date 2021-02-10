from .core import RomanDataModel

__all__ = ["RampModel"]


class RampModel(RomanDataModel):
    """
    A data model for 3D ramps.

    Parameters
    __________
    data : numpy float32 array
         The science data

    pixeldq : numpy uint32 array
         2-D data quality array for all planes

    groupdq : numpy uint8 array
         3-D data quality array for each plane

    err : numpy float32 array
         Error array

    zeroframe : numpy float32 array
         Zeroframe array

    group : numpy table
         group parameters table

    """
    schema_url = "http://stsci.edu/schemas/roman_datamodel/ramp.schema"

    def __init__(self, init=None, **kwargs):
        super(RampModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.pixeldq = self.pixeldq
        self.groupdq = self.groupdq
        self.err = self.err

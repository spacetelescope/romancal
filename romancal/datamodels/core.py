from stdatamodels import DataModel
from astropy.time import Time


class RomanDataModel(DataModel):
    """
    Base class for Roman data models.  Also serves as a generic
    model class for data that does not identify its model type.
    """
    schema_url = "http://stsci.edu/schemas/roman_datamodel/core.schema"

    # CRDS observatory code for Roman.  Used by stpipe when selecting
    # references.
    crds_observatory = "roman"

    def get_crds_parameters(self):
        """
        Get parameters used by CRDS to select references for this model.

        Returns
        -------
        dict
        """
        return {
            key: val for key, val in self.to_flat_dict(include_arrays=False).items()
            if isinstance(val, (str, int, float, complex, bool))
        }

    def on_init(self, init):
        """
        Hook invoked by the base class before returning a newly
        created model instance.
        """
        super().on_init(init)

        if self.meta.date is None:
            self.meta.date = Time(Time.now(), format="isot")

    def on_save(self, init):
        """
        Hook invoked by the base class before writing a model
        to a file.
        """
        super().on_save(init)

        self.meta.date = Time(Time.now(), format="isot")

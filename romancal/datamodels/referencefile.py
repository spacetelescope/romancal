from stdatamodels import DataModel
# Saved for validation
# from .validate import ValidationWarning


class RomanReferenceFileModel(DataModel):
    """
    Data model for reference files.
    """
    schema_url = "http://stsci.edu/schemas/roman_datamodel/reference_files/referencefile.schema"

    def __init__(self, init=None, **kwargs):
        super(RomanReferenceFileModel, self).__init__(init=init, **kwargs)
        # self._no_asdf_extension = True
        self.meta.telescope = "Roman"

    # Saved for validation
    # def validate(self):
    #     """
    #     Convenience function to be run when files are created.
    #     Checks that required reference file keywords are set.
    #     """
    #     to_fix = []
    #     to_check = ['telescope','description', 'reftype', 'author', 'pedigree', 'useafter']
    #     for field in to_check:
    #         if getattr(self.meta, field) is None:
    #             to_fix.append(field)
    #     if self.meta.instrument.name is None:
    #         to_fix.append('instrument.name')
    #     if self.meta.telescope != 'JWST':
    #         to_fix.append('telescope')
    #     if to_fix:
    #         Saved for validation
    #         self.print_err(f'Model.meta is missing values for {to_fix}')
    #     super().validate()

    # Saved for validation
    # def print_err(self, message):
    #     if self._strict_validation:
    #         raise ValueError(message)
    #     else:
    #         warnings.warn(message, ValidationWarning)


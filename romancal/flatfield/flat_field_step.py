from ..stpipe import Step
from .. import datamodels
from . import flat_field


__all__ = ["FlatFieldStep"]


class FlatFieldStep(Step):
    """Flat-field a science image using a flatfield reference image.
    """

    spec = """
        save_interpolated_flat = boolean(default=False) # DG - for Roman ?
    """
    reference_file_types = ["flat"]

    # Define suffix for optional saved output of the interpolated flat for NRS
    # DG - for Roman ?
    flat_suffix = 'interpolatedflat'

    def process(self, input):

        input_model = datamodels.open(input)
        exposure_type = input_model.meta.exposure.type.upper()

        self.log.debug("Input is {} of exposure type {}".format(
            input_model.__class__.__name__, exposure_type))

        # Get reference file paths
        reference_file_names = {}
        reftype = "flat"
        reffile = self.get_reference_file(input_model, "flat")
        reference_file_names[reftype] = reffile if reffile != 'N/A' else None

        # Define mapping between reftype and datamodel type
        model_type = dict(
            flat=datamodels.FlatModel,
            )

        # Open the relevant reference files as datamodels
        reference_file_models = {}

        if reffile is not None:
            reference_file_models[reftype] = model_type[reftype](reffile)
            self.log.debug('Using %s ref file: %s', reftype.upper(), reffile)
        else:
            reference_file_models[reftype] = None

        # Do the flat-field correction
        output_model, interpolated_flats = flat_field.do_correction(
            input_model,
            **reference_file_models,
            )

        # Close the input and reference files
        input_model.close()
        try:
            for model in reference_file_models.values():
                model.close()
        except AttributeError:
            pass

        if self.save_interpolated_flat and interpolated_flats is not None:
            self.log.info("Writing interpolated flat field.")
            self.save_model(interpolated_flats, suffix=self.flat_suffix)
            interpolated_flats.close()

        return output_model

    def skip_step(self, input_model):
        """Set the calibration switch to SKIPPED.

        This method makes a copy of input_model, sets the calibration
        switch for the flat_field step to SKIPPED in the copy, closes
        input_model, and returns the copy.
        """

        result = input_model.copy()
        result.meta.cal_step.flat_field = "SKIPPED"
        input_model.close()
        return result

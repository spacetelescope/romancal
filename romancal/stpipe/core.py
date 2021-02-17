from stpipe import Step, Pipeline

from .. import datamodels


class RomanStep(Step):
    """
    Base class for Roman calibration pipeline steps.
    """
    spec = """
    output_ext =  string(default='.asdf')    # Default type of output
    """

    @classmethod
    def _datamodels_open(cls, init, **kwargs):
        """
        Provide access to this package's datamodels.open function
        so that the stpipe infrastructure knows how to instantiate
        models.
        """
        return datamodels.open(init, **kwargs)

    def finalize_result(self, model, reference_files_used):
        """
        Hook that allows the Step to set metadata on the output model
        before save.

        Parameters
        ----------
        model : stdatamodels.DataModel
            Output model.

        reference_files_used : list of tuple(str, str)
            List of reference files used.  The first element of each tuple
            is the reftype code, the second element is the filename.
        """
        #  JWST uses this to add the cal code and CRDS software versions.
        pass

    def record_step_status(self, model, step_name, success=True):
        """
        Record step completion status in the model's metadata.

        Parameters
        ----------
        model : stdatamodels.DataModel
            Output model.
        step_name : str
            Calibration step name.
        success : bool
            If True, then the step was run successfully.
        """
        # JWST sets model.meta.cal_step.<step name> here.  Roman
        # may do the same, depending on how the metadata format
        # turns out.  Seems like we might be able to combine this
        # with finalize_result somehow.
        pass

    def remove_suffix(self, name):
        """
        Remove any Roman step-specific suffix from the given filename.

        Parameters
        ----------
        name : str
            Filename.

        Returns
        -------
        str
            Filename with step suffix removed.
        """
        # JWST maintains a list of relevant suffixes that is monitored
        # by tests to be up-to-date.  Roman will likely need to do
        # something similar.
        return name, "_"


# RomanPipeline needs to inherit from Pipeline, but also
# be a subclass of RomanStep so that it will pass checks
# when constructing a pipeline using RomanStep class methods.
class RomanPipeline(Pipeline, RomanStep):
    pass

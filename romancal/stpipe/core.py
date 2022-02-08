"""
Roman Calibration Pipeline base class
"""
import logging
import time

from stpipe import Step, Pipeline

import roman_datamodels as rdm
from roman_datamodels.datamodels import ImageModel
from ..lib.suffix import remove_suffix


_LOG_FORMATTER = logging.Formatter(
    "%(asctime)s.%(msecs)03dZ :: %(name)s :: %(levelname)s :: %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S"
)
_LOG_FORMATTER.converter = time.gmtime


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
        return rdm.open(init, **kwargs)

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
        if isinstance(model, ImageModel):
            for log_record in self.log_records:
                model.cal_logs.append(_LOG_FORMATTER.format(log_record))

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
        return remove_suffix(name)


# RomanPipeline needs to inherit from Pipeline, but also
# be a subclass of RomanStep so that it will pass checks
# when constructing a pipeline using RomanStep class methods.
class RomanPipeline(Pipeline, RomanStep):
    pass

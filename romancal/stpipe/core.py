"""
Roman Calibration Pipeline base class
"""

import importlib.metadata
import logging
import time
from pathlib import Path

import roman_datamodels as rdm
from roman_datamodels.datamodels import ImageModel, MosaicModel
from stpipe import Pipeline, Step, crds_client

from romancal.datamodels.library import ModelLibrary

from ..lib.suffix import remove_suffix

_LOG_FORMATTER = logging.Formatter(
    "%(asctime)s.%(msecs)03dZ :: %(name)s :: %(levelname)s :: %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
_LOG_FORMATTER.converter = time.gmtime


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class RomanStep(Step):
    """
    Base class for Roman calibration pipeline steps.
    """

    spec = """
    output_ext =  string(default='.asdf')    # Default type of output
    """

    _log_records_formatter = _LOG_FORMATTER

    @classmethod
    def _datamodels_open(cls, init, **kwargs):
        """
        Provide access to this package's datamodels.open function
        so that the stpipe infrastructure knows how to instantiate
        models and containers.
        """
        if isinstance(init, str):
            init = Path(init)
        if isinstance(init, Path):
            ext = init.suffix.lower()
            if ext == ".asdf":
                return rdm.open(init, **kwargs)
            if ext in (".json", ".yaml"):
                return ModelLibrary(init, **kwargs)
        if isinstance(init, rdm.DataModel):
            return rdm.open(init, **kwargs)
        if isinstance(init, ModelLibrary):
            return ModelLibrary(init)
        raise TypeError(f"Invalid input: {init}")

    @classmethod
    def _get_crds_parameters(cls, dataset):
        crds_parameters, crds_observatory = super()._get_crds_parameters(dataset)
        if "roman.meta.instrument.detector" not in crds_parameters:
            crds_parameters["roman.meta.instrument.detector"] = "WFI02"
        if "roman.meta.exposure.start_time" not in crds_parameters:
            crds_parameters["roman.meta.exposure.start_time"] = crds_parameters[
                "roman.meta.coadd_info.time_first"
            ]
        return crds_parameters, crds_observatory

    def finalize_result(self, model, reference_files_used):
        """
        Hook that allows the Step to set metadata on the output model
        before save.

        Parameters
        ----------
        model : roman_datamodels.datamodels.DataModel
            Output model.

        reference_files_used : list of tuple(str, str)
            List of reference files used.  The first element of each tuple
            is the reftype code, the second element is the filename.
        """

        model.meta.calibration_software_version = importlib.metadata.version("romancal")

        if isinstance(model, ImageModel | MosaicModel):
            # convert to model.cal_logs type to avoid validation errors
            model.meta.cal_logs = type(model.meta.cal_logs)(self.log_records)

        if len(reference_files_used) > 0:
            if not hasattr(model.meta, "ref_file"):
                log.error(
                    f"Model[{model}] is missing meta.ref_file. {reference_files_used} will not be recorded"
                )
            else:
                for ref_name, ref_file in reference_files_used:
                    if hasattr(model.meta.ref_file, ref_name):
                        setattr(model.meta.ref_file, ref_name, ref_file)
                        # getattr(model.meta.ref_file, ref_name).name = ref_file
                model.meta.ref_file.crds.version = crds_client.get_svn_version()
                model.meta.ref_file.crds.context = crds_client.get_context_used(
                    model.crds_observatory
                )

                # this will only run if 'parent' is none, which happens when an individual
                # step is being run or if self is a RomanPipeline and not a RomanStep.
                if self.parent is None:
                    log.info(
                        f"Results used CRDS context: {model.meta.ref_file.crds.context}"
                    )

    def record_step_status(self, model, step_name, success=True):
        """
        Record step completion status in the model's metadata.

        Parameters
        ----------
        model : roman_datamodels.datamodels.DataModel
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

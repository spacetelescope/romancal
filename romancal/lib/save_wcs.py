"""Step-level utility to create/save WfiWcsModel"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from roman_datamodels.datamodels import WfiWcsModel

from romancal.datamodels import ModelLibrary

if TYPE_CHECKING:
    from romancal.stpipe import RomanStep

# Define logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def save_wfiwcs(step: RomanStep, result: ModelLibrary, force: bool = False):
    """Create and save the WfiWcs products

    Parameters
    ----------
    result : ModelLibrary
         The final L2 models

    force : boolean
        Regardless of whether ``save_results`` is `False`
        and no ``output_file`` is specified, try saving.
    """
    if isinstance(result, ModelLibrary):
        list(
            result.map_function(
                lambda model, index: save_wfiwcs(step, model, force), modify=False
            )
        )
    else:
        # this is a datamodel
        try:
            wfiwcs = WfiWcsModel.from_model_with_wcs(result)
        # TODO make this a more specific exception
        except ValueError:
            log.info(
                f"No WCS information for model {result}. Now `_wcs` product will be created."
            )
            return
        step.finalize_result(wfiwcs, [])
        step.save_model(wfiwcs, suffix="wcs", force=force)

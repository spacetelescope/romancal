"""Step-level utility to create/save WfiWcsModel"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from roman_datamodels.datamodels import WfiWcsModel

if TYPE_CHECKING:
    from romancal.datamodels import ModelLibrary
    from romancal.stpipe import RomanStep

# Define logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def save_wfiwcs(step: RomanStep, lib: ModelLibrary, force: bool = False):
    """Create and save the WfiWcs products

    Parameters
    ----------
    lib : ModelLibrary
         The final L2 models

    force : boolean
        Regardless of whether ``save_results`` is `False`
        and no ``output_file`` is specified, try saving.
    """
    log.info("Writing the WCS files...")
    with lib:
        for model in lib:
            try:
                wfiwcs = WfiWcsModel.from_model_with_wcs(model)
            except ValueError:
                log.info(
                    f"No WCS information for model {model}. Now `_wcs` product will be created."
                )
                lib.shelve(model)
                continue
            step.finalize_result(wfiwcs, [])
            step.save_model(wfiwcs, suffix="wcs", force=force)
            lib.shelve(model)

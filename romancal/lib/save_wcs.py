"""Step-level utility to create/save WfiWcsModel"""

import logging

from roman_datamodels.datamodels import WfiWcsModel

# Define logging
log = logging.getLogger()
log.setLevel(logging.DEBUG)


def save_wfiwcs(step, lib):
    """Create and save the WfiWcs products

    Parameters
    ----------
    lib : ModelLibrary
         The final L2 models
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
            step.save_model(wfiwcs, suffix="wcs")
            lib.shelve(model)

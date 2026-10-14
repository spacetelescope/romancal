"""
Roman step for sky matching.
"""

from __future__ import annotations

import logging
from copy import deepcopy
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from astropy.nddata.bitmask import bitfield_to_boolean_mask, interpret_bit_flags
from roman_datamodels.dqflags import pixel
from stcal.skymatch import SkyImage, SkyStats, skymatch

from romancal.datamodels.fileio import open_dataset
from romancal.skymatch.offset_io import (
    load_offset_table,
    lookup_offset,
    resolve_offsets_path_from_asn,
)
from romancal.stpipe import RomanStep

if TYPE_CHECKING:
    from typing import ClassVar

__all__ = ["SkyMatchStep"]

log = logging.getLogger(__name__)


class SkyMatchStep(RomanStep):
    """
    SkyMatchStep: Subtraction or equalization of sky background in science images.
    """

    class_alias = "skymatch"

    spec = """
        # General sky matching parameters:
        skymethod = option('local', 'global', 'match', 'global+match', default='match') # sky computation method
        match_down = boolean(default=True) # adjust sky to lowest measured value?
        subtract = boolean(default=False) # subtract computed sky from image data?

        # User / ASN offset table (parquet). When set, stcal skymatch is skipped.
        offset_file = string(default='')  # parquet offsets table; overrides ASN field when set

        # Image's bounding polygon parameters:
        stepsize = integer(default=None) # Max vertex separation

        # Sky statistics parameters:
        skystat = option('median', 'midpt', 'mean', 'mode', default='mode') # sky statistics
        dqbits = string(default='~DO_NOT_USE') # "good" DQ bits
        lower = float(default=None) # Lower limit of "good" pixel values
        upper = float(default=None) # Upper limit of "good" pixel values
        nclip = integer(min=0, default=5) # number of sky clipping iterations
        lsigma = float(min=0.0, default=4.0) # Lower clipping limit, in sigma
        usigma = float(min=0.0, default=4.0) # Upper clipping limit, in sigma
        binwidth = float(min=0.0, default=0.1) # Bin width for 'mode' and 'midpt' `skystat`, in sigma
    """

    reference_file_types: ClassVar = []

    def process(self, dataset):
        library = open_dataset(
            dataset, update_version=self.update_version, as_library=True
        )

        offsets_path = self._resolve_offset_file(library)
        if offsets_path is not None:
            return self._process_from_offset_file(library, offsets_path)

        self._dqbits = interpret_bit_flags(self.dqbits, flag_name_map=pixel)

        # set sky statistics:
        self._skystat = SkyStats(
            skystat=self.skystat,
            lower=self.lower,
            upper=self.upper,
            nclip=self.nclip,
            lsig=self.lsigma,
            usig=self.usigma,
            binwidth=self.binwidth,
        )

        # create a list of "Sky" Images and/or Groups:
        sky_images = []
        with library:
            for index, model in enumerate(library):
                sky_images.append(self._imodel2skyim(model, index))
                library.shelve(model, index, modify=False)

            # match/compute sky values:
            skymatch(
                sky_images,
                skymethod=self.skymethod,
                match_down=self.match_down,
                subtract=self.subtract,
            )

            # set sky background value in each image's meta:
            for sky_image in sky_images:
                index = sky_image.meta["index"]
                image = library.borrow(index)
                self._set_sky_background(
                    image,
                    level=sky_image.sky if sky_image.sky is not None else 0.0,
                    method=str(self.skymethod),
                    complete=bool(sky_image.is_sky_valid),
                    data=sky_image.image if self.subtract else None,
                )
                library.shelve(image, index)

        return library

    def _resolve_offset_file(self, library) -> Path | None:
        """Return parquet offsets path from step param or ASN, if any."""
        if self.offset_file:
            path = Path(str(self.offset_file)).expanduser()
            if not path.is_file():
                raise FileNotFoundError(f"offset_file not found: {path}")
            return path.resolve()

        asn = getattr(library, "asn", None)
        if not asn:
            return None
        # AbstractModelLibrary stores the ASN directory for member path joins.
        asn_dir = getattr(library, "_asn_dir", None)
        asn_path = Path(asn_dir) / "asn.json" if asn_dir else None
        path = resolve_offsets_path_from_asn(asn, asn_path=asn_path)
        if path is None:
            return None
        if not path.is_file():
            raise FileNotFoundError(f"ASN skymatch_offsets_file not found: {path}")
        return path

    def _process_from_offset_file(self, library, offsets_path: Path):
        """Apply pipeline levels from a parquet offsets table."""
        log.info("Using skymatch offsets from %s", offsets_path)
        table = load_offset_table(offsets_path)
        method = "file"

        with library:
            for index, model in enumerate(library):
                filename = getattr(model.meta, "filename", None) or f"model_{index}"
                try:
                    level = lookup_offset(filename, table)
                except KeyError as err:
                    raise ValueError(
                        f"No offset entry for library member {filename!r} "
                        f"in {offsets_path}"
                    ) from err

                self._set_sky_background(
                    model,
                    level=float(level),
                    method=method,
                    complete=True,
                    data=None,
                )
                if self.subtract:
                    model.data = np.asarray(model.data) - float(level)
                library.shelve(model, index)

        return library

    def _imodel2skyim(self, image_model, index):
        background_meta = image_model.meta.get("background", {})

        if self._dqbits is None:
            dqmask = np.isfinite(image_model.data)
        else:
            dqmask = bitfield_to_boolean_mask(
                image_model.dq, self._dqbits, good_mask_value=True, dtype="bool"
            ) & np.isfinite(image_model.data)

        # see if 'skymatch' was previously run and raise an exception
        # if 'subtract' mode has changed compared to the previous pass:
        if background_meta.get("subtracted", False):
            if "level" not in background_meta:
                raise ValueError(
                    "Background level was subtracted but the "
                    "'level' property is undefined (None)."
                )
            if self.subtract:
                raise ValueError(
                    "'subtract' step's specification is "
                    "inconsistent with background info already "
                    f"present in image '{image_model.meta.filename:s}' meta."
                )

        wcs = deepcopy(image_model.meta.wcs)

        sky_im = SkyImage(
            image=image_model.data,
            wcs_fwd=wcs.forward_transform,
            wcs_inv=wcs.backward_transform,
            mask=dqmask,
            sky_id=image_model.meta.filename,
            skystat=self._skystat,
            stepsize=self.stepsize,
            meta={"index": index},
        )

        if self.subtract and background_meta.get("level") is not None:
            sky_im.sky = background_meta["level"]

        return sky_im

    def _set_sky_background(self, image, *, level, method, complete, data=None):
        if "background" not in image.meta:
            image.meta["background"] = {}

        image.meta.background.method = str(method)
        image.meta.background.subtracted = self.subtract
        image.meta.background.level = float(level)

        if complete and self.subtract and data is not None:
            image.data[...] = data[...]

        image.meta.cal_step.skymatch = "COMPLETE" if complete else "SKIPPED"

"""
Roman step for sky matching.
"""

from __future__ import annotations

import logging
from copy import deepcopy
from typing import TYPE_CHECKING
from urllib.parse import urlparse

import numpy as np
from astropy.nddata.bitmask import bitfield_to_boolean_mask, interpret_bit_flags
from roman_datamodels.dqflags import pixel
from stcal.skymatch import SkyImage, SkyStats, skymatch

from romancal.datamodels.fileio import open_dataset
from romancal.skymatch.offset_io import (
    load_offset_table,
    lookup_offset,
    offsets_path_from_asn,
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
        skymethod = option('local', 'global', 'match', 'global+match', 'offset_file', default='match') # sky computation method
        match_down = boolean(default=True) # adjust sky to lowest measured value?
        subtract = boolean(default=False) # subtract computed sky from image data?

        # Precomputed levels parquet (used only when skymethod='offset_file'):
        offset_file = string(default='')  # parquet levels table; overrides ASN field when set

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

        if self.skymethod == "offset_file":
            offsets_path = self._resolve_offset_file(library)
            if not offsets_path:
                raise ValueError(
                    "skymethod='offset_file' requires offset_file or ASN "
                    "skymatch_offsets_file"
                )
            return self._process_from_offset_file(library, offsets_path)

        if self.offset_file:
            log.info(
                "Ignoring offset_file=%r because skymethod=%r (use "
                "skymethod='offset_file' to apply precomputed levels)",
                self.offset_file,
                self.skymethod,
            )

        self._dqbits = interpret_bit_flags(self.dqbits, flag_name_map=pixel)

        self._skystat = SkyStats(
            skystat=self.skystat,
            lower=self.lower,
            upper=self.upper,
            nclip=self.nclip,
            lsig=self.lsigma,
            usig=self.usigma,
            binwidth=self.binwidth,
        )

        sky_images = []
        with library:
            for index, model in enumerate(library):
                sky_images.append(self._imodel2skyim(model, index))
                library.shelve(model, index, modify=False)

            skymatch(
                sky_images,
                skymethod=self.skymethod,
                match_down=self.match_down,
                subtract=self.subtract,
            )

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

    def _resolve_offset_file(self, library) -> str | None:
        """Return offsets path from step param or top-level ASN field.

        Uses ``ModelLibrary._member_to_filename`` so relative paths and URLs
        follow the same rules as association science members.
        """
        if self.offset_file:
            raw = str(self.offset_file).strip()
        else:
            asn = getattr(library, "asn", None)
            if not asn:
                return None
            raw = offsets_path_from_asn(asn)
            if not raw:
                return None

        return library._member_to_filename({"expname": raw})

    def _process_from_offset_file(self, library, offsets_path: str):
        """Apply pipeline levels from a parquet offsets table."""
        # Path is already library-resolved in _resolve_offset_file.
        path = str(offsets_path)
        log.info("Using skymatch offsets from %s", path)

        try:
            table = load_offset_table(path)
        except FileNotFoundError as err:
            raise FileNotFoundError(f"offset_file not found: {path}") from err
        except OSError as err:
            if not urlparse(str(path)).scheme:
                raise FileNotFoundError(f"offset_file not found: {path}") from err
            raise

        with library:
            for index, model in enumerate(library):
                background_meta = model.meta.get("background", {}) or {}
                if background_meta.get("subtracted", False):
                    raise ValueError(
                        "skymethod='offset_file' cannot run when background was "
                        f"already subtracted on {model.meta.filename!r}"
                    )

                filename = getattr(model.meta, "filename", None) or f"model_{index}"
                try:
                    level = lookup_offset(filename, table)
                except KeyError as err:
                    raise ValueError(
                        f"No offset entry for library member {filename!r} in {path}"
                    ) from err

                self._set_sky_background(
                    model,
                    level=float(level),
                    method="offset_file",
                    complete=True,
                    data=None,
                    offset_file=str(path),
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

    def _set_sky_background(
        self, image, *, level, method, complete, data=None, offset_file=None
    ):
        if "background" not in image.meta:
            image.meta["background"] = {}

        image.meta.background.method = str(method)
        image.meta.background.subtracted = self.subtract
        image.meta.background.level = float(level)
        if offset_file is not None:
            image.meta.background["offset_file"] = str(offset_file)

        if complete and self.subtract and data is not None:
            image.data[...] = data[...]

        image.meta.cal_step.skymatch = "COMPLETE" if complete else "SKIPPED"

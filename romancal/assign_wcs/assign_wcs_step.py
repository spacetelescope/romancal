"""
Assign a gWCS object to a science image.

"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import gwcs.coordinate_frames as cf
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import bind_bounding_box
from gwcs.wcs import WCS, Step
from roman_datamodels import datamodels as rdm

from ..stpipe import RomanStep
from . import pointing
from .utils import add_s_region, wcs_bbox_from_shape

if TYPE_CHECKING:
    from typing import ClassVar

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["AssignWcsStep", "load_wcs"]


class AssignWcsStep(RomanStep):
    """Assign a gWCS object to a science image."""

    class_alias = "assign_wcs"

    reference_file_types: ClassVar = ["distortion"]

    def process(self, input):
        reference_file_names = {}
        if isinstance(input, rdm.DataModel):
            input_model = input
        else:
            input_model = rdm.open(input)

        for reftype in self.reference_file_types:
            log.info(f"reftype, {reftype}")
            reffile = self.get_reference_file(input_model, reftype)
            # Check for a valid reference file
            if reffile == "N/A":
                self.log.warning("No DISTORTION reference file found")
                self.log.warning("Assign WCS step will be skipped")
                result = input_model.copy()
                result.meta.cal_step.assign_wcs = "SKIPPED"
                return result

            reference_file_names[reftype] = reffile if reffile else ""
        log.info("Using reference files: %s for assign_wcs", reference_file_names)
        result = load_wcs(input_model, reference_file_names)

        if self.save_results:
            try:
                self.suffix = "assignwcs"
            except AttributeError:
                self["suffix"] = "assignwcs"

        return result


def load_wcs(input_model, reference_files=None):
    """Create a gWCS object and store it in ``Model.meta``.

    Parameters
    ----------
    input_model : `~roman_datamodels.datamodels.WfiImage`
        The exposure.
    reference_files : dict
        A dict {reftype: reference_file_name} containing all
        reference files that apply to this exposure.

    Returns
    -------
    output_model : `~roman_datamodels.ImageModel`
        The input image file with attached gWCS object.
        The input_model is modified in place.
    """
    output_model = input_model

    if reference_files is not None:
        for ref_type, ref_file in reference_files.items():
            reference_files[ref_type] = (
                ref_file if ref_file not in ["N/A", ""] else None
            )
    else:
        reference_files = {}

    # Frames
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(
        name="v2v3",
        axes_order=(0, 1),
        axes_names=("v2", "v3"),
        unit=(u.arcsec, u.arcsec),
    )
    v2v3vacorr = cf.Frame2D(
        name="v2v3vacorr",
        axes_order=(0, 1),
        axes_names=("v2", "v3"),
        unit=(u.arcsec, u.arcsec),
    )
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name="world")

    # Transforms between frames
    distortion = wfi_distortion(output_model, reference_files)
    tel2sky = pointing.v23tosky(output_model)

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref,
    )

    pipeline = [
        Step(detector, distortion),
        Step(v2v3, va_corr),
        Step(v2v3vacorr, tel2sky),
        Step(world, None),
    ]
    wcs = WCS(pipeline)
    if wcs.bounding_box is None:
        wcs.bounding_box = wcs_bbox_from_shape(output_model.data.shape)

    output_model.meta["wcs"] = wcs

    # update S_REGION
    add_s_region(output_model)

    output_model.meta.cal_step["assign_wcs"] = "COMPLETE"

    return output_model


def wfi_distortion(model, reference_files):
    """
    Create the "detector" to "v2v3" transform for WFI

    Parameters
    ----------
    model : `~roman_datamodels.datamodels.WfiImage`
        The data model for processing
    reference_files : dict
        A dict {reftype: reference_file_name} containing all
        reference files that apply to this exposure.

    Returns
    -------
    The transform model
    """

    dist = rdm.DistortionRefModel(reference_files["distortion"])
    transform = dist.coordinate_distortion_transform

    try:
        bbox = transform.bounding_box.bounding_box(order="F")
    except NotImplementedError:
        # Check if the transform in the reference file has a ``bounding_box``.
        # If not set a ``bounding_box`` equal to the size of the image after
        # assembling all distortion corrections.
        bbox = None
    dist.close()

    bind_bounding_box(
        transform,
        wcs_bbox_from_shape(model.data.shape) if bbox is None else bbox,
        order="F",
    )

    return transform

"""
Assign a gWCS object to a science image.

"""
import logging

from astropy import coordinates as coord
from astropy import units as u
import gwcs.coordinate_frames as cf

from gwcs.wcs import WCS, Step

from roman_datamodels import datamodels as rdm
from ..stpipe import RomanStep
from . import pointing
from .utils import wcs_bbox_from_shape


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["AssignWcsStep", "load_wcs"]


class AssignWcsStep(RomanStep):
    """ Assign a gWCS object to a science image.
    """

    reference_file_types = ['distortion']

    def process(self, input):
        reference_file_names = {}
        with rdm.open(input) as input_model:
            for reftype in self.reference_file_types:
                log.info(f'reftype, {reftype}')
                reffile = self.get_reference_file(input_model, reftype)
                reference_file_names[reftype] = reffile if reffile else ""
            log.debug(f'reference files used in assign_wcs: {reference_file_names}')
            result = load_wcs(input_model, reference_file_names)
        return result


def load_wcs(dmodel, reference_files):
    """ Create a gWCS object and store it in ``Model.meta``.

    Parameters
    ----------
    dmodel : `~roman_datamodels.datamodels.WfiImage`
        The exposure.
    reference_files : dict
        A dict {reftype: reference_file_name} containing all
        reference files that apply to this exposure.

    Returns
    -------
    dmodel : `~roman_datamodels.ImageModel`
        The input image file with attached gWCS object.
        The data is not modified.
    """
    if reference_files:
        for ref_type, ref_file in reference_files.items():
            if ref_file not in ["N/A", ""]:
                reference_files[ref_type] = ref_file
            else:
                reference_files[ref_type] = None

    # Frames
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), axes_names=('v2', 'v3'), unit=(u.arcsec, u.arcsec))
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name='world')

    # Transforms between frames
    distortion = wfi_distortion(dmodel, reference_files)
    tel2sky = pointing.v23tosky(dmodel)

    pipeline = [Step(detector, distortion),
                Step(v2v3, tel2sky),
                Step(world, None)]
    wcs = WCS(pipeline)
    if wcs.bounding_box is None:
        wcs.bounding_box = wcs_bbox_from_shape(dmodel.data.shape)

    dmodel.meta['wcs'] = wcs
    dmodel.meta.cal_step['assign_wcs'] = 'COMPLETE'

    return dmodel

def wfi_distortion(input_model, reference_files):
    """
    Create the "detector" to "v2v3" transform for WFI

    Parameters
    ----------
    input_model : `~roman_datamodels.datamodels.WfiImage`
        The data model for processing
    reference_files : dict
        A dict {reftype: reference_file_name} containing all
        reference files that apply to this exposure.

    Returns
    -------
    The transform model
    """

    dist = rdm.DistortionRefModel(reference_files['distortion'])
    transform = dist.model

    try:
        bbox = transform.bounding_box
    except NotImplementedError:
        # Check if the transform in the reference file has a ``bounding_box``.
        # If not set a ``bounding_box`` equal to the size of the image after
        # assembling all distortion corrections.
        bbox = None
    dist.close()

    if bbox is None:
        transform.bounding_box = wcs_bbox_from_shape(input_model.data.shape)
    else:
        transform.bounding_box = bbox

    return transform
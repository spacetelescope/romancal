from astropy.modeling.models import Mapping

from .assign_wcs_step import AssignWcsStep

__all__ = ["AssignWcsStep", "_distorted_to_undistorted"]


def _distorted_to_undistorted(wcsobj):
    """
    This is a convenience function for Build 17 which returns the transform
    from "detector" to "undistorted pixels".
    TODO: Add an 'undistorted' coordinate frame. The transform from detector to
    'undistorted' is currently in the SIAF from 'detector' to 'Ideal' frame.
    The split is done in the code now. In the future this should be done in the
    distortion reference file and become part of the WCS.
    """
    transform = wcsobj.get_transform("detector", "v2v3")
    shifts = transform[:4]
    poly = transform[5:7]
    return shifts | Mapping((0, 1, 0, 1)) | poly

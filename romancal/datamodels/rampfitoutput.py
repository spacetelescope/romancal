from .core import RomanDataModel

__all__ = ["RampFitOutputModel"]

class RampFitOutputModel(RomanDataModel):
    """
    A data model for the optional output of the ramp fitting step.

    In the parameter definitions below, `max_seg` is the maximum number
    of segments that were fit, and `ny` and `nx` are the height and width
    of the image.

    Parameters
    __________

    slope : numpy float32 array (max_seg, ny, nx)
        Segment-specific slope

    sigslope : numpy float32 array (max_seg, ny, nx)
        Sigma for segment-specific slope

    var_poisson : numpy float32 array (max_seg, ny, nx)
        Variance due to poisson noise for segment-specific slope

    var_rnoise : numpy float32 array (max_seg, ny, nx)
        Variance due to read noise for segment-specific slope

    yint : numpy float32 array (max_seg, ny, nx)
        Segment-specific y-intercept

    sigyint : numpy float32 array (max_seg, ny, nx)
        Sigma for segment-specific y-intercept

    pedestal : numpy float32 array (max_seg, ny, nx)
        Pedestal array

    weights : numpy float32 array (max_seg, ny, nx)
        Weights for segment-specific fits

    crmag : numpy float32 array (max_seg, ny, nx)
        Approximate CR magnitudes
    """
    schema_url = "http://stsci.edu/schemas/roman_datamodel/rampfitoutput.schema"

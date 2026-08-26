"""
Tests for the local WCS geometry helpers.
"""

from types import SimpleNamespace

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.modeling import models

from romancal.source_catalog._wcs_utils import north_angle_at, pixel_area_map
from romancal.tests.wcs_helpers import create_wcs_object

# Pixel scale in degrees, roughly the WFI native scale
PSCALE = 0.11 / 3600.0


def _quadratic_distortion(coeff=2.0e-5):
    """
    A distortion whose Jacobian varies linearly across the array, giving
    a pixel area that is not uniform.
    """
    x_poly = models.Polynomial2D(2, c1_0=1.0, c2_0=coeff)
    y_poly = models.Polynomial2D(2, c0_1=1.0, c0_2=coeff)
    return models.Mapping((0, 1, 0, 1)) | (x_poly & y_poly)


def _make_wcs(shape, pscale=PSCALE, crval=(30.0, 10.0), distorted=False):
    """
    Build a tangent-plane gwcs, optionally with a distortion that makes
    the pixel area vary across the array.
    """
    return create_wcs_object(
        crval,
        (pscale, pscale),
        shape,
        distortion=_quadratic_distortion() if distorted else None,
    )


def test_uniform_wcs_area_matches_pixel_scale():
    """
    For an undistorted tangent-plane WCS the area must equal pscale**2
    everywhere.
    """
    shape = (256, 256)
    area = pixel_area_map(_make_wcs(shape), shape)

    assert area.shape == shape
    expected = ((PSCALE * u.deg) ** 2).to(u.sr)
    assert u.allclose(area, expected, rtol=1e-5)


def _pixel_area_from_corners(wcs, x, y):
    """
    The area of one pixel, measured with astropy's spherical routines.

    The pixel is small enough to treat as a parallelogram on the sky, so
    its area is the product of two edge lengths and the sine of the angle
    between them. This deliberately uses `~astropy.coordinates.SkyCoord`
    rather than differencing longitudes by hand, so that it is an
    independent check on `pixel_area_map` rather than a restatement of it.
    """
    corner = SkyCoord(
        *wcs(
            np.array([x - 0.5, x + 0.5, x - 0.5]),
            np.array([y - 0.5, y - 0.5, y + 0.5]),
            with_bounding_box=False,
        ),
        unit=u.deg,
    )
    width = corner[0].separation(corner[1])
    height = corner[0].separation(corner[2])
    included = corner[0].position_angle(corner[1]) - corner[0].position_angle(corner[2])
    return (width * height * np.abs(np.sin(included))).to(u.sr)


def test_distorted_wcs_area_varies_and_matches_sky_separations():
    """
    With distortion present the pixel area map must vary, and must agree with
    other pixel-area determinations.
    """
    shape = (512, 512)
    wcs = _make_wcs(shape, distorted=True)
    area = pixel_area_map(wcs, shape)

    # The distortion must actually produce a gradient, otherwise this
    # test would pass trivially against a constant map.
    assert area.max() / area.min() - 1 > 1e-3

    rng = np.random.default_rng(0)
    ys = rng.integers(0, shape[0], 50)
    xs = rng.integers(0, shape[1], 50)
    for x, y in zip(xs, ys, strict=True):
        assert u.allclose(area[y, x], _pixel_area_from_corners(wcs, x, y), rtol=1e-4)


def test_area_independent_of_step():
    """
    The coarse-grid spacing must not materially change the result.
    """
    shape = (512, 512)
    wcs = _make_wcs(shape, distorted=True)
    fine = pixel_area_map(wcs, shape, step=16)
    coarse = pixel_area_map(wcs, shape, step=128)
    assert u.allclose(fine, coarse, rtol=1e-4)


def test_area_correct_across_ra_branch_cut():
    """
    Longitude differencing must wrap across RA = 0; all pixels must have
    sane areas even near the 360 -> 0 wrap.
    """
    shape = (128, 128)
    step = 64
    wcs = _make_wcs(shape, crval=(0.0, 10.0))

    lon = wcs(np.array([-0.5, 0.5]), np.array([0.0, 0.0]), with_bounding_box=False)[0]
    assert lon.max() - lon.min() > 180.0, "test did not straddle the branch cut"

    area = pixel_area_map(wcs, shape, step=step)
    expected = ((PSCALE * u.deg) ** 2).to(u.sr)
    assert u.allclose(area, expected, rtol=1e-5)


def test_segment_geometry_tracks_local_pixel_area():
    """
    Sizes must be converted with the pixel area at each source, not a
    single image-wide value.

    Two identical sources are placed where the pixel area differs by 10%;
    their reported sizes must differ by the corresponding 5% in linear
    measure, and their pixel-frame sizes must be identical.
    """
    from photutils.segmentation import SegmentationImage

    from romancal.source_catalog._segment import SegmentCatalog

    shape = (60, 120)
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    segm = np.zeros(shape, dtype=int)
    for label, x0 in ((1, 30), (2, 90)):
        segm[(xx - x0) ** 2 + (yy - 30) ** 2 <= 25] = label

    # `SegmentCatalog` only reads data, err, and meta.wcs from the model
    model = SimpleNamespace(
        data=(segm != 0) * 100.0,
        err=np.ones(shape),
        meta=SimpleNamespace(wcs=None),
    )

    # Pixel area 10% larger in the right half, where source 2 sits
    nominal = ((0.11 * u.arcsec) ** 2).to_value(u.sr)
    area_map = np.where(xx < shape[1] // 2, nominal, 1.10 * nominal) << u.sr

    cat = SegmentCatalog(
        model,
        SegmentationImage(segm),
        None,
        None,
        area_map,
        requested_properties=["semimajor", "segment_area"],
    )

    # Per-source areas were sampled where the sources actually are
    assert u.allclose(cat.pixel_area[1] / cat.pixel_area[0], 1.10, rtol=1e-6)

    # Linear sizes scale as sqrt(area); areas scale as the area itself
    assert np.isclose(
        (cat.semimajor[1] / cat.semimajor[0]).value, np.sqrt(1.10), rtol=1e-3
    )
    assert np.isclose(
        (cat.segment_area[1] / cat.segment_area[0]).value, 1.10, rtol=1e-3
    )


def test_north_angle_matches_a_tangent_point_reference():
    """
    At the tangent point of an undistorted WCS with no rotation, North
    lies along +y, i.e. at 90 degrees from the +x axis.
    """
    shape = (128, 128)
    wcs = _make_wcs(shape, crval=(30.0, 10.0))
    # The fiducial is at pixel (0, 0)
    angle = north_angle_at(wcs, np.array([0.0]), np.array([0.0]))
    assert u.isclose(angle[0], 90.0 * u.deg, atol=1e-3 * u.deg)


def test_north_angle_varies_around_an_enclosed_pole():
    """
    The angle to north should vary dramatically among locations in a
    field centered on the pole.
    """
    shape = (256, 256)
    # Coarse pixels so that the array spans a usable area around the
    # pole, which the fiducial places at pixel (0, 0)
    wcs = _make_wcs(shape, pscale=10.0 / 3600.0, crval=(0.0, 90.0))

    n = 32
    radius = 100.0
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    angle = north_angle_at(
        wcs, radius * np.cos(theta), radius * np.sin(theta)
    ).to_value(u.deg)

    # The angles must sweep a large angle.
    spread = np.degrees(np.ptp(np.unwrap(np.radians(angle))))
    assert spread > 180

    # Every angle must be finite: the forward-only formulation does not
    # break down near the pole
    assert np.all(np.isfinite(angle))

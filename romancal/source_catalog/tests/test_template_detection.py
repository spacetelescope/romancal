"""
Tests for the matched-filter template detection.
"""

import numpy as np
import pytest
import scipy.ndimage
from astropy.modeling.models import Gaussian2D

from romancal.source_catalog._background import RomanBackground
from romancal.source_catalog._detection import (
    ivw_convolve,
    make_gaussian_kernel,
    snr_from_ivw,
)
from romancal.source_catalog._template_detection import (
    _TEMPLATE_SIZE_FACTOR,
    _bkg_box_size,
    _template_fwhms,
    make_segmentation_image_template,
)


@pytest.fixture
def wide_image():
    """A frame with sources spanning the template sizes, plus a close pair."""
    shape = (330, 330)
    yy, xx = np.mgrid[0 : shape[0], 0 : shape[1]]
    sources = (
        Gaussian2D(300.0, 70, 70, 0.85, 0.85),  # point source
        Gaussian2D(60.0, 170, 70, 3.0, 3.0),  # small galaxy
        Gaussian2D(12.0, 70, 190, 12.0, 12.0),  # large galaxy
        Gaussian2D(200.0, 250, 240, 0.85, 0.85),  # close pair, bright
        Gaussian2D(200.0, 266, 240, 0.85, 0.85),  # close pair, companion
    )
    data = sum(src(xx, yy) for src in sources).astype("float32")
    rng = np.random.default_rng(seed=42)
    data += rng.normal(0.0, 1.0, size=shape)
    err = np.ones(shape, dtype=np.float32)
    return data.astype("float32"), err


def _detect(data, err, mask=None, **kwargs):
    bkg = RomanBackground(data, box_size=100, mask=mask)
    params = dict(
        snr_threshold=5.0,
        n_pixels=9,
        kernel_fwhm=2.0,
        deblend=True,
        mask=mask,
        bkg_boxsize=100,
    )
    params.update(kwargs)
    return make_segmentation_image_template(data - bkg.background, err, **params)


def test_background_box_tracks_the_template():
    """Each template removes a background at about its own kernel scale.

    ``RomanBackground`` medians the mesh again over 3x3, so the scale actually
    removed is about three boxes; the factor of four here is what makes that
    ~12 FWHM.
    """
    boxes = [_bkg_box_size(f) for f in _template_fwhms(2.0)]
    assert boxes == sorted(boxes)  # wider templates, wider background
    for fwhm, box in zip(_template_fwhms(2.0), boxes, strict=True):
        assert 3 * box == pytest.approx(12 * fwhm, rel=0.05)
        # never larger than the kernel, which must already fit the image
        assert (
            box
            <= make_gaussian_kernel(fwhm, size_factor=_TEMPLATE_SIZE_FACTOR).shape[0]
        )


def test_kernels_only_reach_the_wings():
    """Kernels are sized for the template, not for a background box."""
    for fwhm in _template_fwhms(2.0):
        kernel = make_gaussian_kernel(fwhm, size_factor=_TEMPLATE_SIZE_FACTOR)
        assert kernel.shape[0] <= 4 * fwhm + 2
        assert kernel.sum() == pytest.approx(1.0, abs=1e-6)


def test_significance_defined_on_masked_pixels():
    """The matched filter is well defined in small masked regions."""
    data = np.zeros((60, 60), dtype=np.float32)
    err = np.ones_like(data)
    mask = np.zeros(data.shape, dtype=bool)
    mask[30, 30] = True
    wht = np.where(mask, 0.0, 1.0 / err**2)
    num, denom2 = ivw_convolve(data, wht, make_gaussian_kernel(2.0), mask=mask)
    assert denom2[30, 30] > 0
    assert np.isfinite(snr_from_ivw(num, denom2)[30, 30])


def test_templates_larger_than_image_are_skipped(wide_image):
    """Too-large kernels for an image are skipped."""
    data, err = wide_image
    small = data[:120, :120].copy(), err[:120, :120].copy()
    segment_img, _, template_index, _ = _detect(*small)
    if segment_img is not None:
        # only the PSF template fits inside a 120 px frame
        assert np.all(template_index == 0)


def test_detects_sources_of_each_size(wide_image):
    """Sources are assigned to the template that fits them."""
    data, err = wide_image
    segment_img, conv_psf, template_index, significance = _detect(data, err)

    assert segment_img is not None
    assert conv_psf.shape == data.shape
    assert len(template_index) == segment_img.n_labels
    assert len(significance) == segment_img.n_labels
    assert np.all(significance >= 5.0)  # the threshold is the floor
    assert len(_template_fwhms(2.0)) == 3

    # the point source and the large galaxy are not assigned the same template
    labels = np.asarray(segment_img.data)
    point = template_index[labels[70, 70] - 1]
    extended = template_index[labels[190, 70] - 1]
    assert point < extended


def test_segments_are_contiguous(wide_image):
    """Every label is a single connected piece.

    Not automatic: narrowing a max-image segment to its template's footprint, and
    losing pixels to a brighter neighbour, can both split it into islands,
    and a label whose pixels lie in two places has its centroid in the gap
    between them.  `_assign_segments` trims each source to the piece holding
    its peak, and this guards that.
    """
    data, err = wide_image
    segment_img, _, _, _ = _detect(data, err)
    labels = np.asarray(segment_img.data)
    structure = np.ones((3, 3), dtype=bool)
    for label in np.unique(labels[labels > 0]):
        _, n_pieces = scipy.ndimage.label(labels == label, structure=structure)
        assert n_pieces == 1


def test_close_pair_is_deblended(wide_image):
    """Two point sources 16 px apart become two sources."""
    data, err = wide_image
    segment_img, _, _, _ = _detect(data, err)
    labels = np.asarray(segment_img.data)
    assert labels[240, 250] != labels[240, 266]
    assert labels[240, 250] > 0 and labels[240, 266] > 0


def test_masked_stripe_does_not_split_a_segment(wide_image):
    """A masked stripe through a source -- a bleed trail, say -- must not cut
    its segment in two.  The significance and the detection threshold are both
    defined on masked pixels, so the segment closes over the stripe."""
    data, err = wide_image
    mask = np.zeros(data.shape, dtype=bool)
    mask[185:196, 70] = True  # a stripe across the large galaxy

    segment_img, _, _, _ = _detect(data, err, mask=mask)
    assert segment_img is not None
    labels = np.asarray(segment_img.data)

    # the masked pixels themselves are not claimed by any source ...
    assert labels[188, 70] == 0
    # ... but the source is one label, not two, on either side of the stripe
    above, below = labels[182, 70], labels[199, 70]
    assert above > 0
    assert above == below


@pytest.mark.parametrize("spoil", [lambda c: -np.abs(c), lambda c: c * np.inf])
def test_size_cut_counts_moment_positive_pixels(wide_image, monkeypatch, spoil):
    """The minimum segment size counts the pixels the moments can weight.

    Photutils zeroes both the non-positive and the non-finite pixels before
    taking moments, so a segment holding only those has no centroid and returns
    NaN, which later stops the step at the PSF fit.  Whether the moment image
    is negative everywhere or infinite everywhere, no source may survive with
    an undefined centroid.
    """
    import romancal.source_catalog._template_detection as td

    real = td.make_template_snr_images

    def spoiled(*args, **kwargs):
        snr_images, conv = real(*args, **kwargs)
        return snr_images, spoil(conv)

    monkeypatch.setattr(td, "make_template_snr_images", spoiled)
    segment_img, _, _, _ = _detect(*wide_image)
    assert segment_img is None


def test_no_template_fits_returns_no_sources():
    """An image smaller than the narrowest kernel yields no sources, rather
    than failing on a detection image that was never built."""
    data = np.zeros((6, 6), dtype="float32")
    err = np.ones((6, 6), dtype="float32")
    result = make_segmentation_image_template(
        data, err, snr_threshold=5.0, n_pixels=9, kernel_fwhm=2.0, deblend=True
    )
    assert result == (None, None, None, None)


@pytest.mark.parametrize(
    ("depth", "pixel_err"),
    [
        (-1000.0, 1.0),  # a deeply negative pixel with an honest error bar
        (-60.0, 0.25),  # milder, but with the error under-estimated 4x
    ],
)
def test_one_bad_pixel_does_not_make_a_swarm(wide_image, depth, pixel_err):
    """A single badly negative pedestal must not bias the SNR image
    positive.

    Negative pixels can lead to negative sky subtraction leading to
    positive SNR; this test verifies that the local background subtraction
    is robust.
    """
    data, err = wide_image
    data, err = data.copy(), err.copy()
    y, x = 150, 260  # empty sky, far from every source in the fixture
    data[y, x] = depth
    err[y, x] = pixel_err

    segment_img, _, _, _ = _detect(data, err)
    labels = (
        np.zeros(data.shape, dtype=int)
        if segment_img is None
        else np.asarray(segment_img.data)
    )
    box = labels[y - 20 : y + 21, x - 20 : x + 21]
    assert len(np.unique(box[box > 0])) == 0

    # and the real sources are all still there
    for sy, sx in [(70, 70), (70, 170), (190, 70), (240, 250), (240, 266)]:
        assert labels[sy, sx] > 0

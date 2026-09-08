"""
Module to detect sources with a bank of matched filters.

The image is convolved with several templates spanning a range of source
sizes.  Each convolution is turned into a significance image, the
maximum-likelihood amplitude of that template divided by its uncertainty,
and the per-pixel maximum over templates becomes a single detection image on
which sources are found via watershed and deblended.  This enables
the deblending and sounce detection to happen on a single image, reducing
challenges surrounding multiply detecting the same sources with different kernels.

However, for segment definition purposes we want to compute segments going
out to some rough isophotal limit that's appropriate for the source in question;
we do not want to convolve a bright star with a large galaxy template and
construct a segment that goes out until that large galaxy cross-correlation
has fallen off significantly.  Because each source's peak in the significance
image corresponds to a peak from a specific template, we can use the max
significance image to define a parent set of segments and then intersect
that with the set of significant pixels belonging to the template that
generated the peak.  So PSFs end up going out to isophotes appropriate
for PSFs and large galaxies go deeper.

Additionally, before each convolution we subtract a local background
estimated at a scale set by that template, so that a small kernel does not
fire merely because it rides on the light of a large source.  The background
is a sigma-clipped median of the pixels rather than a mean, so that a single
corrupt pixel cannot move it.
"""

import logging
import math

import numpy as np
import scipy.ndimage
from photutils.segmentation import SegmentationImage

from romancal.source_catalog._background import RomanBackground
from romancal.source_catalog._detection import (
    fft_convolve,
    ivw_convolve,
    make_gaussian_kernel,
    make_segmentation_image,
    snr_from_ivw,
)

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Template bank.  Every template is a Gaussian; these are the FWHM of the
# larger ones, standing in for galaxies.  The PSF template comes first and
# is sized by the step's ``kernel_fwhm``, so the default bank is 0.2, 0.6
# and 2.4 arcsec.  The sizes are angular rather than in pixels so that the
# same physical scales are searched whatever the pixel scale: an L3 coadd
# at 0.055 arcsec/pixel spreads a galaxy over twice as many pixels as the
# L2 image at 0.11 that went into it.
#
# The spacing is logarithmic.  A source falling between two rungs is worst
# off at their geometric mean, where a template mismatched in width by a
# factor r keeps 2r / (1 + r^2) of the SNR: 87 percent across the first
# gap of three, 80 percent across the second of four.  The tighter gap is
# at the bottom because that is where the small faint galaxies are.
_TEMPLATE_FWHM = (0.6, 2.4)  # arcsec

# Kernel box size for the templates, in units of the template FWHM.  The
# kernel only has to reach the template's wings, so a few FWHM is enough.
_TEMPLATE_SIZE_FACTOR = 4

# Local background scale for each template, as a multiple of its FWHM.
# ``RomanBackground`` medians the pixels in boxes of this size and then
# medians that mesh again over a 3x3 window, so the scale it actually
# removes is about three boxes, a scale of _BKG_BOX_FACTOR * 3 FWHM,
# ~ 12 FWHM with the default 4.
_BKG_BOX_FACTOR = 4

# Kernel box size for the image the moments are measured on, in units of the
# PSF FWHM.  Nothing is being subtracted here, so the box need only enclose
# most of the flux.
_MOMENT_SIZE_FACTOR = 4

# Deblending contrast for the maximum image; children must hold at least this
# fraction of a parent's flux to be deblended.
_DEBLEND_CONTRAST = 0.001

# Grow each final segment by this many pixels.  We are willing to detect
# sources that are 5 sigma in a single pixel, but want to then grow those segments
# to be at least large enough to compute a moment.  We dilate so that isolated
# single pixels grow to 3x3 segments.
_SEGMENT_DILATE = 1


def _template_fwhms(kernel_fwhm, pixel_scale):
    """Template FWHMs in pixels, PSF first.

    ``kernel_fwhm`` and the bank are in arcsec; ``pixel_scale`` is in
    arcsec per pixel.
    """
    return tuple(f / pixel_scale for f in (float(kernel_fwhm), *_TEMPLATE_FWHM))


def _bkg_box_size(fwhm):
    """Background mesh box size, in pixels, for a template of this FWHM."""
    return max(math.ceil(_BKG_BOX_FACTOR * fwhm), 1)


def make_template_snr_images(data, err, kernel_fwhm, pixel_scale, mask=None):
    """
    Compute a matched-filter significance image for each template.

    Parameters
    ----------
    data : 2D `numpy.ndarray`
        Background-subtracted data.
    err : 2D `numpy.ndarray`
        Per-pixel uncertainty.
    kernel_fwhm : float
        FWHM of the PSF template, in arcsec.
    pixel_scale : float
        Pixel scale in arcsec per pixel, used to put the template bank on
        the pixel grid.
    mask : 2D `numpy.ndarray`, optional
        Boolean mask; True values are given zero weight.

    Returns
    -------
    snr_images : list of 2D `numpy.ndarray`
        One matched-filter SNR image per template, each computed after
        subtracting a local background at that template's own scale.
    conv_psf : 2D `numpy.ndarray`
        The PSF-convolved flux image, for centroids and moments.  No
        background is subtracted from this one; see the note at the call
        site.
    """
    # we use single precision here because these arrays and the transforms built
    # from them are the bulk of the step's memory.
    data = np.asarray(data, dtype=np.float32)
    good = ~mask if mask is not None else np.ones(data.shape, dtype=bool)
    wht = np.where(good, 1.0 / np.where(good, err, 1.0) ** 2, 0.0)
    wht = wht.astype(np.float32)

    snr_images = []
    for fwhm in _template_fwhms(kernel_fwhm, pixel_scale):
        kernel = make_gaussian_kernel(fwhm, size_factor=_TEMPLATE_SIZE_FACTOR)
        if kernel.shape[0] > min(data.shape):
            # The kernel is wider than the image, so it is all boundary.
            # Triggered by cutouts and by test frames; real images are far
            # larger than the bank.
            log.info(
                f"Skipping template FWHM={fwhm:.1f} px: kernel "
                f"{kernel.shape[0]} px exceeds the image"
            )
            continue
        # Remove the light this template should not be detecting: anything
        # smooth on a scale well outside the template itself.  Done on the
        # pixels rather than on the convolution, because convolving first
        # would smear a single bad pixel across the whole kernel footprint
        # and leave the median nothing to reject.
        box = _bkg_box_size(fwhm)
        bkg = RomanBackground(data, box_size=box, mask=mask).background
        num, denom2 = ivw_convolve(data - bkg, wht, kernel, mask=mask)
        snr_images.append(snr_from_ivw(num, denom2))
        del bkg
        log.info(
            f"Template FWHM={fwhm:.1f} px: kernel {kernel.shape[0]} px, "
            f"background box {box} px"
        )

    if not snr_images:
        log.warning(
            f"No template fits within an image of shape {data.shape}; "
            "no sources can be detected"
        )
        return [], None

    # Centroids and moments are measured on a PSF-convolved image with no
    # additional background removed.
    # We should investigate subtracting one, but we need to more carefully
    # track what scale the background should be removed on for each different template
    psf_kernel = make_gaussian_kernel(
        kernel_fwhm / pixel_scale, size_factor=_MOMENT_SIZE_FACTOR
    )
    conv_psf = fft_convolve(data, psf_kernel, mask=mask)

    return snr_images, conv_psf


def _max_detection_image(snr_images):
    """
    Per-pixel maximum over the templates' significance images.

    Each input is already in units of its own background noise.

    Parameters
    ----------
    snr_images : list of 2D `numpy.ndarray`
        One significance image per template, in sigma.

    Returns
    -------
    max_image : 2D `numpy.ndarray`
        The per-pixel maximum, in sigma.
    template_at_pixel : 2D `numpy.ndarray`
        Index of the template attaining the maximum at each pixel.
    """
    max_image = None
    template_at_pixel = None
    for index, snr in enumerate(snr_images):
        if max_image is None:
            max_image = snr.copy()
            template_at_pixel = np.zeros(snr.shape, dtype=np.uint8)
        else:
            better = snr > max_image
            max_image[better] = snr[better]
            template_at_pixel[better] = index
    return max_image, template_at_pixel


def _assign_segments(
    deblended,
    max_image,
    template_at_pixel,
    footprints,
    moment_image,
    n_pixels,
    mask=None,
    dilate=_SEGMENT_DILATE,
    max_sources=0,
):
    """
    Turn deblended segments into the final segmentation image.

    Start with a set of segments derived from the max_image, and refine
    those by looping over segments in order of their peak SNR.  For each
    segment, find the template that corresponds to it, intersect the
    max_image segment with the significant regions from the corresponding
    template, essentially narrowing the segment to what would have been found
    using that template only.  Paint the intersection onto the final derived
    segmentation image, after dilating.  Later lower SNR segments
    cannot overwrite pixels claimed by an earlier segment, and are dropped
    if they end up with fewer than ``n_pixels`` pixels the photometry can
    weight.

    Parameters
    ----------
    deblended : `SegmentationImage`
        The deblended segments from the max_image, one label per source,
        covering the whole detection footprint.
    max_image : 2D `numpy.ndarray`
        The maximum-over-templates significance image, in sigma.  Used to
        locate each segment's peak and to read its significance there.
    template_at_pixel : 2D `numpy.ndarray`
        Index of the template attaining the maximum at each pixel.
    footprints : list of 2D `numpy.ndarray`
        Boolean footprint of each template, i.e. where that template's SNR
        exceeds its own threshold.  Used to narrow segments from max_det to
        the segments from the template that fits it.
    moment_image : 2D `numpy.ndarray`
        The image the moments will be measured on.  Only its sign is used, to
        count the pixels a segment contributes to its own centroid.
    n_pixels : int
        Smallest final segment to keep, counting only positive, unmasked pixels
        in the moment image, so moments can be computed.
    mask : 2D `numpy.ndarray`, optional
        Boolean mask indicating bad pixels.
    dilate : int, optional
        Grow each segment by this many pixels before painting.
    max_sources : int, optional
        Stop after this many sources have been painted.  Zero means no limit.

    Returns
    -------
    segment_img : `SegmentationImage` or None
    template_index : 1D `numpy.ndarray`
    significance : 1D `numpy.ndarray`

    Returns
    -------
    segment_img : `SegmentationImage` or None
    template_index : 1D `numpy.ndarray`
    significance : 1D `numpy.ndarray`
    """
    labels = np.asarray(deblended.data)
    shape = labels.shape
    structure = np.ones((3, 3), dtype=bool)

    # for a pixel to contribute to a moment it needs to be in a segment
    # with a finite, positive value; select these
    moment_positive = np.isfinite(moment_image) & (moment_image > 0)

    # Collect each segment's peak and the template that wins there.
    candidates = []
    for label, slices in enumerate(scipy.ndimage.find_objects(labels), start=1):
        if slices is None:
            continue
        inside = labels[slices] == label
        peak = np.unravel_index(
            np.argmax(np.where(inside, max_image[slices], -np.inf)), inside.shape
        )
        ypeak = peak[0] + slices[0].start
        xpeak = peak[1] + slices[1].start
        index = int(template_at_pixel[ypeak, xpeak])
        candidates.append(
            (float(max_image[ypeak, xpeak]), slices, inside, ypeak, xpeak, index)
        )

    # brightest first
    candidates.sort(key=lambda item: -item[0])

    merged = np.zeros(shape, dtype=np.int32)
    # 'claimed' tracks every pixel a source has taken, independent of
    # the mask.  It differs from merged in that merged contains only
    # the unmasked pixels.
    claimed = np.zeros(shape, dtype=bool)
    template_index = []
    significance = []
    n_label = 0
    n_dropped = 0
    n_trimmed = 0

    n_capped = 0
    for sig, slices, inside, ypeak, xpeak, index in candidates:
        if max_sources and n_label >= max_sources:
            # Allow early termination for crowded fields if requested, keeping
            # the bright sources.
            n_capped = len(candidates) - n_label - n_dropped
            # n_capped = the number of sources skipped due to the cap
            break

        # pad to allow segment dilation
        pad = int(dilate) + 1
        y0 = max(slices[0].start - pad, 0)
        y1 = min(slices[0].stop + pad, shape[0])
        x0 = max(slices[1].start - pad, 0)
        x1 = min(slices[1].stop + pad, shape[1])
        window = (slice(y0, y1), slice(x0, x1))

        # construct nominal segment in padded stamp
        local = np.zeros((y1 - y0, x1 - x0), dtype=bool)
        local[
            slices[0].start - y0 : slices[0].stop - y0,
            slices[1].start - x0 : slices[1].stop - x0,
        ] = inside

        # Narrow to the isophote of the template that fits this source, so a
        # star keeps the compact PSF isophote rather than the inflated one from
        # max_image.
        # The peak significance is always inside its template footprint,
        # since being the maximum in the segment defines the winning template
        local &= footprints[index][window]

        # dilate local segment
        if dilate > 0:
            local = scipy.ndimage.binary_dilation(
                local, structure=structure, iterations=int(dilate)
            )
        # restrict to pixels that haven't already been claimed
        local &= ~claimed[window]

        if local.any():
            # intersecting with a template footprint or losing pixels to
            # a brighter neighbour can split a segment into islands.
            # Here we trim the segment to the island containing the peak.
            pieces, n_pieces = scipy.ndimage.label(local, structure=structure)
            if n_pieces > 1:
                keep = int(pieces[ypeak - y0, xpeak - x0])
                if keep == 0:  # peak taken by a brighter neighbour
                    counts = np.bincount(pieces.ravel())
                    counts[0] = 0
                    keep = int(counts.argmax())
                # trimming to the segment containing the peak
                local = pieces == keep
                n_trimmed += 1

        # add the bad pixel mask, which had been left off intentionally
        # so that masked pixels could not split segments into islands
        usable = local if mask is None else local & ~mask[window]
        # count only usable pixels that can be used for computing moments
        # so that moment computation cannot fail
        if int((usable & moment_positive[window]).sum()) < n_pixels:
            n_dropped += 1
            continue

        n_label += 1
        claimed[window] |= local
        merged[window][usable] = n_label
        template_index.append(index)
        significance.append(sig)

    log.info(
        f"Detected {n_label} sources ({n_dropped} dropped with fewer than "
        f"{n_pixels} usable positive pixels, {n_trimmed} trimmed "
        "to their peak component)"
    )
    if n_capped:
        log.info(
            f"Kept the {max_sources} most significant sources; skipped "
            f"{n_capped} fainter candidates (max_sources)"
        )
    if n_label == 0:
        return None, None, None

    return (
        SegmentationImage(merged),
        np.array(template_index, dtype=np.uint8),
        np.array(significance, dtype=np.float32),
    )


def make_segmentation_image_template(
    data,
    err,
    snr_threshold,
    n_pixels,
    kernel_fwhm,
    pixel_scale,
    deblend=True,
    mask=None,
    bkg_boxsize=100,
    max_sources=0,
):
    """
    Detect sources with a bank of matched filters.

    Parameters
    ----------
    data : 2D `numpy.ndarray`
        Background-subtracted data.
    err : 2D `numpy.ndarray`
        Per-pixel uncertainty.
    snr_threshold : float
        Detection threshold in sigma, the same for every template.
    n_pixels : int
        Smallest final segment to keep, counted in the pixels the photometry
        can use for moment computation.  It is applied after narrowing and
        dilation, not to the deblended children of the maximum image.
    kernel_fwhm : float
        FWHM of the PSF template, in arcsec.
    pixel_scale : float
        Pixel scale in arcsec per pixel.
    deblend : bool, optional
        Whether to deblend the maximum image.
    mask : 2D `numpy.ndarray`, optional
        Boolean mask; True values are given zero weight.
    bkg_boxsize : int, optional
        Box size for estimating the noise of each significance image.
    max_sources : int, optional
        Keep at most this many sources, the most significant first.  Zero
        means no limit.

    Returns
    -------
    segment_img : `SegmentationImage` or None
        The segmentation image.
    conv_psf : 2D `numpy.ndarray`
        The PSF-convolved flux image, for centroids and moments.
    template_index : 1D `numpy.ndarray` or None
        Index of the template that detected each source.
    significance : 1D `numpy.ndarray` or None
        Peak significance of each source, in sigma.
    """
    snr_images, conv_psf = make_template_snr_images(
        data, err, kernel_fwhm, pixel_scale, mask=mask
    )
    if not snr_images:
        # Nothing in the bank fits the image; already logged.
        return None, None, None, None

    # Put each template's SNR image in units of its own background noise, so
    # that one threshold means the same thing for every template and the
    # maximum over templates is itself a significance image in sigma.  The
    # noise is not 1 in practice: signal leaks into the background estimate,
    # and resampling correlates the noise, both more so for wider kernels.
    footprints = []
    for i, snr in enumerate(snr_images):
        # mask rather than coverage_mask because snr_image is well defined
        # on masked pixels, and using mask leads bkg_rms to be well defined
        # there
        bkg_rms = RomanBackground(snr, box_size=bkg_boxsize, mask=mask).background_rms
        snr_images[i] = np.where(bkg_rms > 0, snr / bkg_rms, 0.0)
        footprints.append(snr_images[i] > snr_threshold)

    max_image, template_at_pixel = _max_detection_image(snr_images)
    max_image = np.where(np.isfinite(max_image), max_image, 0.0)

    # Detection and deblending on the maximum significance image
    # A unit background RMS is passed because the image has already
    # been divided by the estimated noise.
    # We here allow single-pixel children but impose a n_pixels cut at the end,
    # once each segment has been narrowed to its template's isophote
    # and dilated.  n_pixels = 1 corresponds to a true 5-sigma source in
    # these detections, and so is a natural value.
    deblended = make_segmentation_image(
        max_image,
        snr_threshold,
        1,
        1.0,
        deblend=deblend,
        mask=None,
        deblend_contrast=_DEBLEND_CONTRAST,
    )
    if deblended is None:
        return None, conv_psf, None, None
    log.info(f"{deblended.n_labels} sources after deblending")

    segment_img, template_index, significance = _assign_segments(
        deblended,
        max_image,
        template_at_pixel,
        footprints,
        conv_psf,
        n_pixels,
        mask=mask,
        max_sources=max_sources,
    )

    return segment_img, conv_psf, template_index, significance

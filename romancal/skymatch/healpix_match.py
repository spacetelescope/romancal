"""
HEALPix-based multi-exposure sky level solver.

Romancal-native port of the pairwise inverse-variance least-squares approach
used by ``wfiskymatch.pixmatch``, operating on in-memory skyval samples.

Given overlapping HEALPix sky measurements from multiple exposures, solve for
one sky level per exposure. Downstream outlier/resample **subtract**
``meta.background.level``, so returned values are those levels.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence

import numpy as np
from astropy.stats import biweight_location
from scipy.sparse import csc_array

log = logging.getLogger(__name__)

__all__ = ["match_healpix_skyvals"]

# Weak prior: level ~ 0 +/- PRIOR_SIGMA_DN (before photom scaling).
PRIOR_SIGMA_DN = 10.0


def match_healpix_skyvals(
    samples: Sequence[Mapping[str, np.ndarray]],
    *,
    names: Sequence[str] | None = None,
    min_overlap_cells: int = 1,
    recenter: bool = True,
    prior_sigma_dn: float = PRIOR_SIGMA_DN,
) -> tuple[np.ndarray, list[str]]:
    """
    Solve per-exposure sky levels from HEALPix skyval samples.

    Parameters
    ----------
    samples
        One mapping per exposure with keys:

        * ``pixel`` : (N,) int — nested HEALPix indices (e.g. nside=2**17)
        * ``value`` : (N,) float — sky estimates (DN/s unless already scaled)
        * ``unc`` : (N,) float — uncertainties (must be > 0 for used cells)

        Optional:

        * ``coverage`` : coarse indices used to skip pairs with no coarse
          overlap (performance).
        * ``conversion_megajanskys`` : float — L2 photom factor DN/s → MJy/sr.
          Applied to ``value`` and ``unc`` before solving so returned levels
          match flux-calibrated image units.

    names
        Optional exposure labels (length ``len(samples)``). Defaults to
        ``exp000``, ``exp001``, ...
    min_overlap_cells
        Minimum shared fine healpix cells required to form a pair constraint.
    recenter
        If True (default), subtract biweight location so levels are centered
        (``wfiskymatch.pixmatch`` behavior).
    prior_sigma_dn
        Weak Gaussian prior width in DN/s toward zero level. After photom
        scaling, exposure *i* gets diagonal weight
        ``1/(prior_sigma_dn * conversion_i)**2``. Set ``<= 0`` to disable.

    Returns
    -------
    levels : ndarray of float64, shape (n_exp,)
        Sky levels to store as ``meta.background.level`` (subtract from data).
    names : list of str
        Exposure names aligned with ``levels``.
    """
    nfiles = len(samples)
    if nfiles == 0:
        return np.array([], dtype=np.float64), []
    if names is None:
        names = [f"exp{i:03d}" for i in range(nfiles)]
    if len(names) != nfiles:
        raise ValueError("names must match the number of sample entries")

    pixels = []
    values = []
    unc2 = []
    coverage = []
    prior_weight = np.zeros(nfiles, dtype=np.float64)

    for i, sample in enumerate(samples):
        pix = np.asarray(sample["pixel"], dtype=np.int64).reshape(-1)
        val = np.asarray(sample["value"], dtype=np.float64).reshape(-1)
        unc = np.asarray(sample["unc"], dtype=np.float64).reshape(-1)
        if pix.size != val.size or pix.size != unc.size:
            raise ValueError(
                f"sample {names[i]!r}: pixel/value/unc length mismatch "
                f"({pix.size}, {val.size}, {unc.size})"
            )

        conv = float(sample.get("conversion_megajanskys", 1.0) or 1.0)
        if not np.isfinite(conv) or conv <= 0:
            raise ValueError(
                f"sample {names[i]!r}: conversion_megajanskys must be finite and > 0"
            )
        val = val * conv
        unc = unc * conv
        if prior_sigma_dn and prior_sigma_dn > 0:
            sigma = prior_sigma_dn * conv
            prior_weight[i] = 1.0 / (sigma * sigma)

        order = np.argsort(pix, kind="mergesort")
        pix, val, unc = pix[order], val[order], unc[order]
        bad = ~np.isfinite(unc) | (unc <= 0) | ~np.isfinite(val)
        if np.any(bad):
            pix, val, unc = pix[~bad], val[~bad], unc[~bad]
        pixels.append(pix)
        values.append(val)
        unc2.append(unc * unc)
        if "coverage" in sample and sample["coverage"] is not None:
            coverage.append(np.unique(np.asarray(sample["coverage"], dtype=np.int64)))
        else:
            coverage.append(None)

    xcorr = np.ones((nfiles, nfiles), dtype=bool)
    np.fill_diagonal(xcorr, False)
    if all(c is not None for c in coverage):
        xcorr[:] = False
        for i in range(nfiles - 1):
            for j in range(i + 1, nfiles):
                if np.intersect1d(coverage[i], coverage[j], assume_unique=True).size:
                    xcorr[i, j] = xcorr[j, i] = True

    # Pairwise system in pixmatch "epsilon" gauge: corrected = value + epsilon.
    # Convert to subtractive levels at the end (level = -epsilon).
    row: list[int] = []
    col: list[int] = []
    data: list[float] = []
    B = np.zeros(nfiles, dtype=np.float64)
    D = np.zeros(nfiles, dtype=np.float64)
    n_pairs = 0

    for i in range(nfiles - 1):
        apix, aval, aunc2 = pixels[i], values[i], unc2[i]
        if apix.size == 0:
            continue
        for j in range(i + 1, nfiles):
            if not xcorr[i, j]:
                continue
            bpix, bval, bunc2 = pixels[j], values[j], unc2[j]
            if bpix.size == 0:
                continue
            common, aidx, bidx = np.intersect1d(
                apix, bpix, assume_unique=True, return_indices=True
            )
            if common.size < min_overlap_cells:
                continue
            w = 1.0 / (aunc2[aidx] + bunc2[bidx])
            a_ij = -np.sum(w)
            I_ij = -np.sum((aval[aidx] - bval[bidx]) * w)
            col.extend([i, j])
            row.extend([j, i])
            data.extend([a_ij, a_ij])
            B[i] += I_ij
            B[j] -= I_ij
            D[i] -= a_ij
            D[j] -= a_ij
            n_pairs += 1

    log.info("healpix sky match: %d exposures, %d pairwise overlaps", nfiles, n_pairs)

    if n_pairs == 0 and not np.any(prior_weight > 0):
        log.warning(
            "No overlapping HEALPix skyval pairs; returning zero levels for all exposures"
        )
        return np.zeros(nfiles, dtype=np.float64), list(names)

    ii = np.arange(nfiles)
    row.extend(ii.tolist())
    col.extend(ii.tolist())
    data.extend(D.tolist())
    A = csc_array((data, (row, col)), shape=(nfiles, nfiles))

    A_solve = A.toarray()
    B_solve = B.copy()
    if np.any(prior_weight > 0):
        A_solve = A_solve + np.diag(prior_weight)
    else:
        A_solve[-1, :] = 1.0
        B_solve[-1] = 0.0

    try:
        epsilon = np.linalg.solve(A_solve, B_solve)
    except np.linalg.LinAlgError as err:
        raise RuntimeError(
            "healpix sky match failed: singular solution "
            f"(n_pairs={n_pairs}). Check overlaps and uncertainties."
        ) from err
    epsilon = np.asarray(epsilon, dtype=np.float64)
    if not np.all(np.isfinite(epsilon)):
        raise RuntimeError(
            "healpix sky match failed: non-finite solution "
            f"(n_pairs={n_pairs}). Check overlaps and uncertainties."
        )

    if recenter:
        delta = -float(biweight_location(epsilon))
        epsilon = epsilon + delta
        log.debug("healpix sky match biweight recenter delta=%s", delta)

    levels = -epsilon
    return levels, list(names)

"""
Module to compute nearest neighbors for a catalog of sources.
"""

import numpy as np
from astropy.utils.decorators import lazyproperty
from scipy.spatial import KDTree


class NNCatalog:
    """
    Class to compute nearest neighbors for a catalog of sources.
    """

    # Source-property column names this catalog can produce
    available_properties = ("nn_label", "nn_distance")

    def __init__(
        self, label, xypos, xypos_finite, pixel_scale, *, requested_properties=None
    ):
        self.label = label
        self.xypos = xypos
        self.xypos_finite = xypos_finite
        self.pixel_scale = pixel_scale

        self.nonfinite_mask = ~np.isfinite(xypos).all(axis=1)

        if requested_properties is None:
            self.properties = list(self.available_properties)
        else:
            requested = set(requested_properties)
            self.properties = [
                prop for prop in self.available_properties if prop in requested
            ]

    def __len__(self):
        """
        Return the number of sources in the catalog.
        """
        return len(self.label)

    @lazyproperty
    def _kdtree_query(self):
        """
        The distance in pixels to the nearest neighbor and its index.

        Returns
        -------
        nn_distance : `~numpy.ndarray`
            Length-N array of pixel distances to the nearest neighbor.
            For catalogs with fewer than two sources the distance is
            NaN.

        nn_index : `~numpy.ndarray`
            Length-N array of indices into ``self.label`` for the
            nearest neighbor. For catalogs with fewer than two sources
            the index is 0 (the corresponding label is replaced with -1
            by `nn_label`).
        """
        # Skip the KDTree for degenerate catalogs (0 or 1 sources)
        n_sources = len(self)
        if n_sources <= 1:
            return (
                np.full(n_sources, np.nan, dtype=float),
                np.zeros(n_sources, dtype=np.intp),
            )

        # Non-finite xypos causes memory errors on linux, but not MacOS
        tree = KDTree(self.xypos_finite)
        qdist, qidx = tree.query(self.xypos_finite, k=[2])

        return np.transpose(qdist)[0], np.transpose(qidx)[0]

    @lazyproperty
    def nn_label(self):
        """
        The label number of the nearest neighbor.

        A label value of -1 is returned for catalogs with fewer than two
        sources and for sources with a non-finite centroid.
        """
        if len(self) <= 1:
            return np.full(len(self), -1, dtype=np.int32)

        nn_label = self.label[self._kdtree_query[1]].astype(np.int32)
        # Assign a label of -1 for non-finite xypos
        nn_label[self.nonfinite_mask] = -1

        return nn_label

    @lazyproperty
    def nn_distance(self):
        """
        The distance in arcsec to the nearest neighbor.

        NaN is returned for non-finite centroid positions or for
        catalogs with fewer than two sources.
        """
        nn_distance = self._kdtree_query[0]
        if len(self) > 1:
            # Assign a distance of np.nan for non-finite xypos
            nn_distance[self.nonfinite_mask] = np.nan

        return (nn_distance * self.pixel_scale).astype(np.float32)

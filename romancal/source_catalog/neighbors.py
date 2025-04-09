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

    def __init__(self, label, xypos, xypos_finite, pixel_scale):
        self.label = label
        self.xypos = xypos
        self.xypos_finite = xypos_finite
        self.pixel_scale = pixel_scale

        self.nonfinite_mask = ~np.isfinite(xypos).all(axis=1)

        self.names = ["nn_label", "nn_dist"]

    def __len__(self):
        """
        Return the number of sources in the catalog.
        """
        return len(self.label)

    @lazyproperty
    def _kdtree_query(self):
        """
        The distance in pixels to the nearest neighbor and its index.
        """
        if len(self) == 1:
            return [np.nan], [np.nan]

        # non-finite xypos causes memory errors on linux, but not MacOS
        tree = KDTree(self.xypos_finite)
        qdist, qidx = tree.query(self.xypos_finite, k=[2])
        return np.transpose(qdist)[0], np.transpose(qidx)[0]

    @lazyproperty
    def nn_label(self):
        """
        The label number of the nearest neighbor.

        A label value of -1 is returned if there is only one detected
        source and for sources with a non-finite centroid.
        """
        if len(self) == 1:
            return np.int32(-1)

        nn_label = self.label[self._kdtree_query[1]].astype(np.int32)
        # assign a label of -1 for non-finite xypos
        nn_label[self.nonfinite_mask] = -1

        return nn_label

    @lazyproperty
    def nn_dist(self):
        """
        The distance in arcsec to the nearest neighbor.

        NaN is returned for non-finite centroid positions or when
        the catalog contains only one source.
        """
        nn_dist = self._kdtree_query[0]
        if len(self) != 1:
            # assign a distance of np.nan for non-finite xypos
            nn_dist[self.nonfinite_mask] = np.nan

        return (nn_dist * self.pixel_scale).astype(np.float32)

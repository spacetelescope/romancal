"""
convert old patch names (i.e. `r274dp63x31y80`) from the old patch table
to the names of the closest sky cells (i.e. `270p65x48y69`) from the new sky cell table

remove this script when we no longer have any old files with old patch names
"""

from pathlib import Path

import asdf
import numpy as np
import spherical_geometry.vector as sgv

import romancal.skycell.skymap as sc

PATCH_TABLE_PATH = Path("/grp/roman/scsb/tesselation/patches.asdf")


def patch_names_to_skycell_names(
    patch_names: list[str],
    patch_table_path: Path = PATCH_TABLE_PATH,
    skymap: sc.SkyMap = sc.SKYMAP,
) -> list[str]:
    """
    convert old patch names (i.e. `r274dp63x31y80`) from the old patch table
    to the names of the closest sky cells (i.e. `270p65x48y69`) from the new sky cell table

    Parameters
    ----------
    patch_names : list[str]
        list of patch names, i.e. `r274dp63x31y80`
    patch_table_path : Path
        path to patch table; defaults to /grp/roman/scsb/tesselation/patches.asdf
    skymap : romancal.skycell.skymap.SkyMap
        SkyMap instance; defaults to global `skycells` table from CRDS

    Returns
    -------
    list of sky cell names, i.e. `270p65x48y69`

    Examples
    --------
    >>> patch_names_to_skycell_names(["r274dp63x31y80", "r274dp63x31y81", "r274dp63x32y82", "r274dp63x32y80", "r274dp63x32y81", "r274dp63x32y82"])
    ['225p90x49y67', '225p90x49y69', '225p90x50y70', '225p90x50y67', '225p90x50y69', '225p90x50y70']
    """

    patch_table = asdf.open(patch_table_path, memmap=True)

    skycell_names = []
    projregions = {}
    for patch_name in patch_names:
        patch = patch_table["patches"][
            np.where(patch_table["patches"]["name"] == patch_name)
        ]
        patch_radec_center = np.array(
            sgv.lonlat_to_vector(patch["ra_center"], patch["dec_center"])
        ).T[0]
        neighboring_projregion_indices = skymap.projection_regions_kdtree.query(
            patch_radec_center,
            k=4,
        )[1]

        skycell_indices = {}
        for index in neighboring_projregion_indices:
            if index not in projregions:
                projregions[index] = sc.ProjectionRegion(index)
            distance, skycell_index = projregions[index].skycells_kdtree.query(
                patch_radec_center, k=1
            )
            skycell_indices[distance] = skycell_index

        skycell_names.append(sc.SkyCell(skycell_indices[min(skycell_indices)]).name)

    return np.array(skycell_names).tolist()

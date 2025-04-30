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
    patch_table = asdf.open(patch_table_path)

    skycell_names = []
    projregions = {}
    for patch_name in patch_names:
        patch = patch_table["patches"][
            np.where(patch_table["patches"]["name"] == patch_name)
        ]
        patch_radec_center = np.array(
            sgv.lonlat_to_vector(patch["ra_center"], patch["dec_center"])
        ).T[0]
        neighboring_projregion_indices = sc.SKYMAP.projection_regions_kdtree.query(
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

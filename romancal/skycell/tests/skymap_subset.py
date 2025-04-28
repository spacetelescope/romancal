from pathlib import Path

import romancal.skycell.skymap as sc

LAST_PROJREGION_INDEX = 2

if __name__ == "__main__":
    skymap_subset = sc.SKYMAP.data.copy()

    # to maintain the proper indices, the subset must contain all the previous projection regions up to the specified index
    skymap_subset["roman"]["projection_regions"] = sc.SKYMAP.projection_regions[
        : LAST_PROJREGION_INDEX + 1
    ].copy()
    skymap_subset["roman"]["skycells"] = sc.SKYMAP.skycells[
        : skymap_subset["roman"]["projection_regions"][-1]["skycell_end"] + 1
    ].copy()

    skymap_subset.write_to(Path(__file__).parent / "skymap_subset.asdf")

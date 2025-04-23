from pathlib import Path

import romancal.skycell.skymap as sm

SKYTILES = 23

if __name__ == "__main__":
    skymap_subset = sm.SKYMAP.data.copy()

    # to maintain the proper indices, the subset must contain all the previous skytiles up to the specified index
    skymap_subset["roman"]["projection_regions"] = sm.SKYMAP.skytiles[
        : SKYTILES + 1
    ].copy()
    skymap_subset["roman"]["skycells"] = sm.SKYMAP.skycells[
        : skymap_subset["roman"]["projection_regions"][-1]["skycell_end"] + 1
    ].copy()

    skymap_subset.write_to(Path(__file__).parent / "skymap_subset.asdf")

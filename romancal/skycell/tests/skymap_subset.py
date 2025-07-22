from pathlib import Path
from romancal.skycell import skymap

LAST_PROJREGION_INDEX = 1
#DATA_DIRECTORY = Path(__file__).parent / "data"
DATA_DIRECTORY = Path("/Users/dencheva")

if __name__ == "__main__":
    skymap_subset = skymap.SKYMAP.model.copy()

    skymap_subset.projection_regions = skymap.SKYMAP.model.projection_regions[: LAST_PROJREGION_INDEX + 1].copy()
    skymap_subset.skycells = skymap.SKYMAP.model.skycells[: skymap_subset.projection_regions[-1]["skycell_end"] + 1].copy()
    skymap_subset.save(DATA_DIRECTORY / "skymap_subset.asdf")

from matplotlib import pyplot

import romancal.skycell.plot as sp
import romancal.skycell.skymap as sc
import romancal.skycell.tests.test_skycell_match as tsm

if __name__ == "__main__":
    # corners = tsm.mk_im_corners(
    #     sc.SKYMAP.skycells[3560]["ra_corn3"],
    #     sc.SKYMAP.skycells[3560]["dec_corn3"] + tsm.e,
    #     tsm.cpa,
    #     0.001,
    # )
    # corners = tsm.mk_im_corners(
    #     sc.SKYMAP.skycells[3560]["ra_corn3"],
    #     sc.SKYMAP.skycells[3560]["dec_corn3"],
    #     tsm.cpa,
    #     0.5,
    # )
    corners = tsm.mk_im_corners(
        sc.SKYMAP.skycells[3000]["ra_corn3"],
        sc.SKYMAP.skycells[3000]["dec_corn3"],
        45.0,
        0.5,
    )

    sp.plot_image_footprint_and_skycells(corners)

    pyplot.show()

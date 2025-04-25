from matplotlib import pyplot

import romancal.skycell.plot as sp
import romancal.skycell.skymap as sc
import romancal.skycell.tests.test_skycell_match as tsm

if __name__ == "__main__":
    crecord = sc.SKYMAP.skycells[3560]
    cra = crecord["ra_corn3"]
    cdec = crecord["dec_corn3"]

    # corners = tsm.mk_im_corners(cra, cdec + tsm.e, tsm.cpa, 0.001)
    corners = tsm.mk_im_corners(cra, cdec, tsm.cpa, 0.5)

    sp.plot_image_footprint_and_skycells(corners)

    pyplot.show()

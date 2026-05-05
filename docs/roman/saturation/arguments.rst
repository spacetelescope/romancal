Step Arguments
==============

``n_pix_grow_sat``
  Number of pixels by which to spatially expand each saturated pixel in all
  directions. Default is 0 (no expansion).

``backup``
  Number of times to extend saturation flags backwards to preceding
  multi-read resultants. When a resultant with multiple reads saturates,
  the algorithm may miss it if the per-read average stays below threshold.
  Setting ``backup=1`` conservatively flags the preceding multi-read
  resultant as well. Default is 0 (disabled).

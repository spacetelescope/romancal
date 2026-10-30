Add a ``sigma_stat`` parameter to choose the statistic used for clipping in
the relative and absolute fits, and make the default ``"mad"`` (previously
``"rmse"``), which is more robust to mismatched pairs in crowded fields.
Raise the default ``nclip`` and ``abs_nclip`` from 3 to 10; clipping stops
early once it converges.

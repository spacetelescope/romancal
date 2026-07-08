from romancal.associations import _AssociationPool
from romancal.associations.tests.helpers import t_path

POOL_FILE = t_path("data/jw93060_20150312T160130_pool.csv")


def test_pool(tmp_path):
    pool = _AssociationPool.read(POOL_FILE)
    assert len(pool) == 636

    file_path = tmp_path / __name__
    file_path.mkdir()

    tmp_pool = str(file_path / "tmp_pool.csv")
    pool.write(tmp_pool)

    roundtrip = _AssociationPool.read(tmp_pool)
    assert len(pool) == len(roundtrip)
    assert set(pool.colnames) == set(roundtrip.colnames)

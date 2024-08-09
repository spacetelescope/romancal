import tempfile
from pathlib import Path

import numpy as np


# not inheriting from MutableSequence here as insert is complicated
class ArrayLibrary:
    def __init__(self, tempdir=""):
        self._temp_dir = tempfile.TemporaryDirectory(dir=tempdir)
        self._temp_path = Path(self._temp_dir.name)
        self._filenames = []

    @property
    def closed(self):
        return not hasattr(self, "_temp_dir")

    def close(self):
        if self.closed:
            return
        self._temp_dir.cleanup()
        del self._temp_dir

    def __del__(self):
        self.close()

    def __len__(self):
        if self.closed:
            raise Exception("use after close")
        return len(self._filenames)

    def __getitem__(self, index):
        if self.closed:
            raise Exception("use after close")
        fn = self._filenames[index]
        return np.load(fn)

    def __setitem__(self, index, value):
        if self.closed:
            raise Exception("use after close")
        fn = self._filenames[index]
        if fn is None:
            fn = self._temp_path / f"{index}.npy"
        np.save(fn, value, False)
        self._filenames[index] = fn

    def append(self, value):
        if self.closed:
            raise Exception("use after close")
        index = len(self)
        self._filenames.append(None)
        self.__setitem__(index, value)

    def median(self):
        if self.closed:
            raise Exception("use after close")
        if not len(self):
            raise Exception("can't take median of empty list")
        return np.nanmedian(self, axis=0)
        # TODO make something like get_sections here
        # buffer = None
        # n_arrays = len(self)
        # allowed_memory = 100 << 20  # 100 MB

        # # figure out how big the buffer can be
        # allowed_memory_per_array = allowed_memory / n_arrays
        # # we'll allocated a buffer that is:
        # # [n_arrays, d0, some_subset_of_d1]
        # # but we don't yet know d0 so can't calculate the subset
        # # FIXME for now just load the first array
        # example_array = self[0]
        # dtype = example_array.dtype
        # shape = example_array.shape
        # if example_array.ndim != 2:
        #     raise Exception(f"Only works for 2 dimensions: {example_array.ndim}")
        # n_bytes_per_item = dtype.itemsize
        # n_dim_1 = allowed_memory_per_array // (n_bytes_per_item * shape[0])
        # if n_dim_1 < 1:
        #     raise Exception(f"Not enough memory")
        # buffer = np.empty((n_arrays, shape[0], n_dim_0

        # # first p
        # for i in range(len(
        # pass

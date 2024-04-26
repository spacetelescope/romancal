import glob
import io
import os

S3_BUCKET_NAME = "test-s3-data"


class MockS3Client:
    def __init__(self, s3_test_data_path):
        self.s3_test_data_path = s3_test_data_path

    def get_object(self, bucket_name, key):
        assert self.object_exists(bucket_name, key)

        with open(self._get_path(key), "rb") as f:
            return io.BytesIO(f.read())

    def object_exists(self, bucket_name, key):
        if bucket_name != S3_BUCKET_NAME:
            return False

        return os.path.isfile(self._get_path(key))

    def prefix_exists(self, bucket_name, key_prefix):
        return any(self.iterate_keys(bucket_name, key_prefix))

    def iterate_keys(self, bucket_name, key_prefix):
        if bucket_name != S3_BUCKET_NAME:
            return

        for k in self._list_keys():
            if k.startswith(key_prefix):
                yield k

    def _get_path(self, key):
        return os.path.join(self.s3_test_data_path, key)

    def _list_keys(self):
        paths = glob.glob(self.s3_test_data_path + "/**", recursive=True)
        paths = [p for p in paths if os.path.isfile(p)]
        return [p.replace(self.s3_test_data_path, "") for p in paths]


# Check strings based on words using length precision
def word_precision_check(str1, str2, length=5):
    """Check to strings word-by-word based for word length

    The strings are checked word for word, but only for the first
    `length` characters

    Parameters
    ----------
    str1, str2: str
        The strings to compare

    length: int
        The number of characters in each word to check.

    Returns
    -------
    match: boolean
        True if the strings match

    Raises
    ------
    AssertionError
        Raised if the number of word differ or the beginning components of each word differ.
    """
    words1 = str1.split()
    words2 = str2.split()
    if len(words1) != len(words2):
        raise AssertionError(
            f"str1 has different number of words {len(words1)} than str2 {len(words2)}"
        )
    for w1, w2 in zip(words1, words2):
        if w1[:length] != w2[:length]:
            raise AssertionError(f"str1 word {w1[:length]} != str2 {w2[:length]}")
    return True

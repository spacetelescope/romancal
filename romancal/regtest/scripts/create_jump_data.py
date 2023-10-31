import pickle
from pathlib import Path
from typing import NamedTuple

import astropy.units as u
import numpy as np
from roman_datamodels.datamodels import GainRefModel, RampModel, ReadnoiseRefModel
from roman_datamodels.maker_utils import mk_datamodel

REFERENCE_DATA_PATH = "/user/wjamieson/jump_detection/reference_data/"


class RegressionData(NamedTuple):
    input_model: RampModel
    readnoise_model: ReadnoiseRefModel
    gain_model: GainRefModel


def ma_table_to_read_pattern(ma_table, read_time):
    return [(entry / read_time).astype(int).tolist() for entry in ma_table]


def reference_data_to_input_models(input_data):
    counts = input_data["counts"].astype(np.float32)
    n_resultants, n_counts = counts.shape
    n_rows = np.sqrt(n_resultants).astype(np.int32)
    if n_rows**2 != n_resultants:
        raise ValueError("Number of resultants must be a square number")
    counts = counts.transpose().reshape((n_counts, n_rows, n_rows))

    padded = np.zeros((n_counts, n_rows + 8, n_rows + 8), dtype=np.float32)
    padded[:, 4:-4, 4:-4] = counts

    ma_table = input_data["params"]["ma_table"]
    read_time = ma_table[0][0]
    read_pattern = ma_table_to_read_pattern(ma_table, read_time)

    # output_slopes = input_data["stats"]["est_slope"]
    # print("Making output_model no jump")
    # output_model = mk_datamodel(ImageModel, shape=padded.shape[1:])

    print("Making input_model")
    input_model = mk_datamodel(RampModel, shape=padded.shape)
    input_model.data = padded * u.electron
    input_model.meta.exposure.read_pattern = read_pattern
    input_model.meta.exposure.frame_time = read_time

    print("Making readnoise_model")
    read_noise = input_data["params"]["read_noise"]
    readnoise_model = mk_datamodel(ReadnoiseRefModel, shape=padded.shape[1:])
    readnoise_model.data += read_noise * readnoise_model.data.unit

    print("Making gain_model")
    gain_model = mk_datamodel(GainRefModel, shape=padded.shape[1:])
    gain_model.data += 1 * gain_model.data.unit

    return RegressionData(input_model, readnoise_model, gain_model)


def read_data(file_path):
    print(f"Using {file_path=}")
    data_files = Path(file_path).glob("*.pkl")

    for data_file in data_files:
        print(f"reading: {data_file}")
        with data_file.open("rb") as f:
            input_data = pickle.load(f)

        yield reference_data_to_input_models(input_data), data_file.stem


def main(file_path=None):
    if file_path is None:
        file_path = REFERENCE_DATA_PATH

    for data, file_name in read_data(file_path):
        print(f"writing: {file_name} data files")
        data.input_model.to_asdf(f"{file_name}_input.asdf")
        data.readnoise_model.to_asdf(f"{file_name}_readnoise.asdf")
        data.gain_model.to_asdf(f"{file_name}_gain.asdf")

    print("Done")


if __name__ == "__main__":
    main()

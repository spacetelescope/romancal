from romancal.datamodels.container import ModelContainer
from roman_datamodels import datamodels as rdm
from roman_datamodels import maker_utils as utils


def test_model_container_input_as_list_of_datamodels(tmp_path):
    filepath1 = (
        tmp_path / "test_model_container_input_as_list_of_filepaths_01.asdf"
    )
    filepath2 = (
        tmp_path / "test_model_container_input_as_list_of_filepaths_02.asdf"
    )
    # create L2 file using filepath
    utils.mk_level2_image(filepath=filepath1)
    utils.mk_level2_image(filepath=filepath2)
    # parse ASDF file as RDM
    datamodel1 = rdm.open(str(filepath1))
    datamodel2 = rdm.open(str(filepath2))
    # append datamodel to datamodel list
    datamodel_list = [datamodel1, datamodel2]

    # provide datamodel list as input to ModelContainer
    model_container = ModelContainer(datamodel_list)

    assert all(isinstance(x, rdm.DataModel) for x in model_container)


def test_modelcontainer_init():
    img = utils.mk_level1_science_raw()
    m = rdm.ImageModel(img)
    mc1 = ModelContainer([m])

    # initialize with an instance of ModelContainer
    mc2 = ModelContainer(mc1)

    assert len(mc1) == len(mc2)

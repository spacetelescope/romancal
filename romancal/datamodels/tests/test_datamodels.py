import pytest
from roman_datamodels import datamodels as rdm
from roman_datamodels import maker_utils as utils

from romancal.datamodels.container import ModelContainer


@pytest.fixture(name="set_up_list_of_l2_files")
def set_up_list_of_l2_files(tmp_path, request):
    # generate a list of n filepaths and files to be read later on by ModelContainer
    marker = request.node.get_closest_marker("set_up_list_of_l2_files_data")
    number_of_files_to_create = marker.args[0]
    type_of_returned_object = marker.args[1]

    result_list = []
    for i in range(number_of_files_to_create):
        filepath = (
            tmp_path / f"test_model_container_input_as_list_of_filepaths_{i:02}.asdf"
        )
        # create L2 file using filepath
        utils.mk_level2_image(filepath=filepath)
        if type_of_returned_object == "asdf":
            # append filepath to filepath list
            result_list.append(str(filepath))
        elif type_of_returned_object == "datamodel":
            # parse ASDF file as RDM
            datamodel = rdm.open(str(filepath))
            # append datamodel to datamodel list
            result_list.append(datamodel)

    return result_list


@pytest.mark.set_up_list_of_l2_files_data(2, "asdf")
def test_model_container_input_as_list_of_filepaths(set_up_list_of_l2_files):
    filepath_list = set_up_list_of_l2_files
    # provide filepath list as input to ModelContainer
    model_container = ModelContainer(filepath_list)

    assert len(model_container) == 2
    # check if all model_container elements are instances of DataModel
    assert all(isinstance(x, rdm.DataModel) for x in model_container)


@pytest.mark.set_up_list_of_l2_files_data(2, "datamodel")
def test_model_container_input_as_list_of_datamodels(set_up_list_of_l2_files):
    filepath_list = set_up_list_of_l2_files
    # provide filepath list as input to ModelContainer
    model_container = ModelContainer(filepath_list)

    assert len(model_container) == 2
    # check if all model_container elements are instances of DataModel
    assert all(isinstance(x, rdm.DataModel) for x in model_container)


@pytest.mark.set_up_list_of_l2_files_data(2, "asdf")
def test_imagemodel_set_item(set_up_list_of_l2_files):
    filepath_list = set_up_list_of_l2_files
    # provide filepath list as input to ModelContainer
    model_container = ModelContainer(filepath_list)
    # grab first datamodel for testing
    image_model = model_container[0]

    image_model["test_attr"] = "test_attr_value"
    assert hasattr(image_model, "test_attr")
    assert image_model.test_attr == "test_attr_value"
    image_model["test_attr"] = "test_attr_new_value"
    assert image_model.test_attr == "test_attr_new_value"
    # test ValueError
    with pytest.raises(Exception) as e:
        image_model["_test_attr"] = "test_attr_some_other_value"
    assert e.type == ValueError


def test_modelcontainer_init():
    img = utils.mk_level1_science_raw()
    m = rdm.ImageModel(img)
    mc1 = ModelContainer([m])

    # initialize with an instance of ModelContainer
    mc2 = ModelContainer(mc1)

    assert len(mc1) == len(mc2)

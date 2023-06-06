from io import StringIO
from pathlib import Path
import pytest
from roman_datamodels import datamodels as rdm
from roman_datamodels import maker_utils as utils

from romancal.datamodels.container import ModelContainer


def create_asn_file(tmp_path):
    asn_content = """
        {
            "asn_type": "None",
            "asn_rule": "DMS_ELPP_Base",
            "version_id": null,
            "code_version": "0.9.1.dev28+ge987cc9.d20230106",
            "degraded_status": "No known degraded exposures in association.",
            "program": "noprogram",
            "constraints": "No constraints",
            "asn_id": "a3001",
            "target": "none",
            "asn_pool": "test_pool_name",
            "products": [
                {
                    "name": "files.asdf",
                    "members": [
                        {
                            "expname": "img_1.asdf",
                            "exptype": "science"
                        },
                        {
                            "expname": "img_2.asdf",
                            "exptype": "science"
                        }
                    ]
                }
            ]
        }
"""
    asn_file_path = str(tmp_path / "sample_asn.json")
    asn_file = StringIO()
    asn_file.write(asn_content)
    with open(asn_file_path, mode="w") as f:
        print(asn_file.getvalue(), file=f)

    return asn_file_path


@pytest.fixture()
def setup_list_of_l2_files():
    def _setup_list_of_l2_files(n, obj_type, tmp_path):
        """
        Generate a list of `n` ASDF files (and their corresponding path) or datamodels.

        Parameters
        ----------
        n : int
            The number of ASDF files or datamodels to be generated.
        obj_type : str
            The type of object to be generated. Allowed values: "asdf" or "datamodel".
        tmp_path : _type_
            The dir path where the generated ASDF files will be temporarily saved to.

        Returns
        -------
        list
            A list containing either the full path to an ASDF file or datamodels.
        """
        number_of_files_to_create = n
        type_of_returned_object = obj_type

        result_list = []
        for i in range(number_of_files_to_create):
            filepath = (
                tmp_path
                / f"test_model_container_input_as_list_of_filepaths_{i:02}.asdf"
            )
            # create an ASDF file with an L2 model
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

    return _setup_list_of_l2_files


def test_model_container_init_with_modelcontainer_instance():
    img = utils.mk_level2_image()
    m = rdm.ImageModel(img)
    mc1 = ModelContainer([m])

    # initialize with an instance of ModelContainer
    mc2 = ModelContainer(mc1)

    assert isinstance(mc2, ModelContainer)
    assert all(isinstance(x, rdm.DataModel) for x in mc1)
    assert len(mc1) == len(mc2)


@pytest.mark.parametrize("n, obj_type", [(3, "asdf"), (2, "datamodel")])
def test_model_container_init_path_to_asdf_or_datamodels(
    n, obj_type, tmp_path, request
):
    filepath_list = request.getfixturevalue("setup_list_of_l2_files")(
        n, obj_type, tmp_path
    )
    mc = ModelContainer(filepath_list)

    assert all(x in filepath_list for x in mc._models)


def test_model_container_init_with_path_to_asn_file(tmp_path):
    # create ASDF files with L2 datamodel
    utils.mk_level2_image(filepath=tmp_path / "img_1.asdf")
    utils.mk_level2_image(filepath=tmp_path / "img_2.asdf")
    # create ASN file that points to the ASDF files
    asn_filepath = create_asn_file(tmp_path)
    mc = ModelContainer(asn_filepath)

    assert all(hasattr(x.meta, "asn") for x in mc)


@pytest.mark.parametrize(
    "input_object",
    [
        "invalid_object",
        "",
        [1, 2, 3],
        ModelContainer(),
        Path(),
    ],
)
def test_imagemodel_init_error(input_object):
    with pytest.raises(Exception) as e:
        mc = ModelContainer(input_object)

    assert e.type == TypeError


@pytest.mark.parametrize("n, obj_type", [(4, "asdf"), (5, "datamodel")])
def test_imagemodel_slice_n_dice(n, obj_type, tmp_path, request):
    filepath_list = request.getfixturevalue("setup_list_of_l2_files")(
        n, obj_type, tmp_path
    )
    # provide filepath list as input to ModelContainer
    mc = ModelContainer(filepath_list)

    x1 = mc[:]
    x2 = mc[0]
    x3 = mc[2:]
    x4 = mc[-2]
    x5 = mc[-2:]

    assert type(x1) == list
    assert len(x1) == len(filepath_list)

    assert type(x2) == rdm.ImageModel

    assert type(x3) == list
    assert len(x3) == n - 2

    assert type(x4) == rdm.ImageModel

    assert type(x5) == list
    assert len(x5) == 2


def test_imagemodel_set_item(setup_list_of_l2_files, tmp_path):
    filepath_list = setup_list_of_l2_files(4, "datamodel", tmp_path)
    # provide filepath list as input to ModelContainer
    mc1 = ModelContainer(filepath_list[:2])
    mc2 = ModelContainer(filepath_list[-2:])

    mc1[0] = mc2[-2]
    mc1[1] = mc2[-1]

    assert all(id(l) == id(r) for l, r in zip(mc1[:], mc2[:]))


@pytest.mark.parametrize(
    "n, obj_type, input_object",
    [
        (2, "datamodel", "invalid_object"),
        (2, "datamodel", ""),
        (2, "datamodel", [1, 2, 3]),
        (2, "datamodel", ModelContainer()),
    ],
)
def test_imagemodel_set_item_error(
    n, obj_type, input_object, tmp_path, request
):
    filepath_list = request.getfixturevalue("setup_list_of_l2_files")(
        n, obj_type, tmp_path
    )
    # provide filepath list as input to ModelContainer
    mc = ModelContainer(filepath_list)

    with pytest.raises(Exception) as e:
        mc[0] = input_object

    assert e.type == ValueError


@pytest.mark.parametrize("n, obj_type", [(3, "datamodel")])
def test_model_container_insert(n, obj_type, tmp_path, request):
    filepath_list = request.getfixturevalue("setup_list_of_l2_files")(
        n, obj_type, tmp_path
    )
    mc = ModelContainer(filepath_list[:2])

    mc.insert(1, filepath_list[-1])

    assert len(mc) == n
    assert id(mc[1]) == id(filepath_list[-1])


@pytest.mark.parametrize("n, obj_type", [(3, "datamodel")])
def test_model_container_insert_error(n, obj_type, tmp_path, request):
    filepath_list = request.getfixturevalue("setup_list_of_l2_files")(
        n, obj_type, tmp_path
    )
    mc = ModelContainer(filepath_list)

    with pytest.raises(Exception) as e:
        # try to insert a ModelContainer
        mc.insert(1, ModelContainer())

    assert e.type == ValueError


@pytest.mark.parametrize("n, obj_type", [(3, "datamodel")])
def test_model_container_append(n, obj_type, tmp_path, request):
    filepath_list = request.getfixturevalue("setup_list_of_l2_files")(
        n, obj_type, tmp_path
    )
    mc = ModelContainer(filepath_list[:2])

    mc.append(filepath_list[-1])

    assert len(mc) == n
    assert id(mc[-1]) == id(filepath_list[-1])


@pytest.mark.parametrize("n, obj_type", [(3, "datamodel")])
def test_model_container_append_error(n, obj_type, tmp_path, request):
    filepath_list = request.getfixturevalue("setup_list_of_l2_files")(
        n, obj_type, tmp_path
    )
    mc = ModelContainer(filepath_list)

    with pytest.raises(Exception) as e:
        # try to append a ModelContainer
        mc.append(ModelContainer())

    assert e.type == ValueError


@pytest.mark.parametrize("n, obj_type", [(4, "datamodel")])
def test_model_container_extend(n, obj_type, tmp_path, request):
    filepath_list = request.getfixturevalue("setup_list_of_l2_files")(
        n, obj_type, tmp_path
    )
    mc1 = ModelContainer(filepath_list[:2])
    mc2 = ModelContainer(filepath_list[-2:])

    mc1.extend(mc2)

    assert len(mc1) == n
    assert all(id(l) == id(r) for l, r in zip(mc1[-2:], filepath_list[-2:]))


@pytest.mark.parametrize(
    "n, obj_type, input_object",
    [
        (3, "datamodel", ["trying_to_sneak_in", 1, 2, 3]),
        (3, "datamodel", ""),
        (3, "datamodel", "trying_to_sneak_in"),
    ],
)
def test_model_container_extend_error(
    n, obj_type, input_object, tmp_path, request
):
    filepath_list = request.getfixturevalue("setup_list_of_l2_files")(
        n, obj_type, tmp_path
    )
    mc = ModelContainer(filepath_list)

    with pytest.raises(Exception) as e:
        # try to insert a string
        mc.extend(input_object)

    assert e.type == ValueError


@pytest.mark.parametrize("n, obj_type", [(3, "datamodel")])
def test_model_container_pop_last_item(n, obj_type, tmp_path, request):
    filepath_list = request.getfixturevalue("setup_list_of_l2_files")(
        n, obj_type, tmp_path
    )
    mc = ModelContainer(filepath_list)

    mc.pop()

    assert len(mc) == n - 1
    assert filepath_list[-1] not in mc[:]
    assert all(id(l) == id(r) for l, r in zip(mc[:], filepath_list[:]))


@pytest.mark.parametrize("n, obj_type", [(3, "datamodel")])
def test_model_container_pop_with_index(n, obj_type, tmp_path, request):
    filepath_list = request.getfixturevalue("setup_list_of_l2_files")(
        n, obj_type, tmp_path
    )
    mc = ModelContainer(filepath_list)

    index_to_be_removed = 1
    mc.pop(index_to_be_removed)

    assert len(mc) == n - 1
    assert filepath_list[index_to_be_removed] not in mc[:]


@pytest.mark.parametrize("n, obj_type", [(2, "asdf"), (3, "datamodel")])
def test_model_container_copy(n, obj_type, tmp_path, request):
    filepath_list = request.getfixturevalue("setup_list_of_l2_files")(
        n, obj_type, tmp_path
    )

    mc = ModelContainer(filepath_list)
    mc_copy_dict = mc.copy()

    mc_dict = mc.__dict__
    mc_copy_dict = mc_copy_dict.__dict__

    assert id(mc_dict) != id(mc_copy_dict)
    assert all(x in mc_dict for x in mc_copy_dict)
    assert all(
        type(o) == type(c)
        for o, c in zip(mc_copy_dict.values(), mc_dict.values())
    )

import json
from contextlib import nullcontext

import pytest
import roman_datamodels.datamodels as dm
from roman_datamodels.maker_utils import mk_level2_image

from romancal.associations import load_asn
from romancal.associations.asn_from_list import asn_from_list
from romancal.datamodels.library import BorrowError, ClosedLibraryError, ModelLibrary

# for the example association, set 2 different observation numbers
# so the association will have 2 groups (since all other group_id
# determining meta is the same, see `example_asn_path`)
_OBSERVATION_NUMBERS = [1, 1, 2]
_N_MODELS = len(_OBSERVATION_NUMBERS)
_N_GROUPS = len(set(_OBSERVATION_NUMBERS))
_PRODUCT_NAME = "foo_out"


@pytest.fixture
def example_asn_path(tmp_path):
    """
    Fixture that creates a simple association, saves it (and the models)
    to disk, and returns the path of the saved association
    """
    fns = []
    for i in range(_N_MODELS):
        m = dm.ImageModel(mk_level2_image(shape=(2, 2)))
        m.meta.observation.program = "0001"
        m.meta.observation.observation = _OBSERVATION_NUMBERS[i]
        m.meta.observation.visit = 1
        m.meta.observation.visit_file_group = 1
        m.meta.observation.visit_file_sequence = 1
        m.meta.observation.visit_file_activity = "01"
        m.meta.observation.exposure = 1
        base_fn = f"{i}.asdf"
        m.meta.filename = base_fn
        m.save(str(tmp_path / base_fn))
        fns.append(base_fn)
    asn = asn_from_list(fns, product_name=_PRODUCT_NAME)
    base_fn, contents = asn.dump(format="json")
    asn_filename = tmp_path / base_fn
    with open(asn_filename, "w") as f:
        f.write(contents)
    return asn_filename


@pytest.fixture
def example_library(example_asn_path):
    """
    Fixture that builds off of `example_asn_path` and returns a
    library created from the association with default options
    """
    return ModelLibrary(example_asn_path)


def _set_custom_member_attr(example_asn_path, member_index, attr, value):
    """
    Helper function to modify the association at `example_asn_path`
    by adding an attribute `attr` to the member list (at index
    `member_index`) with value `value`. This is used to modify
    the `group_id` or `exptype` of a certain member for some tests.
    """
    with open(example_asn_path) as f:
        asn_data = load_asn(f)
    asn_data["products"][0]["members"][member_index][attr] = value
    with open(example_asn_path, "w") as f:
        json.dump(asn_data, f)


def test_load_asn(example_library):
    """
    Test that __len__ returns the number of models/members loaded
    from the association (and does not require opening the library)
    """
    assert len(example_library) == _N_MODELS


@pytest.mark.parametrize("asn_n_members", range(_N_MODELS))
def test_asn_n_members(example_asn_path, asn_n_members):
    """
    Test that creating a library with a `asn_n_members` filter
    includes only the first N members
    """
    library = ModelLibrary(example_asn_path, asn_n_members=asn_n_members)
    assert len(library) == asn_n_members


def test_asn_exptypes(example_asn_path):
    """
    Test that creating a library with a `asn_exptypes` filter
    includes only the members with a matching `exptype`
    """
    _set_custom_member_attr(example_asn_path, 0, "exptype", "background")
    library = ModelLibrary(example_asn_path, asn_exptypes="science")
    assert len(library) == _N_MODELS - 1
    library = ModelLibrary(example_asn_path, asn_exptypes="background")
    assert len(library) == 1


def test_group_names(example_library):
    """
    Test that `group_names` returns appropriate names
    based on the inferred group ids and that these names match
    the `model.meta.group_id` values
    """
    assert len(example_library.group_names) == _N_GROUPS
    group_names = set()
    with example_library:
        for index, model in enumerate(example_library):
            group_names.add(model.meta.group_id)
            example_library.discard(index, model)
    assert group_names == set(example_library.group_names)


def test_group_indices(example_library):
    """
    Test that `group_indices` returns appropriate model indices
    based on the inferred group ids
    """
    group_indices = example_library.group_indices
    assert len(group_indices) == _N_GROUPS
    with example_library:
        for group_name in group_indices:
            indices = group_indices[group_name]
            for index in indices:
                model = example_library[index]
                assert model.meta.group_id == group_name
                example_library.discard(index, model)


@pytest.mark.parametrize("attr", ["group_names", "group_indices"])
def test_group_with_no_datamodels_open(example_asn_path, attr, monkeypatch):
    """
    Test that the "grouping" methods do not call datamodels.open
    """

    # patch datamodels.open to always raise an exception
    # this will serve as a smoke test to see if any of the attribute
    # accesses (or instance creation) attempts to open models
    def no_open(*args, **kwargs):
        raise Exception()

    monkeypatch.setattr(dm, "open", no_open)

    # use example_asn_path here to make the instance after we've patched
    # datamodels.open
    library = ModelLibrary(example_asn_path)
    getattr(library, attr)


# @pytest.mark.parametrize(
#     "asn_group_id, meta_group_id, expected_group_id", [
#         ('42', None, '42'),
#         (None, '42', '42'),
#         ('42', '26', '42'),
#     ])
# def test_group_id_override(example_asn_path, asn_group_id, meta_group_id, expected_group_id):
#     """
#     Test that overriding a models group_id via:
#         - the association member entry
#         - the model.meta.group_id
#     overwrites the automatically calculated group_id (with the asn taking precedence)
#     """
#     if asn_group_id:
#         _set_custom_member_attr(example_asn_path, 0, 'group_id', asn_group_id)
#     if meta_group_id:
#         model_filename = example_asn_path.parent / '0.fits'
#         with dm.open(model_filename) as model:
#             model.meta.group_id = meta_group_id
#             model.save(model_filename)
#     library = ModelLibrary(example_asn_path)
#     group_names = library.group_names
#     assert len(group_names) == 3
#     assert expected_group_id in group_names
#     with library:
#         model = library[0]
#         assert model.meta.group_id == expected_group_id
#         library.discard(0, model)


@pytest.mark.parametrize("return_method", ("__setitem__", "discard"))
def test_model_iteration(example_library, return_method):
    """
    Test that iteration through models and returning (or discarding) models
    returns the appropriate models
    """
    with example_library:
        for i, model in enumerate(example_library):
            assert int(model.meta.filename.split(".")[0]) == i
            getattr(example_library, return_method)(i, model)


@pytest.mark.parametrize("return_method", ("__setitem__", "discard"))
def test_model_indexing(example_library, return_method):
    """
    Test that borrowing models (using __getitem__)  and returning (or discarding)
    models returns the appropriate models
    """
    with example_library:
        for i in range(_N_MODELS):
            model = example_library[i]
            assert int(model.meta.filename.split(".")[0]) == i
            getattr(example_library, return_method)(i, model)


def test_closed_library_model_getitem(example_library):
    """
    Test that indexing a library when it is not open triggers an error
    """
    with pytest.raises(ClosedLibraryError, match="ModelLibrary is not open"):
        example_library[0]


def test_closed_library_model_iter(example_library):
    """
    Test that attempting to iterate a library that is not open triggers an error
    """
    with pytest.raises(ClosedLibraryError, match="ModelLibrary is not open"):
        for model in example_library:
            pass


def test_double_borrow_by_index(example_library):
    """
    Test that double-borrowing a model (using __getitem__) results in an error
    """
    with pytest.raises(BorrowError, match="1 un-returned models"):
        with example_library:
            model0 = example_library[0]  # noqa: F841
            with pytest.raises(BorrowError, match="Attempt to double-borrow model"):
                model1 = example_library[0]  # noqa: F841


def test_double_borrow_during_iter(example_library):
    """
    Test that double-borrowing a model (once via iter and once via __getitem__)
    results in an error
    """
    with pytest.raises(BorrowError, match="1 un-returned models"):
        with example_library:
            for index, model in enumerate(example_library):
                with pytest.raises(BorrowError, match="Attempt to double-borrow model"):
                    model1 = example_library[index]  # noqa: F841
                break


def test_non_borrowed_setitem(example_library):
    """
    Test that attempting to return a non-borrowed item results in an error
    """
    with example_library:
        with pytest.raises(BorrowError, match="Attempt to return non-borrowed model"):
            example_library[0] = None


def test_non_borrowed_discard(example_library):
    """
    Test that attempting to discard a non-borrowed item results in an error
    """
    with example_library:
        with pytest.raises(BorrowError, match="Attempt to discard non-borrowed model"):
            example_library.discard(0, None)


@pytest.mark.parametrize("n_borrowed", (1, 2))
def test_no_return_getitem(example_library, n_borrowed):
    """
    Test that borrowing and not returning models results in an
    error noting the number of un-returned models.
    """
    with pytest.raises(
        BorrowError, match=f"ModelLibrary has {n_borrowed} un-returned models"
    ):
        with example_library:
            for i in range(n_borrowed):
                example_library[i]


def test_exception_while_open(example_library):
    """
    Test that the __exit__ implementation for the library
    passes exceptions that occur in the context
    """
    with pytest.raises(Exception, match="test"):
        with example_library:
            raise Exception("test")


def test_exception_with_borrow(example_library):
    """
    Test that an exception while the library is open and has a borrowed
    model results in a chained exception containing both:
        - the original exception (as the __context__)
        - an exception about the un-returned model
    """
    with pytest.raises(BorrowError, match="1 un-returned models") as exc_info:
        with example_library:
            model = example_library[0]  # noqa: F841
            raise Exception("test")
    # check that Exception above is the __context__ (in the chain)
    assert exc_info.value.__context__.__class__ is Exception
    assert exc_info.value.__context__.args == ("test",)


def test_asn_data(example_library):
    """
    Test that `asn` returns the association information
    """
    assert example_library.asn["products"][0]["name"] == _PRODUCT_NAME


def test_asn_readonly(example_library):
    """
    Test that modifying the product (dict) in the `asn` result triggers an exception
    """
    with pytest.raises(TypeError, match="object does not support item assignment"):
        example_library.asn["products"][0]["name"] = f"{_PRODUCT_NAME}_new"


def test_asn_members_readonly(example_library):
    """
    Test that modifying members (list) in the `asn` result triggers an exception
    """
    with pytest.raises(TypeError, match="object does not support item assignment"):
        example_library.asn["products"][0]["members"][0]["group_id"] = "42"


def test_asn_members_tuple(example_library):
    """
    Test that even nested items in `asn` (like `members`) are immutable
    """
    assert isinstance(example_library.asn["products"][0]["members"], tuple)


# def test_members(example_library):
#     assert example_library.asn['products'][0]['members'] == example_library.members
#
#
# def test_members_tuple(example_library):
#     assert isinstance(example_library.members, tuple)


@pytest.mark.parametrize("n, err", [(1, False), (2, True)])
def test_stpipe_models_access(example_asn_path, n, err):
    """
    stpipe currently reaches into _models (but only when asn_n_members
    is 1) so we support this `_models` attribute (with a loaded model)
    only under that condition until stpipe can be updated to not reach
    into `_models`.
    """
    library = ModelLibrary(example_asn_path, asn_n_members=n)
    if err:
        ctx = pytest.raises(AttributeError, match="object has no attribute '_models'")
    else:
        ctx = nullcontext()
    with ctx:
        assert library._models[0].get_crds_parameters()


@pytest.mark.parametrize("discard", [True, False])
def test_on_disk_model_modification(example_asn_path, discard):
    """
    Test that modifying a model in a library that is on_disk
    does not persist if the model is discarded (instead of
    returned via __setitem__)
    """
    library = ModelLibrary(example_asn_path, on_disk=True)
    with library:
        model = library[0]
        model.meta["foo"] = "bar"
        if discard:
            library.discard(0, model)
        else:
            library[0] = model
        model = library[0]
        if discard:
            # since the model was 'discarded' and the library is 'on_disk'
            # the modification should not persist
            assert getattr(model.meta, "foo", None) is None
        else:
            # if instead, we used __setitem__ the modification should be saved
            assert getattr(model.meta, "foo") == "bar"
        library.discard(0, model)


@pytest.mark.parametrize("on_disk", [True, False])
def test_on_disk_no_overwrite(example_asn_path, on_disk):
    """
    Test that modifying a model in a library does not overwrite
    the input file (even if on_disk==True)
    """
    library = ModelLibrary(example_asn_path, on_disk=on_disk)
    with library:
        model = library[0]
        model.meta["foo"] = "bar"
        library[0] = model

    library2 = ModelLibrary(example_asn_path, on_disk=on_disk)
    with library2:
        model = library2[0]
        assert getattr(model.meta, "foo", None) is None
        library2[0] = model


# TODO container conversion
# TODO index
# TODO memmap?

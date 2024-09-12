import json
from pathlib import Path

TEST_REQUIREMENTS_FILENAME = Path(__file__).parent / "dms_requirement_tests.json"


def test_requirements(all_test_names):
    test_requirements_filename = TEST_REQUIREMENTS_FILENAME.expanduser().absolute()

    with open(test_requirements_filename) as test_requirements_file:
        requirements = json.load(test_requirements_file)

    required_test_names = {
        test
        for requirement_tests in requirements.values()
        for test in requirement_tests
    }

    missing_test_names = required_test_names - all_test_names
    assert not missing_test_names, f"Missing tests {missing_test_names} required by DMS"

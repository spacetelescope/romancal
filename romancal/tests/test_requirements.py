import json
import re
from pathlib import Path

TEST_REQUIREMENTS_FILENAME = (
    Path(__file__).parent.parent.parent / "test_requirements.json"
)
TEST_DIRECTORY = Path(__file__).parent.parent


def test_requirements():
    test_requirements_filename = TEST_REQUIREMENTS_FILENAME.expanduser().absolute()
    test_directory = TEST_DIRECTORY.expanduser().absolute()

    with open(test_requirements_filename) as test_requirements_file:
        requirements = json.load(test_requirements_file)

    required_tests = sorted(
        {
            test
            for requirement_tests in requirements.values()
            for test in requirement_tests
        }
    )

    existing_tests = []
    test_regex = re.compile(r"def (test_[^\(]+)\(.*\):")
    for test_filename in test_directory.glob("**/test_*.py"):
        with open(test_filename) as test_file:
            test_file_contents = test_file.read()

        for match in re.finditer(test_regex, test_file_contents):
            test = f"{test_directory.stem}.{str(test_filename.relative_to(test_directory).parent).replace('/', '.')}.{test_filename.stem}.{match.group(1)}"
            if test in required_tests:
                existing_tests.append(test)

    existing_tests = sorted(existing_tests)

    assert existing_tests == required_tests

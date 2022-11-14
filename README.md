> # Roman Calibration Pipeline


[![Documentation Status](https://readthedocs.org/projects/roman-pipeline/badge/?version=latest)](https://roman-pipeline.readthedocs.io/en/latest/?badge=latest)
[![CI](https://github.com/spacetelescope/romancal/actions/workflows/roman_ci.yml/badge.svg)](https://github.com/spacetelescope/romancal/actions/workflows/roman_ci.yml)
[![Powered by STScI Badge](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)
[![Powered by Astropy Badge](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

![STScI Logo](docs/_static/stsci_logo.png)

**Roman requires Python 3.8 or above and a C compiler for dependencies.**

**Linux and MacOS platforms are tested and supported. Windows is not currently supported.**

## Installation

The easiest way to install the latest `romancal` release into a fresh virtualenv or conda environment is

    pip install romancal

### Detailed Installation

The `romancal` package can be installed into a virtualenv or conda environment via `pip`. We recommend that for each
installation you start by creating a fresh environment that only has Python installed and then install the `romancal`
package and its dependencies into that bare environment. If using conda environments, first make sure you have a recent
version of Anaconda or Miniconda installed. If desired, you can create multiple environments to allow for switching
between different versions of the `romancal` package (e.g. a released version versus the current development version).

In all cases, the installation is generally a 3-step process:

* Create a conda environment
* Activate that environment
* Install the desired version of the `romancal` package into that environment

Details are given below on how to do this for different types of installations, including tagged releases, DMS builds
used in operations, and development versions. Remember that all conda operations must be done from within a bash shell.

### Installing latest releases

You can install the latest released version via `pip`. From a bash shell:

    conda create -n <env_name> python
    conda activate <env_name>
    pip install romancal

You can also install a specific version (from `romancal 0.1.0` onward):

    conda create -n <env_name> python
    conda activate <env_name>
    pip install romancal==0.5.0

### Installing the development version from Github

You can install the latest development version (not as well tested) from the Github main branch:

    conda create -n <env_name> python
    conda activate <env_name>
    pip install git+https://github.com/spacetelescope/romancal

### Installing for Developers

If you want to be able to work on and test the source code with the `romancal` package, the high-level procedure to do
this is to first create a conda environment using the same procedures outlined above, but then install your personal
copy of the code overtop of the original code in that environment. Again, this should be done in a separate conda
environment from any existing environments that you may have already installed with released versions of the `romancal`
package.

As usual, the first two steps are to create and activate an environment:

    conda create -n <env_name> python
    conda activate <env_name>

To install your own copy of the code into that environment, you first need to fork and clone the `romancal` repo:

    cd <where you want to put the repo>
    git clone https://github.com/spacetelescope/romancal
    cd romancal

*Note: `python setup.py install` and `python setup.py develop` commands do not work.*

Install from your local checked-out copy as an "editable" install:

    pip install -e .

If you want to run the unit or regression tests and/or build the docs, you can make sure those dependencies are
installed too:

    pip install -e .[test]
    pip install -e .[docs]
    pip install -e .[test,docs]

Need other useful packages in your development environment?

    pip install ipython pytest-xdist

## Calibration References Data System (CRDS) Setup

CRDS is the system that manages the reference files needed to run the pipeline. Inside the STScI network, the pipeline
works with default CRDS setup with no modifications. To run the pipeline outside the STScI network, CRDS must be
configured by setting two environment variables:

    export CRDS_PATH=$HOME/crds_cache
    export CRDS_SERVER_URL=https://roman-crds-test.stsci.edu

## Documentation

Documentation (built daily from the Github `main` branch) is available at:

https://roman-pipeline.readthedocs.io/en/latest/

To build the docs yourself, clone this repository and build the documentation with:

    pip install -e .[docs]
    cd docs
    make html

## Contributions and Feedback

We welcome contributions and feedback on the project. Please follow the
[contributing guidelines](CONTRIBUTING.md) to submit an issue or a pull request.

We strive to provide a welcoming community to all of our users by abiding with the [Code of Conduct](CODE_OF_CONDUCT.md)
.

If you have questions or concerns regarding the software, please open
an [issue](https://github.com/spacetelescope/romancal/issues).

## Software vs DMS build version map

| roman tag | DMS build | CRDS_CONTEXT | Date      | Notes                                 |
|-----------|-----------|--------------|-----------|---------------------------------------|
| 0.1.0     | 0.0       | 003          | Nov  2020 | Release for Build 0.0                 |
| 0.2.0     | 0.1       | 004          | Mar  2021 | Release for Build 0.1                 |
| 0.3.0     | 0.2       | 007          | May  2021 | Release for Build 0.2                 |
| 0.3.1     | 0.2       | 007          | Jun  2021 | Release for Build 0.2 CRDS tests      |
| 0.4.2     | 0.3       | 011          | Sep  2021 | Release for Build 0.3                 |
| 0.5.0     | 0.4       | 023          | Dec  2021 | Release for Build 0.4                 |
| 0.6.0     | 0.5       | 030          | Mar  2022 | Release for Build 0.5                 |
| 0.7.0     | 22Q3_B6   | 032          | May  2022 | Release for Build 22Q3_B6 (Build 0.6) |
| 0.7.1     | 22Q3_B6   | 032          | May  2022 | Release for Build 22Q3_B6 (Build 0.6) |
| 0.8.0     | 22Q4_B7   | 038          | Aug  2022 | Release for Build 22Q4_B7 (Build 0.7) |
| 0.8.1     | 22Q4_B7   | 038          | Aug  2022 | Release for Build 22Q4_B7 (Build 0.7) |
| 0.9.0     | 23Q1_B8   | 039          | Nov  2022 | Release for Build 23Q1_B8 (Build 8)   |

Note: CRDS_CONTEXT values flagged with an asterisk in the above table are estimates
(formal CONTEXT deliveries are only provided with final builds).

## Unit Tests

### Setup

The test suite require access to a CRDS cache, but currently (2021-02-09) the shared /grp/crds cache does not include
Roman files. Developers inside the STScI network can sync a cache from roman-crds-test.stsci.edu (if working from home,
be sure to connect to the VPN first):

```bash
$ export CRDS_SERVER_URL=https://roman-crds-test.stsci.edu
$ export CRDS_PATH=$HOME/roman-crds-test-cache
$ crds sync --contexts roman-edit
```

The CRDS_READONLY_CACHE variable should not be set, since references will need to be downloaded to your local cache as
they are requested.

### Running tests

Unit tests can be run via `pytest`. Within the top level of your local `roman` repo checkout:

    pip install -e .[test]
    pytest

Need to parallelize your test runs over 8 cores?

    pip install pytest-xdist
    pytest -n 8

## Regression Tests

Latest regression test results can be found here (STScI staff only):

https://plwishmaster.stsci.edu:8081/job/RT/job/romancal/

To run the regression tests on your local machine, get the test dependencies and set the environment variable
TEST_BIGDATA to our Artifactory server
(STSci staff members only):

    pip install -e ".[test]"
    export TEST_BIGDATA=https://bytesalad.stsci.edu/artifactory

To run all the regression tests (except the very slow ones):

    pytest --bigdata romancal/regtest

You can control where the test results are written with the
`--basetemp=<PATH>` arg to `pytest`.  _NOTE that `pytest` will wipe this directory clean for each test session, so make
sure it is a scratch area._

If you would like to run a specific test, find its name or ID and use the `-k` option:

    pytest --bigdata romancal/regtest -k test_flat

If developers need to update the truth files in our nightly regression tests, there are instructions in this wiki.

https://github.com/spacetelescope/jwst/wiki/Maintaining-Regression-Tests

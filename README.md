# Roman Calibration Pipeline


[![Documentation Status](https://readthedocs.org/projects/roman-cal-pipeline/badge/?version=latest)](https://roman-cal-pipeline.readthedocs.io/en/latest/?badge=latest)
![Roman CI](https://github.com/spacetelescope/romancal/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/spacetelescope/romancal/branch/main/graph/badge.svg?token=S6KW6J7FZP)](https://codecov.io/gh/spacetelescope/romancal)
[![Powered by STScI Badge](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)
[![Powered by Astropy Badge](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

![STScI Logo](docs/_static/stsci_logo.png)

**Roman requires Python 3.6 or above and a C compiler for dependencies.**

**Linux and MacOS platforms are tested and supported.  Windows is not currently supported.**


## Installation

The easiest way to install the latest `roman` release into a fresh virtualenv or conda environment is

    pip install romancal

### Detailed Installation

The `romancal` package can be installed into a virtualenv or conda environment via `pip`.
We recommend that for each installation you start by creating a fresh
environment that only has Python installed and then install the `romancal` package and
its dependencies into that bare environment.
If using conda environments, first make sure you have a recent version of Anaconda
or Miniconda installed.
If desired, you can create multiple environments to allow for switching between different
versions of the `romancal` package (e.g. a released version versus the current development version).

In all cases, the installation is generally a 3-step process:
* Create a conda environment
* Activate that environment
* Install the desired version of the `romancal` package into that environment

Details are given below on how to do this for different types of installations,
including tagged releases, DMS builds used in operations, and development versions.
Remember that all conda operations must be done from within a bash shell.


### Installing latest releases

You can install the latest released version via `pip`.  From a bash shell:

    conda create -n <env_name> python
    conda activate <env_name>
    pip install romancal

You can also install a specific version (from `romancal 0.17.0` onward):

    conda create -n <env_name> python
    conda activate <env_name>
    pip install romancal==0.17.1

Installing specific versions before `romancal 0.17.0` need to be installed from Github:

    conda create -n <env_name> python
    conda activate <env_name>
    pip install git+https://github.com/spacetelescope/romancal@0.16.2


### Installing the development version from Github

You can install the latest development version (not as well tested) from the
Github master branch:

    conda create -n <env_name> python
    conda activate <env_name>
    pip install git+https://github.com/spacetelescope/romancal


### Installing a DMS Operational Build

There may be occasions where an exact copy of an operational DMS build is
desired (e.g. for validation testing or debugging operational issues).
We package releases for DMS builds via environment snapshots that specify the
exact versions of all packages to be installed. This method may result in more
stable processing than what was outlined above for installing a particular
tagged release, because that method installs the latest versions of dependency
packages, while this method installs dependencies pinned to particular versions
that have been well tested.

To install a particular DMS build, consult the
[Software vs DMS build version map](https://github.com/spacetelescope/romancal#software-vs-dms-build-version-map)
table shown below to determine the correct roman tag. For example, to install the
version of `romancal` used in DMS build 7.5, use romancal tag 0.16.1. The overall
procedure is similar to the 3-step process outlined in the previous section, but the
details of each command vary, due to the use of environment snapshot files that specify
all of the particular packages to install. Also note that different snapshot files are
used for Linux and Mac OS systems.

Linux:

    conda create -n <env_name> --file https://ssb.stsci.edu/releases/romandp/0.16.1/conda_python_stable-deps.txt
    conda activate <env_name>
    pip install -r https://ssb.stsci.edu/releases/romandp/0.16.1/reqs_stable-deps.txt

MacOS:

    conda create -n <env_name> --file https://ssb.stsci.edu/releases/romandp/0.16.1/conda_python_macos-stable-deps.txt
    conda activate <env_name>
    pip install -r https://ssb.stsci.edu/releases/romandp/0.16.1/reqs_macos-stable-deps.txt

Each DMS delivery has its own installation instructions, which may be found in
the corresponding release documentation linked from this page:
https://github.com/astroconda/astroconda-releases/tree/master/romandp
The installation procedures may change from time to time, so consulting the
documentation page for the specific version in question is the best way to get
that version installed.


### Installing for Developers

If you want to be able to work on and test the source code with the `romancal` package,
the high-level procedure to do this is to first create a conda environment using
the same procedures outlined above, but then install your personal copy of the
code overtop of the original code in that environment. Again, this should be done
in a separate conda environment from any existing environments that you may have
already installed with released versions of the `romancal` package.

As usual, the first two steps are to create and activate an environment:

    conda create -n <env_name> python
    conda activate <env_name>

To install your own copy of the code into that environment, you first need to
fork and clone the `romancal` repo:

    cd <where you want to put the repo>
    git clone https://github.com/spacetelescope/romancal
    cd romancal

*Note: `python setup.py install` and `python setup.py develop` commands do not work.*

Install from your local checked-out copy as an "editable" install:

    pip install -e .

If you want to run the unit or regression tests and/or build the docs, you can make
sure those dependencies are installed too:

    pip install -e .[test]
    pip install -e .[docs]
    pip install -e .[test,docs]

Need other useful packages in your development environment?

    pip install ipython pytest-xdist


## Calibration References Data System (CRDS) Setup

CRDS is the system that manages the reference files needed to run the pipeline.
Inside the STScI network, the pipeline works with default CRDS setup with no modifications.
To run the pipeline outside the STScI network, CRDS must be configured by setting
two environment variables:

    export CRDS_PATH=$HOME/crds_cache
    export CRDS_SERVER_URL=https://roman-crds.stsci.edu


## Documentation

Documentation (built daily from the Github `master` branch) is available at:

https://roman-pipeline.readthedocs.io/en/latest/

To build the docs yourself, clone this repository and build the documentation with:

    pip install -e .[docs]
    cd docs
    make html
    make latexpdf


## Contributions and Feedback

We welcome contributions and feedback on the project. Please follow the
[contributing guidelines](CONTRIBUTING.md) to submit an issue or a pull request.

We strive to provide a welcoming community to all of our users by abiding with
the [Code of Conduct](CODE_OF_CONDUCT.md).

If you have questions or concerns regarding the software, please open an issue
at https://github.com/spacetelescope/roman/issues or
contact the [ROMAN Help Desk](https://romanhelp.stsci.edu).


## Software vs DMS build version map

| roman tag | DMS build | CRDS_CONTEXT |   Date     |          Notes                                |
| -------- | --------- | ------------ | ---------- | ----------------------------------------------|
| 0.1.0    | 0.0       |  003         | Nov  2020  | Release for Build 0.0

Note: CRDS_CONTEXT values flagged with an asterisk in the above table are estimates
(formal CONTEXT deliveries are only provided with final builds).


## Unit Tests

### Setup

The test suite require access to a CRDS cache, but currently (2021-02-09) the shared /grp/crds
cache does not include Roman files.  Developers inside the STScI network can sync a cache from
roman-crds-test.stsci.edu (if working from home, be sure to connect to the VPN first):

```bash
$ export CRDS_SERVER_URL=https://roman-crds-test.stsci.edu
$ export CRDS_PATH=$HOME/roman-crds-test-cache
$ crds sync --contexts roman-edit
```

The CRDS_READONLY_CACHE variable should not be set, since references will need to be downloaded
to your local cache as they are requested.

### Running tests

Unit tests can be run via `pytest`.  Within the top level of your local `roman` repo checkout:

    pip install -e .[test]
    pytest

Need to parallelize your test runs over 8 cores?

    pip install pytest-xdist
    pytest -n 8



## Regression Tests

TBD

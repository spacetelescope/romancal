Installation
============
.. warning::
    Roman requires Python 3.9 or above and a C compiler for dependencies.

.. warning::
    Linux and MacOS platforms are tested and supported. Windows is not currently supported.

Stable releases of the ``romancal`` package are registered at
`PyPI <https://pypi.org/project/romancal/>`_. The development version of `romancal` is
installable from the
`Github repository <https://github.com/spacetelescope/romancal>`_.

The basic method of installing the roman calibration pipeline is to setup your environment and
issue the command,
::

    $ pip install romancal

Detailed Installation
---------------------

The `romancal` package can be installed into a virtualenv or conda environment via `pip`. We recommend that for each
installation you start by creating a fresh environment that only has Python installed and then install the `romancal`
package and its dependencies into that bare environment. If using conda environments, first make sure you have a recent
version of Anaconda or Miniconda installed. If desired, you can create multiple environments to allow for switching
between different versions of the `romancal` package (e.g. a released version versus the current development version).

In all cases, the recommended installation is generally a 3-step process:

- create a virtual environment;
- activate that environment;
- install the desired version of the `romancal` package into that environment.

Details are given below on how to do this for different types of installations, including tagged releases, DMS builds
used in operations, and development versions. Note that, although you can use any python environment management system that you are familiar with,
below we will be using `conda` (remember that all conda operations must be done from within a bash shell).

Installing latest releases
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can install the latest released version via `pip`. From a bash shell:
::

    $ conda create -n <env_name> python
    $ conda activate <env_name>
    $ pip install romancal

You can also install a specific version (from `romancal 0.1.0` onward):
::

    $ conda create -n <env_name> python
    $ conda activate <env_name>
    $ pip install romancal==0.5.0

Installing the development version from Github
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can install the latest development version (not as well tested) from the Github main branch:
::

    $ conda create -n <env_name> python
    $ conda activate <env_name>
    $ pip install git+https://github.com/spacetelescope/romancal


Installing for Developers
^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to be able to work on and test the source code with the `romancal` package, the high-level procedure to do
this is to first create a conda environment using the same procedures outlined above, but then install your personal
copy of the code overtop of the original code in that environment. Again, this should be done in a separate conda
environment from any existing environments that you may have already installed with released versions of the `romancal`
package.

As usual, the first two steps are to create and activate an environment:
::

    $ conda create -n <env_name> python
    $ conda activate <env_name>

To install your own copy of the code into that environment, you first need to fork and clone the `romancal` repo:
::

    $ cd <where you want to put the repo>
    $ git clone https://github.com/spacetelescope/romancal
    $ cd romancal

.. note::
    `python setup.py install` and `python setup.py develop` commands do not work.

Install from your local checked-out copy as an "editable" install:
::

    $ pip install -e .

If you want to run the unit or regression tests and/or build the docs, you can make sure those dependencies are
installed as well:
::

    $ pip install -e '.[test]'
    $ pip install -e '.[docs]'
    $ pip install -e '.[test,docs]'

Need other useful packages in your development environment?
::

    $ pip install ipython pytest-xdist

Calibration References Data System (CRDS) Setup
-----------------------------------------------

CRDS is the system that manages the reference files needed to run the pipeline. Inside the STScI network, the pipeline
works with default CRDS setup with no modifications. To run the pipeline outside the STScI network, CRDS must be
configured by setting two environment variables:
::

    $ export CRDS_PATH=$HOME/crds_cache
    $ export CRDS_SERVER_URL=https://roman-crds.stsci.edu

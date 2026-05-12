.. _romancal-version-scheme:

Versioning Scheme
-----------------

.. note::

   This section addresses the release version to PyPI.
   It is not to be confused with the DMS quarterly build version,
   which has a different numbering scheme (see
   `software vs. DMS builds <https://github.com/spacetelescope/romancal/blob/main/README.md#software-vs-dms-build-version-map>`_).

The ``romancal`` package follows `semantic versioning <https://semver.org/>`_ with a few minor
exceptions noted below. In brief this means that backwards incompatible changes are allowed
in major version changes, minor versions can contain new features and patch versions
can contain only bug fixes.

.. _romancal-public-vs-private-api:

API: Public vs Private
----------------------

As per Python convention, any API name that starts with underscore
(e.g., ``_my_private_function``) is considered private.
Any API not officially documented (i.e., you only found it after some extensive
code-diving) is also considered private.

Additionally test code is considered private. This includes:

* all ``conftest.py`` files
* modules that start with ``test_*`` or are named ``tests``
* everything under ``romancal.regtest``

If there is code that you would like to be public, please search
the open `issues <https://github.com/spacetelescope/romancal/issues>`_ and
open a new one describing what you would like to be made public.

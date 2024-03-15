"""Base classes for Roman tests"""

import os
import sys
from glob import glob as _sys_glob
from pathlib import Path

import requests

# Define location of default Artifactory API key, for Jenkins use only
ARTIFACTORY_API_KEY_FILE = "/eng/ssb2/keys/svc_rodata.key"


# Define a request timeout in seconds
TIMEOUT = 30


def _data_glob_local(*glob_parts):
    """Perform a glob on the local path

    Parameters
    ----------
    glob_parts: (path-like,[...])
        List of components that will be built into a single path

    Returns
    -------
    file_paths: [str[, ...]]
        Full file paths that match the glob criterion
    """
    full_glob = Path().joinpath(*glob_parts)
    return _sys_glob(str(full_glob))


def _data_glob_url(*url_parts, root=None):
    """
    Parameters
    ----------
    url: (str[,...])
        List of components that will be used to create a URL path

    root: str
        The root server path to the Artifactory server.
        Normally retrieved from `get_bigdata_root`.

    Returns
    -------
    url_paths: [str[, ...]]
        Full URLS that match the glob criterion
    """
    # Fix root root-ed-ness
    if root.endswith("/"):
        root = root[:-1]

    # Access
    try:
        envkey = os.environ["API_KEY_FILE"]
    except KeyError:
        envkey = ARTIFACTORY_API_KEY_FILE

    try:
        with open(envkey) as fp:
            headers = {"X-JFrog-Art-Api": fp.readline().strip()}
    except (PermissionError, FileNotFoundError):
        print(
            "Warning: Anonymous Artifactory search requests are limited "
            "to 1000 results. Use an API key and define API_KEY_FILE "
            "environment variable to get full search results.",
            file=sys.stderr,
        )
        headers = None

    search_url = "/".join([root, "api/search/pattern"])

    # Join and re-split the url so that every component is identified.
    url = "/".join([root] + [idx for idx in url_parts])
    all_parts = url.split("/")

    # Pick out "roman-pipeline", the repo name
    repo = all_parts[4]

    # Format the pattern
    pattern = repo + ":" + "/".join(all_parts[5:])

    # Make the query
    params = {"pattern": pattern}
    with requests.get(search_url, params=params, headers=headers, timeout=TIMEOUT) as r:
        url_paths = r.json()["files"]

    return url_paths

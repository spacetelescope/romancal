"""
Utilities
"""
import importlib.util
from importlib import import_module
import inspect
import logging
import os
import re

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Step classes that are not user-api steps
NON_STEPS = [
    'EngDBLogStep',
    'FunctionWrapper',
    'RomancalPipeline',
    'RomanPipeline',
    'ExposurePipeline',
    'RomanStep',
    'Pipeline',
    'Step',
    'SystemCall',
]


def all_steps():
    """List all classes subclassed from Step

    Returns
    -------
    steps : dict
        Key is the classname, value is the class
    """
    from romancal.stpipe import RomanStep as Step

    romancal = import_module('romancal')
    romancal_fpath = os.path.split(romancal.__file__)[0]

    steps = {}
    for module in load_local_pkg(romancal_fpath):
        more_steps = {
            klass_name: klass
            for klass_name, klass in inspect.getmembers(
                module,
                lambda o: inspect.isclass(o) and issubclass(o, Step)
            )
            if klass_name not in NON_STEPS
        }
        steps.update(more_steps)

    return steps


def load_local_pkg(fpath):
    """Generator producing all modules under fpath

    Parameters
    ----------
    fpath: string
        File path to the package to load.

    Returns
    -------
    generator
        `module` for each module found in the package.
    """
    package_fpath, package = os.path.split(fpath)
    package_fpath_len = len(package_fpath) + 1
    try:
        for module_fpath in folder_traverse(
            fpath, basename_regex=r'[^_].+\.py$', path_exclude_regex='tests'
        ):
            folder_path, fname = os.path.split(module_fpath[package_fpath_len:])
            module_path = folder_path.split('/')
            module_path.append(os.path.splitext(fname)[0])
            module_path = '.'.join(module_path)
            try:
                spec = importlib.util.spec_from_file_location(module_path,
                                                              module_fpath)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
            except Exception as err:
                logger.debug(f'Cannot load module "{module_path}": {str(err)}')
            else:
                yield module
    except Exception as err:
        logger.debug(f'Cannot complete package loading: Exception occurred: "{str(err)}"')


def folder_traverse(folder_path, basename_regex='.+', path_exclude_regex='^$'):
    """Generator of full file paths for all files
    in a folder.

    Parameters
    ----------
    folder_path: str
        The folder to traverse

    basename_regex: str
        Regular expression that must match
        the `basename` part of the file path.

    path_exclude_regex: str
        Regular expression to exclude a path.

    Returns
    -------
    generator
        A generator, return the next file.
    """
    basename_regex = re.compile(basename_regex)
    path_exclude_regex = re.compile(path_exclude_regex)
    for root, dirs, files in os.walk(folder_path):
        if path_exclude_regex.search(root):
            continue
        for file in files:
            if basename_regex.match(file):
                yield os.path.join(root, file)

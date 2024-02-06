"""
Utilities
"""

import inspect
import logging
from importlib import import_module
from pkgutil import walk_packages

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Step classes that are not user-api steps
NON_STEPS = [
    "EngDBLogStep",
    "FunctionWrapper",
    "RomancalPipeline",
    "RomanPipeline",
    "ExposurePipeline",
    "HighLevelPipeline",
    "RomanStep",
    "Pipeline",
    "Step",
    "SystemCall",
]


def all_steps():
    """List all classes subclassed from Step

    Returns
    -------
    steps : dict
        Key is the classname, value is the class
    """
    from romancal.stpipe import RomanStep as Step

    romancal = import_module("romancal")

    steps = {}
    for module in load_sub_modules(romancal):
        more_steps = {
            klass_name: klass
            for klass_name, klass in inspect.getmembers(
                module, lambda o: inspect.isclass(o) and issubclass(o, Step)
            )
            if klass_name not in NON_STEPS
        }
        steps.update(more_steps)

    return steps


def load_sub_modules(module):
    """
    Recursively loads all submodules of a module (this is not a local import).

    Parameters
    ----------
    module : module
        A python module to walk, load

    Returns
    -------
    generator
        A generator of all submodules of module recursively until no more sub modules are found
    """

    for package_info in walk_packages(module.__path__):
        if package_info.module_finder.path.startswith(module.__path__[0]):
            package = import_module(f"{module.__name__}.{package_info.name}")

            if package_info.ispkg:
                yield from load_sub_modules(package)

            yield package

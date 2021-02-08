from os.path import basename
from setuptools import setup, find_packages
from glob import glob


NAME = 'romancal'

SCRIPTS = [s for s in glob('scripts/*') if basename(s) != '__pycache__']

PACKAGE_DATA = {
    '': [
        '*.fits',
        '*.txt',
        '*.inc',
        '*.cfg',
        '*.csv',
        '*.yaml',
        '*.json',
        '*.asdf'
    ]
}

setup(
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    scripts=SCRIPTS,
    packages=find_packages(),
    package_data=PACKAGE_DATA,
)

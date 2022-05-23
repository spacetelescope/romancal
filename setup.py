from pathlib import Path

from setuptools import setup

SCRIPTS = [str(s) for s in Path('scripts').iterdir() if s.name != '__pycache__' and s.is_file()]

setup(scripts=SCRIPTS)

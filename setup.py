import itertools
import re
import os

from setuptools import find_namespace_packages, setup

dependencies = ["biopython", "marshmallow_dataclass[enum,union]", "marshmallow", "methodtools"]

with open(os.path.join(os.path.dirname(__file__), "inscripta", "biocantor", "__init__.py")) as v_file:
    VERSION = re.compile(r""".*__version__ = ["'](.*?)['"]""", re.S).match(v_file.read()).group(1)

extra_dependencies = {
    "io": ["gffutils"],
    "test": ["black", "flake8", "pytest", "pytest-cov", "pytest-error-for-skips"],
    "docs": [
        "Sphinx",
        "sphinx_rtd_theme",
        "sphinx-autoapi",
        "graphviz",
        "nbconvert",
        "nbsphinx",
        "ipykernel",
        "pandoc",
        "recommonmark",
    ],
    "optional": ["tornado>=5.1"],
}

all_dependencies = list(itertools.chain.from_iterable(extra_dependencies.values()))
extra_dependencies["all"] = all_dependencies

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="BioCantor",
    description="Flexible feature arithmetic, seamlessly integrated with nested coordinate systems.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Inscripta, Inc.",
    url="https://github.com/InscriptaLabs/BioCantor",
    test_suite="pytest",
    packages=find_namespace_packages(include=["inscripta.*"]),
    include_package_data=True,
    tests_require=extra_dependencies["test"],
    extras_require=extra_dependencies,
    install_requires=dependencies,
    version=VERSION,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

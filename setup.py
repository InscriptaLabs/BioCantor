import itertools

from setuptools import find_namespace_packages, setup

import versioneer

dependencies = ["versioneer", "biopython", "marshmallow_dataclass[enum,union]", "marshmallow", "methodtools"]

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

setup(
    name="BioCantor",
    long_description="BioCantor biocantor arithmetic library",
    description="BioCantor",
    author="Inscripta, Inc.",
    test_suite="pytest",
    packages=find_namespace_packages(include=["inscripta.*"]),
    include_package_data=True,
    tests_require=extra_dependencies["test"],
    extras_require=extra_dependencies,
    install_requires=dependencies,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)

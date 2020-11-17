# BioCantor

## Overview

Flexible feature arithmetic, seamlessly integrated with nested coordinate systems.

## What does the name BioCantor mean?

[Georg Cantor](https://en.wikipedia.org/wiki/Georg_Cantor) (1845-1918) 
was a German mathematician who created set theory. BioCantor uses set 
theoretic concepts, though Cantor's work extended into much more 
interesting realms than are used here.

Cantor devised the elegant [diagonal argument](https://en.wikipedia.org/wiki/Cantor%27s_diagonal_argument) 
proving that there is more than one level of infinity.

Also check out the [Cantor set](https://en.wikipedia.org/wiki/Cantor_set)
, which could exist as a genomic feature in a universe where the 
building blocks of genetic material are uncountably infinite.

## Installation

BioCantor can be added to a Python environment as follows.

Latest release from GitHub:

```
pip install git+https://github.com/InscriptaLabs/BioCantor
```

From source code:

```
cd BioCantor
pip install .
```
## Documentation

Full package documentation is available on [Read The Docs](https://biocantor.readthedocs.io/en/latest/).

#### To build documentation HTML pages locally

Install as above but with `libraries` and `docs` extras:

```
pip install -e .[libraries,docs]
```

You will also need to have [pandoc](https://pandoc.org/) installed. This can
be installed using `conda` with

```
conda install -y pandoc
```

You can now build the docs with:

```
cd docs
make html
```

## Support

Bug reports, support requests and feature requests should be submitted 
as issues on GitHub. We will make a reasonable effort to address issues 
dependent on available resources.

## Contributing

Contributors implicitly agree to the terms in [Inscripta_Contributor_License_Agreement.pdf](https://github.com/InscriptaLabs/biocantor/blob/master/Inscripta_Contributor_License_Agreement.pdf).

Users are invited to submit pull requests. Maintainers reserve the right 
to close PRs if the requested change is deemed inappropriate for the 
library or if resources are not available to sufficiently review the PR.

Submissions should follow the following guidelines:

- Code must adhere to [PEP 8](https://www.python.org/dev/peps/pep-0008/) 
(as enforced by [Flake8](https://flake8.pycqa.org/en/latest/)) and 
[Black](https://pypi.org/project/black/) style conventions with maximum 
line length of 120.
- New public elements must be documented with [NumPy style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html)
- Comprehensive unit tests must be added for new features.
- Full test suite must pass. It is recommended that you run the full test suite locally before pushing to GitHub.
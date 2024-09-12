# variant-filtration-tool

This repository contains the source code used in the VCF variant filtration workflows within the GDC. A single CLI is generated with multiple subcommands.

## Installation

```sh
pip install .
```

## Development

* Clone this repository
* Requirements:
  * Python >= 3.9
  * Tox
* `make venv` to create a virtualenv
* `source .venv/bin/activate` to activate new virtualenv
* `make init` to install dependencies and pre-commit hooks

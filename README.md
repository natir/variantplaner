# VariantPlanner

[![ci](https://github.com/natir/variantplanner/workflows/ci/badge.svg)](https://github.com/natir/variantplanner/actions?query=workflow%3Aci)
[![documentation](https://img.shields.io/badge/docs-mkdocs%20material-blue.svg?style=flat)](https://natir.github.io/variantplanner/)
[![pypi version](https://img.shields.io/pypi/v/variantplanner.svg)](https://pypi.org/project/variantplanner/)
[![gitpod](https://img.shields.io/badge/gitpod-workspace-blue.svg?style=flat)](https://gitpod.io/#https://github.com/natir/variantplanner)
[![gitter](https://badges.gitter.im/join%20chat.svg)](https://gitter.im/variantplanner/community)

A tool kit to manage many variant on desktop computer

## Installation

With `pip`:
```bash
pip install variantplanner
```

With [`pipx`](https://github.com/pipxproject/pipx):
```bash
python3.8 -m pip install --user pipx
pipx install variantplanner
```

## Developement setup

Initialisation step:

```
git clone git@github.com:natir/VariantPlanner.git
cd VariantPlanner
pyenv install 3.9 3.10 3.11 # take care of solving build warning
make setup
```

Working step:

```
<write code>
make format  # to auto-format the code

<write tests>
make test  # to run the test suite

make check  # to check if everything is OK

<commit your changes>
```

Please try to follow this commit message convention:

```
<type>(<scope>): <short summary>
  │       │             │
  │       │             └─⫸ Summary in present tense. Not capitalized. No period at the end.
  │       │
  │       └─⫸ Commit Scope: io|format_conversion|data_manipulation|(other)
  │
  └─⫸ Commit Type: build|ci|docs|feat|fix|perf|refactor|test|(other)
```

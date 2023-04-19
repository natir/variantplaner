"""All custom exception could be generate by VariantPlanner."""

# std import
from __future__ import annotations

import typing

if typing.TYPE_CHECKING:
    import pathlib

# 3rd party import

# project import


class NotAVCFError(Exception):
    """Exception raise if file read seems not be a vcf."""

    def __init__(self, path: pathlib.Path):
        """Initilize not a vcf error."""
        super().__init__(f"File {path} seems not be a vcf file.")


class NoGenotypeError(Exception):
    """Exception raise if file read seems not be a vcf."""

    def __init__(self):
        """Intilize no genotype error."""
        super().__init__(f"LazyFrame seems not contains genotypes.")

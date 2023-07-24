"""Exception could be generate by VariantPlanner."""

# std import
from __future__ import annotations

import typing

if typing.TYPE_CHECKING:  # pragma: no cover
    import pathlib

# 3rd party import

# project import


class NotAVCFError(Exception):
    """Exception raise if file read seems not be a vcf, generally not contains a line starts with '#CHROM'."""

    def __init__(self, path: pathlib.Path):
        """Initialize not a vcf error."""
        super().__init__(f"File {path} seems not be a valid vcf file.")


class NoGenotypeError(Exception):
    """Exception raise if vcf file seems not contains genotypes information."""

    def __init__(self):
        """Initialize no genotype error."""
        super().__init__("LazyFrame seems not contains genotypes.")


class NoGTError(Exception):
    """Exception raise if genotype polars.LazyFrame not contains gt column."""

    def __init__(self, message: str):
        """Initialize no gt error."""
        super().__init__(f"In {message} gt column is missing.")

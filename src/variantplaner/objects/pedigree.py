"""Declare Pedigree object."""

# std import
from __future__ import annotations

# 3rd party import
import polars

# project import


class Pedigree(polars.LazyFrame):
    """Object to manage lazyframe as Variants."""

    def __init__(self):
        """Initialize a Variants object."""
        self.lf = polars.LazyFrame(schema=Pedigree.minimal_schema())

    @classmethod
    def minimal_schema(cls) -> dict[str, type]:
        """Get schema of variants polars.LazyFrame."""
        return {
            "family_id": polars.String,
            "personal_id": polars.String,
            "father_id": polars.String,
            "mother_id": polars.String,
            "sex": polars.String,
            "affected": polars.Boolean,
        }

"""Declare Genotypes object."""

# std import
from __future__ import annotations

# 3rd party import
import polars

# project import


class Annotations(polars.LazyFrame):
    """Object to manage lazyframe as Annotations."""

    def __init__(self):
        """Initialize a Annotations object."""
        self.lf = polars.LazyFrame(schema=Annotations.minimal_schema())

    @classmethod
    def minimal_schema(cls) -> dict[str, type]:
        """Get minimal schema of genotypes polars.LazyFrame."""
        return {
            "id": polars.UInt64,
        }

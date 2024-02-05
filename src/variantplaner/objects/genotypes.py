"""Declare Genotypes object."""

# std import
from __future__ import annotations

# 3rd party import
import polars

# project import


class Genotypes(polars.LazyFrame):
    """Object to manage lazyframe as Genotypes."""

    def __init__(self):
        """Initialize a Genotypes object."""
        self.lf = polars.LazyFrame(schema=Genotypes.minimal_schema())

    @classmethod
    def minimal_schema(cls) -> dict[str, type]:
        """Get minimal schema of genotypes polars.LazyFrame."""
        return {
            "id": polars.UInt64,
            "samples": polars.String,
        }

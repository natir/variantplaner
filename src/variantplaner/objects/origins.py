"""Declare Origin object."""

# std import
from __future__ import annotations

# 3rd party import
import polars

# project import


class Origins(polars.LazyFrame):
    """Object to manage lazyframe as Origins."""

    def __init__(self):
        """Initialize a Origins object."""
        self.lf = polars.LazyFrame(schema=Origins.minimal_schema())

    @classmethod
    def minimal_schema(cls) -> dict[str, type]:
        """Get minimal schema of genotypes polars.LazyFrame."""
        return {
            "id": polars.UInt64,
            "index_gt": polars.UInt8,
            "mother_gt": polars.UInt8,
            "father_gt": polars.UInt8,
            "origin": polars.String,
        }

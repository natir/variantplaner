"""Declare Variants object."""

# std import
from __future__ import annotations

# 3rd party import
import polars

# project import


class Variants(polars.LazyFrame):
    """Object to manage lazyframe as Variants."""

    def __init__(self, data: polars.LazyFrame | None = None):
        """Initialize a Variants object."""
        if data is None:
            self.lf = polars.LazyFrame(schema=Variants.minimal_schema())
        else:
            self.lf = data

    @classmethod
    def minimal_schema(cls) -> dict[str, type]:
        """Get schema of variants polars.LazyFrame."""
        return {
            "id": polars.UInt64,
            "chr": polars.String,
            "pos": polars.UInt64,
            "ref": polars.String,
            "alt": polars.String,
        }

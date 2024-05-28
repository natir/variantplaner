"""Declare Genotypes object."""

# std import
from __future__ import annotations

# 3rd party import
import polars

# project import


class Genotypes(polars.LazyFrame):
    """Object to manage lazyframe as Genotypes."""

    def __init__(self, data: polars.LazyFrame | None = None):
        """Initialize a Genotypes object."""
        if data is None:
            self.lf = polars.LazyFrame(schema=Genotypes.minimal_schema())
        else:
            self.lf = data


    def samples_names(self) -> list[str]:
        """Get list of sample name."""
        return self.lf.select("sample").collect().get_column("sample").to_list()

    @classmethod
    def minimal_schema(cls) -> dict[str, type]:
        """Get minimal schema of genotypes polars.LazyFrame."""
        return {
            "id": polars.UInt64,
            "sample": polars.String,
        }

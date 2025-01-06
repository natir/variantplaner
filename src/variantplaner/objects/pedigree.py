"""Declare Pedigree object."""

# std import
from __future__ import annotations

import typing

if typing.TYPE_CHECKING:  # pragma: no cover
    import collections
    import pathlib

# 3rd party import
import polars


class Pedigree(polars.LazyFrame):
    """Object to manage lazyframe as Variants."""

    def __init__(self):
        """Initialize a Variants object."""
        self.lf = polars.LazyFrame(schema=Pedigree.minimal_schema())

    def from_path(self, input_path: pathlib.Path) -> None:
        """Read a pedigree file in [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html).

        Args:
            input_path: Path to pedigree file.

        Returns:
            A [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) that contains ped information ('family_id', 'personal_id', 'father_id', 'mother_id', 'sex', 'affected')
        """
        self.lf = polars.scan_csv(
            input_path,
            separator="\t",
            has_header=False,
            null_values=["None", "unknown"],
            new_columns=[
                "family_id",
                "personal_id",
                "father_id",
                "mother_id",
                "sex",
                "affected",
            ],
            schema_overrides=Pedigree.minimal_schema(),
        )

    def to_path(self, output_path: pathlib.Path) -> None:
        """Write pedigree [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) in ped format.

        Warning: This function performs [polars.LazyFrame.collect][] before write csv, this can have a significant impact on memory usage

        Args:
            lf: LazyFrame contains pedigree information.
            output_path: Path where write pedigree information.

        Returns:
            None
        """
        self.lf.collect().write_csv(output_path, include_header=False, separator="\t")

    @classmethod
    def minimal_schema(
        cls,
    ) -> collections.abc.Mapping[str, polars._typing.PolarsDataType]:
        """Get schema of variants polars.LazyFrame."""
        return {
            "family_id": polars.String,
            "personal_id": polars.String,
            "father_id": polars.String,
            "mother_id": polars.String,
            "sex": polars.String,
            "affected": polars.Boolean,
        }

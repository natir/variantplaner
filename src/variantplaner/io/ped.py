"""Read and write ped file."""

# std import
from __future__ import annotations

import typing

if typing.TYPE_CHECKING:  # pragma: no cover
    import pathlib

# 3rd party import
import polars

# project import


def into_lazyframe(input_path: pathlib.Path) -> polars.LazyFrame:
    """Read a pedigree file in [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html).

    Args:
        input_path: Path to pedigree file.

    Returns:
        A [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) that contains ped information ('family_id', 'personal_id', 'father_id', 'mother_id', 'sex', 'affected')
    """
    return polars.scan_csv(
        input_path,
        separator="\t",
        has_header=False,
        null_values="None",
        new_columns=["family_id", "personal_id", "father_id", "mother_id", "sex", "affected"],
        dtypes={
            "family_id": polars.Utf8,
            "personal_id": polars.Utf8,
            "father_id": polars.Utf8,
            "mother_id": polars.Utf8,
            "sex": polars.Utf8,
            "affected": polars.Boolean,
        },
    )


def from_lazyframe(lf: polars.LazyFrame, output_path: pathlib.Path) -> None:
    """Write pedigree [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) in ped format.

    Warning: This function performs [polars.LazyFrame.collect][] before write csv, this can have a significant impact on memory usage

    Args:
        lf: LazyFrame contains pedigree information.
        output_path: Path where write pedigree information.

    Returns:
        None
    """
    lf.collect().write_csv(output_path, include_header=False, separator="\t")

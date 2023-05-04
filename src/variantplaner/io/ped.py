"""Function to manage input and output of ped file."""

# std import
from __future__ import annotations

import typing

if typing.TYPE_CHECKING:  # pragma: no cover
    import pathlib
    import sys

    if sys.version_info >= (3, 11):
        from typing import ParamSpec
    else:
        from typing_extensions import ParamSpec

    T = typing.TypeVar("T")
    P = ParamSpec("P")


# 3rd party import
import polars

# project import


def into_lazyframe(input_path: pathlib.Path) -> polars.LazyFrame:
    """Read a pedigre file and convert it in lazyframe.

    Args:
        input_path: Path to pedigre file.

    Returns:
        A lazyframe that containt ped information ('family_id', 'personal_id', 'father_id', 'mother_id', 'sex', 'affected')
    """
    return polars.scan_csv(
        input_path,
        separator="\t",
        has_header=False,
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
    """Write pedigre lazyframe in ped format.

    Warning: This function perform LazyFrame.collect() before write csv, this can have a significant impact on memory usage

    Args:
        lf: LazyFrame contains pedigre information.
        output_path: Path where write pedigre information.

    Returns:
        None
    """
    lf.collect().write_csv(output_path, has_header=False, separator="\t")

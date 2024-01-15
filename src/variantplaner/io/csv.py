"""Read and write csv file."""

# std import
from __future__ import annotations

import typing

# 3rd party import
import polars

if typing.TYPE_CHECKING:  # pragma: no cover
    import pathlib
    import sys
    from collections.abc import Sequence

    if sys.version_info >= (3, 11):
        from typing import Unpack
    else:
        from typing_extensions import Unpack

    if sys.version_info >= (3, 11):

        class ScanCsv(typing.TypedDict, total=False):
            """A struct to check type of parameter give to [polars.scan_csv][]."""

            has_header: bool
            separator: str
            comment_prefix: str | None
            quote_char: str | None
            skip_rows: int
            dtypes: polars.type_aliases.SchemaDict | Sequence[polars.type_aliases.PolarsDataType] | None
            null_values: str | Sequence[str] | dict[str, str] | None
            missing_utf8_is_empty_string: bool
            ignore_errors: bool
            cache: bool
            with_column_names: typing.Callable[[list[str]], list[str]] | None
            infer_schema_length: int | None
            n_rows: int | None
            encoding: polars.type_aliases.CsvEncoding
            low_memory: bool
            rechunk: bool
            skip_rows_after_header: int
            row_index_name: str | None
            row_index_offset: int
            try_parse_dates: bool
            eol_char: str
            new_columns: Sequence[str] | None
    else:

        class ScanCsv(typing.TypedDict, total=False):
            """A struct to check type of parameter give to [polars.scan_csv][]."""

            has_header: bool
            separator: str
            comment_prefix: str | None
            quote_char: str | None
            skip_rows: int
            dtypes: polars.type_aliases.SchemaDict | Sequence[polars.type_aliases.PolarsDataType] | None
            null_values: str | Sequence[str] | dict[str, str] | None
            missing_utf8_is_empty_string: bool
            ignore_errors: bool
            cache: bool
            with_column_names: typing.Callable[[list[str]], list[str]] | None
            infer_schema_length: int | None
            n_rows: int | None
            encoding: polars.type_aliases.CsvEncoding
            low_memory: bool
            rechunk: bool
            skip_rows_after_header: int
            row_count_name: str | None
            row_count_offset: int
            try_parse_dates: bool
            eol_char: str
            new_columns: Sequence[str] | None


# project import
from variantplaner import normalization


def into_lazyframe(
    input_path: pathlib.Path,
    chr2len_path: pathlib.Path,
    chromosome_col: str,
    position_col: str,
    reference_col: str,
    alternative_col: str,
    info_cols: list[str],
    /,
    **scan_csv_args: Unpack[ScanCsv],
) -> polars.LazyFrame:
    """Read a csv file and convert it in lazyframe.

    Args:
        input_path: Path to csv.
        chr2len_path: Path to chr2length csv.
        chromosome_col: Name of the column that holds the chromosomes.
        position_col: Name of the column that holds the positions.
        reference_col: Name of the column that holds the reference sequence.
        alternative_col: Name of the column that hold the alternative sequence.
        scan_csv_args: polars.scan_csv parameter.

    Returns:
        A lazyframe that contain csv information

    """
    lf = polars.scan_csv(
        input_path,
        **scan_csv_args,
    )
    chr2len = chr2length_into_lazyframe(chr2len_path)

    lf = lf.rename(
        {
            chromosome_col: "chr",
            position_col: "pos",
            reference_col: "ref",
            alternative_col: "alt",
        },
    )

    lf = lf.cast({"pos": polars.UInt64})

    if info_cols:
        lf = lf.select(["chr", "pos", "ref", "alt", *info_cols])

    return normalization.add_variant_id(lf, chr2len)


def chr2length_into_lazyframe(input_path: pathlib.Path) -> polars.LazyFrame:
    """Read a csv file with two column chr and length and perform some percomputation.

    Args:
        input_path: Path to csv.

    Returns:
        A lazyframe with chromosome name associate to length, offset information
    """
    lf = polars.scan_csv(input_path, schema={"chr": polars.Utf8, "length": polars.UInt64})
    return lf.with_columns(
        offset=polars.col("length").cum_sum() - polars.col("length"),
    )

"""Declare CSV object."""

# std import
from __future__ import annotations

import dataclasses
import typing

# 3rd party import
import polars

# project import
from variantplaner.exception import NotAVariantCsvError

if typing.TYPE_CHECKING:  # pragma: no cover
    import pathlib
    import sys
    from collections.abc import Sequence

    if sys.version_info >= (3, 11):
        from typing import Unpack
    else:
        from typing_extensions import Unpack

    class ScanCsv(typing.TypedDict, total=False):
        """A struct to check type of parameter give to [polars.scan_csv][]."""

        has_header: bool
        separator: str
        comment_prefix: str | None
        quote_char: str | None
        skip_rows: int
        schema_overrides: polars._typing.SchemaDict | Sequence[polars._typing.PolarsDataType] | None
        null_values: str | Sequence[str] | dict[str, str] | None
        missing_utf8_is_empty_string: bool
        ignore_errors: bool
        cache: bool
        with_column_names: typing.Callable[[list[str]], list[str]] | None
        infer_schema_length: int | None
        n_rows: int | None
        encoding: polars._typing.CsvEncoding
        low_memory: bool
        rechunk: bool
        skip_rows_after_header: int
        row_index_name: str | None
        row_index_offset: int
        try_parse_dates: bool
        eol_char: str
        new_columns: Sequence[str] | None


@dataclasses.dataclass
class ColRename:
    """A struct to store rename parameter."""

    chr: str = "chr"
    ref: str = "ref"
    alt: str = "alt"
    other: dict[str, str] = dataclasses.field(default_factory=dict)


class Csv(polars.LazyFrame):
    """Object to manage lazyframe as Csv."""

    def __init__(self):
        """Initialize a Csv object."""
        self.lf = polars.LazyFrame()

    def from_path(self, path: pathlib.Path, /, **scan_csv_args: Unpack[ScanCsv]) -> None:
        """Populate Csv obejct with csv file content."""
        self.lf = polars.scan_csv(path, **scan_csv_args)

    def variants_from_path(
        self,
        path: pathlib.Path,
        col_rename: ColRename,
        /,
        **scan_csv_args: Unpack[ScanCsv],
    ) -> None:
        """Populate Csv object with csv file."""
        self.from_path(path, **scan_csv_args)

        self.lf = self.lf.rename(dataclasses.asdict(col_rename))

        if any(elt not in super().columns for elt in ["chr", "pos", "ref", "alt"]):
            raise NotAVariantCsvError(path)

        self.lf = self.lf.cast({"pos": polars.UInt64})

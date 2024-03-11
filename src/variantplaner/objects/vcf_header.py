"""Declare Vcf object."""

# std import
from __future__ import annotations

import functools
import re
import typing

# 3rd party import
import polars

# project import
from variantplaner.exception import NotVcfHeaderError

# type checking block
if typing.TYPE_CHECKING:
    import pathlib

MINIMAL_COL_NUMBER: int = 8
SAMPLE_COL_BEGIN: int = 9


class VcfHeader:
    """Object that parse and store vcf information."""

    def __init__(self):
        """Initialise VcfHeader."""
        self._header = []

    def from_files(self, path: pathlib.Path) -> None:
        """Populate VcfHeader object with content of only header file.

        Args:
        path: Path of file

        Returns:
        None
        """
        with open(path) as fh:
            for full_line in fh:
                line = full_line.strip()
                self._header.append(line)

    def from_lines(self, lines: typing.Iterator[str]) -> None:
        """Extract all header information of vcf lines.

        Line between start of file and first line start with '#CHROM' or not start with '#'

        Args:
        lines: Iterator of line

        Returns: None

        Raises:
        NotAVcfHeader: If a line not starts with '#'
        NotAVcfHeader: If no line start by '#CHROM'
        """
        for full_line in lines:
            line = full_line.strip()

            if not line.startswith("#"):
                raise NotVcfHeaderError

            if line.startswith("#CHROM"):
                self._header.append(line)
                return

            self._header.append(line)

        raise NotVcfHeaderError

    def info_parser(self, select_info: set[str] | None = None) -> list[polars.Expr]:
        """Generate a list of [polars.Expr](https://pola-rs.github.io/polars/py-polars/html/reference/expressions/index.html) to extract variants information.

        Args:
        header: Line of vcf header
        input_path: Path to vcf file.
        select_info: List of target info field

        Returns:
        List of [polars.Expr](https://pola-rs.github.io/polars/py-polars/html/reference/expressions/index.html) to parse info columns.

        Raises:
        NotVcfHeaderError: If all line not start by '#CHR'
        """
        info_re = re.compile(
            r"ID=(?P<id>([A-Za-z_][0-9A-Za-z_.]*|1000G)),Number=(?P<number>[ARG0-9\.]+),Type=(?P<type>Integer|Float|String|Character)",
        )

        expressions: list[polars.Expr] = []

        for line in self._header:
            if line.startswith("#CHROM"):
                return expressions

            if not line.startswith("##INFO"):
                continue

            if (search := info_re.search(line)) and (not select_info or search["id"] in select_info):
                regex = rf"{search['id']}=([^;]+);?"

                local_expr = polars.col("info").str.extract(regex, 1)

                if search["number"] == "1":
                    if search["type"] == "Integer":
                        local_expr = local_expr.cast(polars.Int64)
                    elif search["type"] == "Float":
                        local_expr = local_expr.cast(polars.Float64)
                    elif search["type"] in {"String", "Character"}:
                        pass  # Not do anything on string or character
                    else:
                        pass  # Not reachable

                else:
                    local_expr = local_expr.str.split(",")
                    if search["type"] == "Integer":
                        local_expr = local_expr.cast(polars.List(polars.Int64))
                    elif search["type"] == "Float":
                        local_expr = local_expr.cast(polars.List(polars.Float64))
                    elif search["type"] in {"String", "Character"}:
                        pass  # Not do anything on string or character
                    else:
                        pass  # Not reachable

                expressions.append(local_expr.alias(search["id"]))

        raise NotVcfHeaderError

    def format_parser(
        self,
        select_format: set[str] | None = None,
    ) -> dict[str, typing.Callable[[polars.Expr, str], polars.Expr]]:
        """Generate a list of [polars.Expr](https://pola-rs.github.io/polars/py-polars/html/reference/expressions/index.html) to extract genotypes information.

        **Warning**: Float values can't be converted for the moment they are stored as String to keep information

        Args:
        header: Line of vcf header.
        input_path: Path to vcf file.
        select_format: List of target format field.

        Returns:
        A dict to link format id to pipeable function with Polars.Expr

        Raises:
        NotVcfHeaderError: If all line not start by '#CHR'
        """
        format_re = re.compile(
            "ID=(?P<id>[A-Za-z_][0-9A-Za-z_.]*),Number=(?P<number>[ARG0-9\\.]+),Type=(?P<type>Integer|Float|String|Character)",
        )

        expressions: dict[str, typing.Callable[[polars.Expr, str], polars.Expr]] = {}

        for line in self._header:
            if line.startswith("#CHROM"):
                return expressions

            if not line.startswith("##FORMAT"):
                continue

            if (search := format_re.search(line)) and (not select_format or search["id"] in select_format):
                name = search["id"]
                number = search["number"]
                format_type = search["type"]

                if name == "GT":
                    expressions["GT"] = VcfHeader.__format_gt
                    continue

                if number == "1":
                    if format_type == "Integer":
                        expressions[name] = VcfHeader.__format_one_int
                    elif format_type == "Float":  # noqa: SIM114 Float isn't already support but in future
                        expressions[name] = VcfHeader.__format_one_str
                    elif format_type in {"String", "Character"}:
                        expressions[name] = VcfHeader.__format_one_str
                    else:
                        pass  # Not reachable

                elif format_type == "Integer":
                    expressions[name] = VcfHeader.__format_list_int
                elif format_type == "Float":  # noqa: SIM114 Float isn't already support but in future
                    expressions[name] = VcfHeader.__format_list_str
                elif format_type in {"String", "Character"}:
                    expressions[name] = VcfHeader.__format_list_str
                else:
                    pass  # Not reachable

        raise NotVcfHeaderError

    @functools.cached_property
    def samples_index(self) -> dict[str, int] | None:
        """Read vcf header to generate an association map between sample name and index.

        Args:
        header: Header string.

        Returns:
        Map that associate a sample name to is sample index.

        Raises:
        NotVcfHeaderError: If all line not start by '#CHR'
        """
        for line in reversed(self._header):
            if line.startswith("#CHR"):
                split_line = line.strip().split("\t")
                if len(split_line) <= MINIMAL_COL_NUMBER:
                    return None

                return {sample: i for (i, sample) in enumerate(split_line[SAMPLE_COL_BEGIN:])}

        raise NotVcfHeaderError

    @functools.cached_property
    def contigs(self) -> typing.Iterator[str]:
        """Get an iterator of line contains chromosomes information.

        Returns: String iterator
        """
        for full_line in self._header:
            if full_line.startswith("##contig"):
                yield full_line.strip()

    def column_name(self, number_of_column: int = MINIMAL_COL_NUMBER) -> typing.Iterator[str]:
        """Get an iterator of correct column name.

        Returns: String iterator
        """
        base_col_name = ["chr", "pos", "vid", "ref", "alt", "qual", "filter", "info"]

        yield from base_col_name

        if number_of_column > MINIMAL_COL_NUMBER and (samples := self.samples_index):
            yield "format"
            yield from (sample for (sample, _) in samples.items())

    @staticmethod
    def __format_gt(expr: polars.Expr, /, col_name: str) -> polars.Expr:
        """Manage gt field."""
        return expr.str.count_matches("1").cast(polars.UInt8).alias(col_name.lower())

    @staticmethod
    def __format_one_int(expr: polars.Expr, /, col_name: str) -> polars.Expr:
        """Manage integer field."""
        return expr.str.to_integer(base=10, strict=False).cast(polars.UInt32).alias(col_name.lower())

    @staticmethod
    def __format_one_str(expr: polars.Expr, /, col_name: str) -> polars.Expr:
        """Manage string field."""
        return expr.alias(col_name.lower())

    @staticmethod
    def __format_list_int(expr: polars.Expr, /, col_name: str) -> polars.Expr:
        """Manage list of integer field."""
        return (
            expr.str.split(",")
            .list.eval(polars.element().str.to_integer(base=10, strict=False).cast(polars.UInt32))
            .alias(col_name.lower())
        )

    @staticmethod
    def __format_list_str(expr: polars.Expr, /, col_name: str) -> polars.Expr:
        """Manag list string field."""
        return expr.str.split(",").alias(col_name.lower())

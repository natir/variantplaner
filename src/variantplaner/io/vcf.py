"""Function to manage input and output of vcf file."""

# std import
from __future__ import annotations

import enum
import logging
import re
import typing
from typing import Callable

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
from variantplaner import normalization
from variantplaner.exception import NotAVCFError

MINIMAL_COL_NUMBER: int = 8
SAMPLE_COL_BEGIN: int = 9

logger = logging.getLogger("io.vcf")


class IntoLazyFrameExtension(enum.Enum):
    """Enumration use to control behavior of IntoLazyFrame."""

    NOTHING = 0
    MANAGE_SV = 1


def extract_header(input_path: pathlib.Path) -> list[str]:
    """Extract all header information of vcf file.

    Line between start of file and first line start with '#CHROM' or not start with '#'

    Args:
        input_path: Path to vcf file.

    Returns:
        List of header string.

    Raises:
        NotAVCFError: If a line not starts with '#'
        NotAVCFError: If all line not start by '#CHROM'
    """
    headers = []
    with open(input_path) as fh:
        for line in (full_line.strip() for full_line in fh):
            if not line.startswith("#"):
                raise NotAVCFError(input_path)

            if line.startswith("#CHROM"):
                headers.append(line)
                return headers

            headers.append(line)

    raise NotAVCFError(input_path)


def info2expr(header: list[str], input_path: pathlib.Path, select_info: set[str] | None = None) -> list[polars.Expr]:
    """Read vcf header to generate a list of polars.Expr to extract variants informations.

    Args:
        header: Line of vcf header
        input_path: Path to vcf file.
        select_info: List of target info field

    Returns:
        List of polars expr to parse info columns.

    Raises:
        NotAVCFError: If all line not start by '#CHR'
        NotSupportType: If header line indicate a not support type
    """
    info_re = re.compile(
        r"ID=(?P<id>([A-Za-z_][0-9A-Za-z_.]*|1000G)),Number=(?P<number>[ARG0-9\.]+),Type=(?P<type>Integer|Float|String|Character)",
    )

    expressions: list[polars.Expr] = []

    for line in header:
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
                elif search["type"] == "String" or search["type"] == "Character":
                    pass  # Not do anything on string or character
                else:
                    pass  # Not reachable

            else:
                local_expr = local_expr.str.split(",")
                if search["type"] == "Integer":
                    local_expr = local_expr.cast(polars.List(polars.Int64))
                elif search["type"] == "Float":
                    local_expr = local_expr.cast(polars.List(polars.Float64))
                elif search["type"] == "String" or search["type"] == "Character":
                    pass  # Not do anything on string or character
                else:
                    pass  # Not reachable

            expressions.append(local_expr.alias(search["id"]))

    raise NotAVCFError(input_path)


def __format_gt(expr: polars.Expr, /, name: str) -> polars.Expr:
    """Manage gt field."""
    return expr.str.count_match("1").cast(polars.UInt8).alias(name.lower())


def __format_one_int(expr: polars.Expr, /, name: str) -> polars.Expr:
    """Manage integer field."""
    return expr.str.parse_int(10, strict=False).cast(polars.UInt16).alias(name.lower())


def __format_one_str(expr: polars.Expr, /, name: str) -> polars.Expr:
    """Manage string field."""
    return expr.alias(name.lower())


def __format_list_int(expr: polars.Expr, /, name: str) -> polars.Expr:
    """Manage list of integer field."""
    return (
        expr.str.split(",")
        .list.eval(polars.element().str.parse_int(10, strict=False).cast(polars.UInt16))
        .alias(name.lower())
    )


def __format_list_str(expr: polars.Expr, /, name: str) -> polars.Expr:
    """Manag list string field."""
    return expr.str.split(",").alias(name.lower())


def format2expr(
    header: list[str],
    input_path: pathlib.Path,
    select_format: set[str] | None = None,
) -> dict[str, Callable[[polars.Expr, str], polars.Expr]]:
    """Read vcf header to generate a list of polars.Expr to extract genotypes informations.

    **Warning**: Float values can't be converted for the moment they are stored as String to keep information

    Args:
        header: Line of vcf header.
        input_path: Path to vcf file.
        select_format: List of target format field.

    Returns:
        A dict to link format id to pipeable function with Polars.Expr

    Raises:
        NotAVCFError: If all line not start by '#CHR'
    """
    format_re = re.compile(
        "ID=(?P<id>[A-Za-z_][0-9A-Za-z_.]*),Number=(?P<number>[ARG0-9\\.]+),Type=(?P<type>Integer|Float|String|Character)",
    )

    expressions: dict[str, Callable[[polars.Expr, str], polars.Expr]] = {}

    for line in header:
        if line.startswith("#CHROM"):
            return expressions

        if not line.startswith("##FORMAT"):
            continue

        if (search := format_re.search(line)) and (not select_format or search["id"] in select_format):
            name = search["id"]
            number = search["number"]
            format_type = search["type"]

            if name == "GT":
                expressions["GT"] = __format_gt
                continue

            if number == "1":
                if format_type == "Integer":
                    expressions[name] = __format_one_int
                elif format_type == "Float":  # noqa: SIM114 Float isn't already support but in future
                    expressions[name] = __format_one_str
                elif format_type == "String" or format_type == "Character":
                    expressions[name] = __format_one_str
                else:
                    pass  # Not reachable

            else:
                if format_type == "Integer":  # noqa: PLR5501 All other number are consider as list
                    expressions[name] = __format_list_int
                elif format_type == "Float":  # noqa: SIM114 Float isn't already support but in future
                    expressions[name] = __format_list_str
                elif format_type == "String" or format_type == "Character":
                    expressions[name] = __format_list_str
                else:
                    pass  # Not reachable

    raise NotAVCFError(input_path)


def sample_index(header: list[str], input_path: pathlib.Path) -> dict[str, int] | None:
    """Read vcf header to generate an association map between sample name and index.

    Args:
        header: Header string.

    Returns:
        Map that associate a sample name to is sample index.

    Raises:
        NotAVCFError: If all line not start by '#CHR'
    """
    for line in reversed(header):
        if line.startswith("#CHR"):
            split_line = line.strip().split("\t")
            if len(split_line) <= MINIMAL_COL_NUMBER:
                return None

            return {sample: i for (i, sample) in enumerate(split_line[SAMPLE_COL_BEGIN:])}

    raise NotAVCFError(input_path)


def __column_name(header: list[str], input_path: pathlib.Path) -> list[str]:
    """Read vcf header to generate list of column name.

    Args:
        header: Header string.

    Returns:
        List of lazyframe columns name.

    Raises:
        NotAVCFError: If all line not start by '#CHR'
    """
    for line in reversed(header):
        if line.startswith("#CHR"):
            split_line = line.strip().split("\t")
            cols_name = ["chr", "pos", "vid", "ref", "alt", "qual", "filter", "info"]
            if len(split_line) > MINIMAL_COL_NUMBER:
                cols_name.append("format")
                for sample in split_line[SAMPLE_COL_BEGIN:]:
                    cols_name.append(sample)

            return cols_name

    raise NotAVCFError(input_path)


def into_lazyframe(
    input_path: pathlib.Path,
    extension: IntoLazyFrameExtension = IntoLazyFrameExtension.NOTHING,
) -> polars.LazyFrame:
    """Read a vcf file and convert it in lazyframe.

    Args:
        input_path: Path to vcf file.

    Returns:
        A lazyframe that containt vcf information ('chr', 'pos', 'vid', 'ref', 'alt', 'qual', 'filter', 'info', ['format'], ['genotypes',…], 'id').
    """
    header = extract_header(input_path)

    col_name = {f"column_{i}": name for (i, name) in enumerate(__column_name(header, input_path), start=1)}

    lf = polars.scan_csv(
        input_path,
        separator="\t",
        comment_char="#",
        has_header=False,
        dtypes={"column_1": polars.Utf8},
        ignore_errors=True,
    )

    lf = lf.rename(col_name)

    if extension == IntoLazyFrameExtension.MANAGE_SV:
        lf = lf.with_columns(info2expr(header, input_path, {"SVTYPE", "SVLEN"}))

    lf = normalization.chromosome2integer(lf)

    lf = normalization.add_variant_id(lf)

    if extension == IntoLazyFrameExtension.MANAGE_SV:
        drop_column = {"SVTYPE", "SVLEN"}
        lf = lf.collect().select([col for col in lf.columns if col not in drop_column]).lazy()

    return lf


if typing.TYPE_CHECKING:  # pragma: no cover
    RenameCol = typing.TypedDict(
        "RenameCol",
        {
            "#CHROM": str,
            "POS": str,
            "ID": str,
            "REF": str,
            "ALT": str,
            "QUAL": str,
            "FILTER": str,
            "INFO": list[tuple[str, str]],
            "FORMAT": str,
            "sample": dict[str, dict[str, str]],
        },
    )

DEFAULT_RENAME: RenameCol = {
    "#CHROM": "chr",
    "POS": "pos",
    "ID": "id",
    "REF": "ref",
    "ALT": "alt",
    "QUAL": ".",
    "FILTER": ".",
    "INFO": [],
    "FORMAT": "",
    "sample": {},
}


def build_rename_column(
    chromosome: str,
    pos: str,
    identifier: str,
    ref: str,
    alt: str,
    qual: str | None = ".",
    filter_col: str | None = ".",
    info: list[tuple[str, str]] | None = None,
    format_string: str | None = None,
    sample: dict[str, dict[str, str]] | None = None,
) -> RenameCol:
    """An helper function to generate rename column dict for io.vcf.from_lazyframe function parameter.

    Returns:
        A rename column dictionary.
    """
    return {
        "#CHROM": chromosome,
        "POS": pos,
        "ID": identifier,
        "REF": ref,
        "ALT": alt,
        "QUAL": "." if qual is None else qual,
        "FILTER": "." if filter_col is None else filter_col,
        "INFO": [] if info is None else info,
        "FORMAT": "" if format_string is None else format_string,
        "sample": {} if sample is None else sample,
    }


def __lazy2format(sample_name: str, format_string: str, col2type: dict[str, polars.PolarsDataType]) -> polars.Expr:
    """Use sample name and column name to build sample column."""
    expression = []

    for col_name in format_string.split(":"):
        lazy_name = f"{sample_name}_{col_name.lower()}"
        if col_name == "GT":
            expression.append(
                polars.col(lazy_name)
                .cast(polars.Utf8)
                .fill_null("./.")
                .str.replace("1", "0/1")
                .str.replace("2", "1/1"),
            )
        elif isinstance(col2type[lazy_name], polars.List):
            expression.append(
                polars.col(lazy_name).cast(polars.List(polars.Utf8)).fill_null(["."]).list.join(","),
            )
        else:
            expression.append(
                polars.col(lazy_name).cast(polars.Utf8).fill_null("."),
            )

    return polars.concat_str(expression, separator=":")


def from_lazyframe(
    lf: polars.LazyFrame,
    output_path: pathlib.Path,
    renaming: RenameCol = DEFAULT_RENAME,
) -> None:
    """Write polars.LazyFrame in vcf format.

    Chromosome name mapping table:
      - 23: X
      - 24: Y
      - 25: MT

    All other chromosome number isn't change.

    Warning: This function perform LazyFrame.collect() before write csv, this can have a significant impact on memory usage.

    Args:
        lf: LazyFrame contains information.
        output_path: Path to where vcf to write.

    Returns:
        None
    """
    select_column: list[str] = []

    lf = lf.with_columns(
        [
            polars.col(renaming["#CHROM"])
            .cast(polars.Utf8)
            .str.replace("23", "X")
            .str.replace("24", "Y")
            .str.replace("25", "MT")
            .alias("#CHROM"),
            polars.col(renaming["POS"]).alias("POS"),
            polars.col(renaming["ID"]).alias("ID"),
            polars.col(renaming["REF"]).alias("REF"),
            polars.col(renaming["ALT"]).alias("ALT"),
        ],
    )

    select_column.extend(["#CHROM", "POS", "ID", "REF", "ALT"])

    header = __generate_header(lf, renaming["INFO"], list(renaming["sample"].keys()), renaming["FORMAT"])

    if renaming["QUAL"] != ".":
        lf = lf.with_columns([polars.col(renaming["QUAL"]).alias("QUAL")])
    else:
        lf = lf.with_columns([polars.lit(".").alias("QUAL")])

    select_column.append("QUAL")

    if renaming["FILTER"] != ".":
        lf = lf.with_columns([polars.col(renaming["FILTER"]).alias("FILTER")])
    else:
        lf = lf.with_columns([polars.lit(".").alias("FILTER")])

    select_column.append("FILTER")

    lf = add_info_column(lf, renaming["INFO"]) if renaming["INFO"] else lf.with_columns(polars.lit(".").alias("INFO"))

    select_column.append("INFO")

    if renaming["FORMAT"]:
        lf = lf.with_columns(polars.lit(renaming["FORMAT"]).alias("FORMAT"))
        select_column.append("FORMAT")

    if renaming["FORMAT"] and renaming["sample"]:
        for sample_name in renaming["sample"]:
            lf = lf.with_columns(
                [
                    __lazy2format(sample_name, renaming["FORMAT"], dict(zip(lf.columns, lf.dtypes))).alias(sample_name),
                ],
            )
            select_column.append(sample_name)

    lf = lf.select([polars.col(col) for col in select_column])

    with open(output_path, "wb") as fh:
        fh.write(header.encode())
        fh.write(lf.collect().write_csv(separator="\t").encode())


def add_info_column(lf: polars.LazyFrame, vcfinfo2parquet_name: list[tuple[str, str]]) -> polars.LazyFrame:
    """Construct an INFO column from multiple columns of lf.

    Args:
        lf: A dataframe.
        vcfinfo2parquet_name: List of vcf column name and lf column name.

    Returns:
        LazyFrame with INFO column and remove lf column use.
    """
    lf = lf.with_columns(
        [
            polars.col(name).list.join(",").fill_null(".").alias(name)
            for name, dtype in zip(lf.columns, lf.dtypes)
            if isinstance(dtype, polars.List)
        ],
    )

    lf = lf.with_columns(
        [
            polars.col(name).cast(str).fill_null(".").alias(name)
            for name, dtype in zip(lf.columns, lf.dtypes)
            if not isinstance(dtype, polars.List)
        ],
    )

    lf = lf.with_columns(
        [
            polars.concat_str(
                [
                    polars.concat_str(
                        [
                            polars.lit(vcf_name),
                            polars.lit("="),
                            polars.col(parquet_name),
                        ],
                    )
                    for vcf_name, parquet_name in vcfinfo2parquet_name
                ],
                separator=";",
            ).alias("INFO"),
        ],
    )

    lf = lf.drop([p for (v, p) in vcfinfo2parquet_name])

    return lf


def __generate_header(
    lf: polars.LazyFrame,
    vcfinfo2parquet_name: list[tuple[str, str]] | None = None,
    samples: list[str] | None = None,
    format_string: str | None = None,
) -> str:
    """Generate header of vcf file.

    Args:
        lf: Dataframe
        vcfinfo2parquet_name: vcf column name link to column name

    Returns:
        The header of vcf file
    """
    header = """##fileformat=VCFv4.3
##source=VariantPlanner
"""
    col2type = dict(zip(lf.columns, lf.dtypes))

    if vcfinfo2parquet_name:
        for vcf_name, col_name in vcfinfo2parquet_name:
            number = "." if isinstance(col2type[col_name], polars.List) else "1"
            type_ = "String"
            type_ = "Float" if isinstance(col2type[col_name], polars.Float64) else type_
            type_ = "Integer" if isinstance(col2type[col_name], polars.Int64) else type_

            header += f'##INFO=<ID={vcf_name},Number={number},Type={type_},Description="Unknow">\n'

    if samples and format_string:
        sample = samples[0]
        for col in format_string.split(":"):
            col_name = f"{sample}_{col.lower()}"

            number = "." if isinstance(col2type[col_name], polars.List) else "1"
            type_ = "String"
            type_ = "Float" if isinstance(col2type[col_name], polars.Float64) else type_
            type_ = "Integer" if isinstance(col2type[col_name], polars.Int64) else type_

            header += f'##FORMAT=<ID={col.upper()},Number={number},Type={type_},Description="Unknow">\n'

    return header

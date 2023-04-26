"""Function to manage input and output of vcf file."""

# std import
from __future__ import annotations

import logging
import re
import typing

if typing.TYPE_CHECKING:  # pragma: no cover
    import pathlib

# 3rd party import
import polars

# project import
from variantplanner import normalization
from variantplanner.exception import NotAVCFError

MINIMAL_COL_NUMBER: int = 8
SAMPLE_COL_BEGIN: int = 9

logger = logging.getLogger("io.vcf")


def info2expr(input_path: pathlib.Path, select_info: set[str] | None = None) -> list[polars.Expr]:
    """Read vcf header to generate a list of polars.Expr to extract info.

    Args:
        input_path: Path to vcf file.
        select_info: List of info field

    Returns:
        List of polars expr to parse info columns.

    Raises:
        NotAVCFError: If a line not starts with '#'
        NotAVCFError: If all line not start by '#CHR'
    """
    info_re = re.compile(
        r"ID=(?P<id>([A-Za-z_][0-9A-Za-z_.]*|1000G)),Number=(?P<number>[ARG0-9\.]+),Type=(?P<type>Integer|Float|String)",
    )

    expressions: list[polars.Expr] = []

    with open(input_path) as fh:
        for line in fh:
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
                else:
                    local_expr = local_expr.str.split(",")
                    if search["type"] == "Integer":
                        local_expr = local_expr.cast(polars.List(polars.Int64))
                    elif search["type"] == "Float":
                        local_expr = local_expr.cast(polars.List(polars.Float64))

                expressions.append(local_expr.alias(search["id"]))

    raise NotAVCFError(input_path)


def sample_index(input_path: pathlib.Path) -> dict[str, int] | None:
    """Read vcf header to generate an association map between sample name and index.

    Args:
        input_path: Path to vcf file.

    Returns:
        Map that associate a sample name to is sample index.

    Raises:
        NotAVCFError: If a line not starts with '#'
        NotAVCFError: If all line not start by '#CHR'
    """
    with open(input_path) as fh:
        for line in fh:
            if not line.startswith("#"):
                raise NotAVCFError(input_path)

            if line.startswith("#CHR"):
                split_line = line.strip().split("\t")
                if len(split_line) <= MINIMAL_COL_NUMBER:
                    return None

                return {sample: i for (i, sample) in enumerate(split_line[SAMPLE_COL_BEGIN:])}

    raise NotAVCFError(input_path)


def __column_name(input_path: pathlib.Path) -> list[str]:
    """Read vcf header to generate list of column name.

    Args:
        input_path: Path to vcf file.

    Returns:
        List of lazyframe columns name.

    Raises:
        NotAVCFError: If a line not starts with '#'
        NotAVCFError: If all line not start by '#CHR'
    """
    with open(input_path) as fh:
        for line in fh:
            if line.startswith("#CHR"):
                split_line = line.strip().split("\t")
                cols_name = ["chr", "pos", "vid", "ref", "alt", "qual", "filter", "info"]
                if len(split_line) > MINIMAL_COL_NUMBER:
                    cols_name.append("format")
                    for sample in split_line[SAMPLE_COL_BEGIN:]:
                        cols_name.append(sample)

                return cols_name
            if not line.startswith("#"):
                raise NotAVCFError(input_path)

    raise NotAVCFError(input_path)


def into_lazyframe(input_path: pathlib.Path) -> polars.LazyFrame:
    """Read a vcf file and convert it in lazyframe.

    Args:
        input_path: Path to vcf file.

    Returns:
        A lazyframe that containt vcf information ('chr', 'pos', 'vid', 'ref', 'alt', 'qual', 'filter', 'info', ['format'], ['genotypes',…], 'id').
    """
    col_name = {f"column_{i}": name for (i, name) in enumerate(__column_name(input_path), start=1)}

    lf = polars.scan_csv(
        input_path,
        separator="\t",
        comment_char="#",
        has_header=False,
        dtypes={"column_1": polars.Utf8},
        ignore_errors=True,
    )

    lf = lf.rename(col_name)

    lf = normalization.chromosome2integer(lf)

    lf = normalization.add_variant_id(lf)

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
            "INFO": dict[str, str],
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
    "INFO": {},
}


def build_rename_column(
    chromosome: str,
    pos: str,
    identifier: str,
    ref: str,
    alt: str,
    qual: str | None = ".",
    filter_col: str | None = ".",
    info: dict[str, str] | None = None,
) -> RenameCol:
    """Helper to generate rename column for from_lazyframe.

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
        "INFO": info if info is not None else {},
    }


def from_lazyframe(
    lf: polars.LazyFrame,
    output_path: pathlib.Path,
    renaming: RenameCol = DEFAULT_RENAME,
) -> None:
    """Write lazyframe in vcf format.

    Only variants information is write no info or genotype.

    Chromosome name mapping table:
        - 23: X
        - 24: Y
        - 25: MT

    All other number isn't change.

    Warning: This function perform LazyFrame.collect() before write csv, this can have a significant impact on memory usage

    Args:
        lf: LazyFrame contains information.
        output_path: Path to where vcf to write.

    Returns:
        None
    """
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

    header = generate_header(lf, renaming["INFO"])

    if renaming["QUAL"] != ".":
        lf = lf.with_columns([polars.col(renaming["QUAL"]).alias("QUAL")])
    else:
        lf = lf.with_columns([polars.lit(".").alias("QUAL")])

    if renaming["FILTER"] != ".":
        lf = lf.with_columns([polars.col(renaming["FILTER"]).alias("FILTER")])
    else:
        lf = lf.with_columns([polars.lit(".").alias("FILTER")])

    if "INFO" in renaming and renaming["INFO"]:
        lf = add_info_column(lf, renaming["INFO"])
    else:
        lf = lf.with_columns(
            polars.lit(".").alias("INFO"),
        )

    lf = lf.select([polars.col(colname) for colname in renaming])

    with open(output_path, "wb") as fh:
        fh.write(header.encode())
        fh.write(lf.collect().write_csv(separator="\t").encode())


def add_info_column(lf: polars.LazyFrame, vcfinfo2parquet_name: dict[str, str]) -> polars.LazyFrame:
    """Add INFO column in polars lazyframe.

    Args:
        lf: The dataframe
        vcfinfo2parquet_name: vcf column name link to column name

    Returns:
        LazyFrame with INFO column and remove select column
    """
    polars.Config().set_tbl_cols(100)

    lf = lf.with_columns(
        [
            polars.col(name).arr.join(",").fill_null(".").alias(name)
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
                    for vcf_name, parquet_name in vcfinfo2parquet_name.items()
                ],
                separator=";",
            ).alias("INFO"),
        ],
    )

    lf = lf.drop(list(vcfinfo2parquet_name.values()))

    return lf


def generate_header(lf: polars.LazyFrame, vcfinfo2parquet_name: dict[str, str] | None = None) -> str:
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

    if vcfinfo2parquet_name:
        col2type = dict(zip(lf.columns, lf.dtypes))
        for vcf_name, col_name in vcfinfo2parquet_name.items():
            number = "." if isinstance(col2type[col_name], polars.List) else "1"
            type_ = "String"
            type_ = "Float" if isinstance(col2type[col_name], polars.Float64) else type_
            type_ = "Integer" if isinstance(col2type[col_name], polars.Int64) else type_

            header += f'##INFO=<ID={vcf_name},Number={number},Type={type_},Description="Unknow">\n'

    return header

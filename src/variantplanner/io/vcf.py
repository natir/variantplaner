"""Function to manage input and output of vcf file."""

# std import
from __future__ import annotations

import re
import typing

if typing.TYPE_CHECKING:
    import pathlib

# 3rd party import
import polars

# project import
from variantplanner import normalization
from variantplanner.exception import NotAVCFError

MINIMAL_COL_NUMBER: int = 8
SAMPLE_COL_BEGIN: int = 9


def info2expr(input_path: pathlib.Path, select_info: set[str] | None) -> list[polars.Expr]:
    """Read vcf header to generate a list of polars.Expr to extract info.

    Args:
        input_path: Path to vcf file.
        select_info: List of info field

    Returns:
        List of polars expr to parse info columns.
    """
    info_re = re.compile(
        r"ID=(?P<id>([A-Za-z_][0-9A-Za-z_.]*|1000G)),Number=(?P<number>[ARG0-9\.]+),Type=(?P<type>Integer|Float|String)",
    )

    type2regex = {
        "Integer": r"((\d+),?)+",
        "Float": r"((\d+.\d+),?)+",
        "String": r"(([^,;]+),?)+",
    }

    expressions: list[polars.Expr] = []

    with open(input_path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return expressions
            if (search := info_re.search(line)) is not None and (select_info is None or search["id"] in select_info):
                regex = rf"{search['id']}="
                regex += type2regex[search["type"]]
                regex += r";?"

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

    raise NotAVCFError


def sample_index(input_path: pathlib.Path) -> dict[str, int] | None:
    """Read vcf header to generate an association map between sample name and index.

    Args:
        input_path: Path to vcf file.

    Returns:
        Map that associate a sample name to is sample index.
    """
    with open(input_path) as fh:
        for line in fh:
            if line.startswith("#CHR"):
                split_line = line.split("\t")
                if len(split_line) > MINIMAL_COL_NUMBER:
                    return None

                return {sample: i for (i, sample) in enumerate(split_line[SAMPLE_COL_BEGIN:])}
            if not line.startswith("#"):
                raise NotAVCFError

    raise NotAVCFError


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
                split_line = line.split("\t")
                cols_name = ["chr", "pos", "vid", "ref", "alt", "qual", "filter", "info"]
                if len(split_line) > MINIMAL_COL_NUMBER:
                    cols_name.append("format")
                    for sample in split_line[SAMPLE_COL_BEGIN:]:
                        cols_name.append(sample)

                return cols_name
            if not line.startswith("#"):
                raise NotAVCFError

    raise NotAVCFError


def into_lazyframe(input_path: pathlib.Path) -> polars.LazyFrame:
    """Read a vcf file and convert it in lazyframe.

    Args:
        input_path: Path to vcf file.

    Returns:
        A lazyframe that containt vcf information.
    """
    lf = polars.scan_csv(
        input_path,
        separator="\t",
        comment_char="#",
        has_header=False,
        dtypes={"column_1": polars.Utf8},
        ignore_errors=True,
    )

    lf = lf.rename({f"column_{i}": name for (i, name) in enumerate(__column_name(input_path))})

    lf = normalization.chromosome2integer(lf)

    lf = normalization.add_variant_id(lf)

    return lf


def from_lazyframe(lf: polars.LazyFrame, output_path: pathlib.Path) -> None:
    """Write lazyframe in vcf format.

    Only variants information is write no info or genotype.

    Chromosome name mapping table:
        - 23: X
        - 24: Y
        - 25: MT

    All other numbre isn't change.

    Warning: This function perform LazyFrame.collect() before write csv, this can have a significant impact on memory usage

    Args:
        lf: LazyFrame contains information.
        output_path: Path to where vcf to write.

    Returns:
        None
    """
    lf = lf.select(
        [
            polars.col("chr")
            .cast(polars.Utf8)
            .str.replace("23", "X")
            .str.replace("24", "Y")
            .str.replace("25", "MT")
            .alias("#CHROM"),
            polars.col("pos").alias("POS"),
            polars.col("id").alias("ID"),
            polars.col("ref").alias("REF"),
            polars.col("alt").alias("ALT"),
            polars.lit(".").alias("QUAL"),
            polars.lit(".").alias("FILTER"),
            polars.lit(".").alias("INFO"),
        ],
    )

    header = """##fileformat=VCFv4.3
    ##source=VariantPlanner
    """

    with open(output_path, "wb") as fh:
        fh.write(header.encode())

    lf.collect().write_csv(output_path, separator="\t")

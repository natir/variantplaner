"""Module contains command line entry point function."""

# Why does this file exist, and why not put this in `__main__`?
#
# You might be tempted to import things from `__main__` later,
# but that will cause problems: the code will get executed twice:
#
# - When you run `python -m variantplaner` python will execute
#   `__main__.py` as a script. That means there won't be any
#   `variantplaner.__main__` in `sys.modules`.
# - When you import `__main__` it will get executed again (as a module) because
#   there's no `variantplaner.__main__` in `sys.modules`.

# std import
from __future__ import annotations

import logging
import os
import pathlib
import sys

# 3rd party import
import click
import polars

# project import
from variantplaner import debug, exception, extract, generate, io, struct, cli

@cli.main.group("vcf2parquet", chain=True)
@click.pass_context
@click.option(
    "-i",
    "--input-path",
    help="Path to vcf input file",
    type=click.Path(exists=True, dir_okay=False, readable=True, allow_dash=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-c",
    "--chrom2length-path",
    help="CSV file that associates a chromosome name with its size",
    type=click.Path(dir_okay=False, writable=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-a",
    "--append",
    help="vcf2parquet switch in append mode",
    type=bool,
    is_flag=True,
)
def vcf2parquet(
        ctx: click.Context,
        input_path: pathlib.Path,
        chrom2length_path: pathlib.Path | None,
        append: bool,
) -> None:
    logger = logging.getLogger("vcf2parquet")

    logger.debug(f"parameter: {input_path=} {chrom2length_path=} {append=}")

    try:
        logger.debug(f"Start extract header")
        header = io.vcf.extract_header(input_path)
        logger.debug(f"End extract header")
    except exception.NotAVCFError:
        logger.exception("")
        sys.exit(11)

    # Read vcf and manage structural variant
    logger.debug(f"Start read vcf")
    lf = io.vcf.into_lazyframe(input_path, chrom2length_path, extension=io.vcf.IntoLazyFrameExtension.MANAGE_SV)
    logger.debug(f"End read vcf")

    ctx.obj = {
        "vcf_path": input_path,
        "lazyframe": lf,
        "append": append,
        "header": header,
    }


@vcf2parquet.command("variants")
@click.pass_context
@click.option(
    "-o",
    "--output-path",
    help="Path where variants will be written",
    type=click.Path(writable=True, path_type=pathlib.Path),
    required=True,
)
def variants(
        ctx: click.Context,
        output_path: pathlib.Path,
) -> None:
    """Write variants"""

    logger = logging.getLogger("vcf2parquet.variants")

    ctx.ensure_object(dict)

    lf = ctx.obj["lazyframe"]
    append = ctx.obj["append"]

    logger.debug(f"parameter: {output_path=}")

    logger.info(f"Start write variants in {output_path}")
    extract.variants(lf).sink_parquet(output_path, maintain_order=False)
    logger.info(f"End write variants in {output_path}")


@vcf2parquet.command("genotypes")
@click.pass_context
@click.option(
    "-o",
    "--output-path",
    help="Path where genotypes will be written",
    type=click.Path(writable=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-f",
    "--format-string",
    help="Value of FORMAT column, line not match with this are ignored",
    type=str,
    default="GT:AD:DP:GQ",
    show_default=True,
)
def genotypes(
        ctx: click.Context,
        output_path: pathlib.Path,
        format_string: str = "GT:AD:DP:GQ",
) -> None:
    """Write variants"""

    logger = logging.getLogger("vcf2parquet.variants")

    ctx.ensure_object(dict)

    lf = ctx.obj["lazyframe"]
    append = ctx.obj["append"]
    header = ctx.obj["header"]
    input_path = ctx.obj["vcf_path"]

    logger.debug(f"parameter: {output_path=} {format_string=}")

    try:
        genotypes_lf = extract.genotypes(lf, io.vcf.format2expr(header, input_path), format_string)
    except exception.NoGenotypeError:
        logger.exception("")
        sys.exit(12)

    logger.info(f"Start write genotypes in {output_path}")
    genotypes_lf.sink_parquet(output_path, maintain_order=False)
    logger.info(f"End write genotypes in {output_path}")


@vcf2parquet.command("annotations")
@click.pass_context
@click.option(
    "-o",
    "--output-path",
    help="Path where variants will be written",
    type=click.Path(writable=True, path_type=pathlib.Path),
    required=True,
)
@click.pass_context
@click.option(
    "-i",
    "--info",
    multiple=True,
    help="List of info fields that are kept if this list is empty all fields are kept only the first vcf file header is read",
    type=str,
)
@click.option(
    "-r",
    "--rename-id",
    help="Set column name of variant id",
    type=str,
)
def annotations(
        ctx: click.Context,
        output_path: pathlib.Path,
        info: set[str] | None = None,
        rename_id: str | None = None,
) -> None:
    """Write annotations"""

    logger = logging.getLogger("vcf2parquet.annotations")

    ctx.ensure_object(dict)

    lf = ctx.obj["lazyframe"]
    append = ctx.obj["append"]
    header = ctx.obj["header"]
    input_path = ctx.obj["vcf_path"]

    logger.debug(f"parameter: {output_path=}")

    logger.info(f"Start extract annotations")
    annotations_lf = lf.with_columns(io.vcf.info2expr(header, input_path, info))
    annotations_lf = annotations_lf.drop(["chr", "pos", "ref", "alt", "filter", "qual", "info"])
    if rename_id:
        logger.info(f"Rename vcf variant id in {rename_id}")
        lf = lf.rename({"vid": rename_id})
    logger.info(f"End extract annotations")

    logger.info(f"Start write annotations in {output_path}")
    annotations_lf.sink_parquet(output_path, maintain_order=False)
    logger.info(f"End write annotations in {output_path}")

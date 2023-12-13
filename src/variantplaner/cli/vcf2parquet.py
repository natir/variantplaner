"""Module contains vcf2parquet subcommand entry point function."""

# std import
from __future__ import annotations

import logging
import pathlib
import sys

# 3rd party import
import click

# project import
from variantplaner import cli, exception, extract, io


@cli.main.group("vcf2parquet", chain=True) # type: ignore[has-type]
@click.pass_context
@click.option(
    "-i",
    "--input-path",
    help="Path to vcf input file.",
    type=click.Path(
        exists=True,
        dir_okay=False,
        readable=True,
        allow_dash=True,
        path_type=pathlib.Path,
    ),
    required=True,
)
@click.option(
    "-c",
    "--chrom2length-path",
    help="CSV file that associates a chromosome name with its size.",
    type=click.Path(dir_okay=False, writable=True, path_type=pathlib.Path),
)
@click.option(
    "-a",
    "--append",
    help="Switch in append mode.",
    type=bool,
    is_flag=True,
)
def vcf2parquet(
    ctx: click.Context,
    input_path: pathlib.Path,
    chrom2length_path: pathlib.Path | None,
    *,
    append: bool,
) -> None:
    """Convert a vcf in parquet."""
    logger = logging.getLogger("vcf2parquet")

    logger.debug(f"parameter: {input_path=} {chrom2length_path=} {append=}")

    if not chrom2length_path:
        logger.error("--chrom2length-path argument is required")

    try:
        logger.debug("Start extract header")
        headers = io.vcf.extract_header(input_path)
        logger.debug("End extract header")
    except exception.NotAVCFError:
        logger.error("Input file seems no be a vcf")  # noqa: TRY400  we are in cli exception isn't readable
        sys.exit(11)

    # Read vcf and manage structural variant
    logger.debug("Start read vcf")
    lf = io.vcf.into_lazyframe(input_path, chrom2length_path, extension=io.vcf.IntoLazyFrameExtension.MANAGE_SV)
    logger.debug("End read vcf")

    ctx.obj["vcf_path"] = input_path
    ctx.obj["lazyframe"] = lf
    ctx.obj["append"] = append
    ctx.obj["headers"] = headers


@vcf2parquet.command("variants")
@click.pass_context
@click.option(
    "-o",
    "--output-path",
    help="Path where variants will be written.",
    type=click.Path(writable=True, path_type=pathlib.Path),
    required=True,
)
def variants(
    ctx: click.Context,
    output_path: pathlib.Path,
) -> None:
    """Write variants."""
    logger = logging.getLogger("vcf2parquet.variants")

    lf = ctx.obj["lazyframe"]
    append = ctx.obj["append"]  # noqa: F841 not used now

    logger.debug(f"parameter: {output_path=}")

    logger.info(f"Start write variants in {output_path}")
    extract.variants(lf).sink_parquet(output_path, maintain_order=False)
    logger.info(f"End write variants in {output_path}")


@vcf2parquet.command("genotypes")
@click.pass_context
@click.option(
    "-o",
    "--output-path",
    help="Path where genotypes will be written.",
    type=click.Path(writable=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-f",
    "--format-string",
    help="Value of FORMAT column, line not match with this are ignored.",
    type=str,
    default="GT:AD:DP:GQ",
    show_default=True,
)
def genotypes(
    ctx: click.Context,
    output_path: pathlib.Path,
    format_string: str = "GT:AD:DP:GQ",
) -> None:
    """Write genotypes."""
    logger = logging.getLogger("vcf2parquet.genotypes")

    lf = ctx.obj["lazyframe"]
    append = ctx.obj["append"]  # noqa: F841 not used now
    headers = ctx.obj["headers"]
    input_path = ctx.obj["vcf_path"]

    logger.debug(f"parameter: {output_path=} {format_string=}")

    try:
        genotypes_lf = extract.genotypes(lf, io.vcf.format2expr(headers, input_path), format_string)
    except exception.NoGenotypeError:
        logger.error("It's seems vcf not contains genotypes information.")  # noqa: TRY400  we are in cli exception isn't readable
        sys.exit(12)

    logger.info(f"Start write genotypes in {output_path}")
    genotypes_lf.sink_parquet(output_path, maintain_order=False)
    logger.info(f"End write genotypes in {output_path}")


@vcf2parquet.command("annotations")
@click.pass_context
@click.option(
    "-o",
    "--output-path",
    help="Path where genotypes will be written.",
    type=click.Path(writable=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-i",
    "--info",
    help="List of info fields that are kept if this list is empty all fields are kept only the first vcf file header is read.",
    cls=cli.MultipleValueOption,
    type=str,
)
@click.option(
    "-r",
    "--rename-id",
    help="Set column name of variant id.",
    type=str,
)
def annotations_subcommand(
    ctx: click.Context,
    output_path: pathlib.Path,
    info: set[str] | None = None,
    rename_id: str | None = None,
) -> None:
    """Write annotations."""
    logger = logging.getLogger("vcf2parquet.annotations")

    lf = ctx.obj["lazyframe"]
    append = ctx.obj["append"]  # noqa: F841 not used now
    headers = ctx.obj["headers"]
    input_path = ctx.obj["vcf_path"]

    logger.debug(f"parameter: {output_path=}")

    logger.info("Start extract annotations")
    annotations_lf = lf.with_columns(io.vcf.info2expr(headers, input_path, info))
    annotations_lf = annotations_lf.drop(["chr", "pos", "ref", "alt", "filter", "qual", "info"])

    if rename_id:
        logger.info(f"Rename vcf variant id in {rename_id}")
        annotations_lf = annotations_lf.rename({"vid": rename_id})
    logger.info("End extract annotations")

    logger.info(f"Start write annotations in {output_path}")
    annotations_lf.sink_parquet(output_path, maintain_order=False)
    logger.info(f"End write annotations in {output_path}")


@vcf2parquet.command("headers")
@click.pass_context
@click.option(
    "-o",
    "--output-path",
    help="Path where header will be written.",
    type=click.Path(writable=True, path_type=pathlib.Path),
    required=True,
)
def headers(
    ctx: click.Context,
    output_path: pathlib.Path,
) -> None:
    """Write vcf headers."""
    logger = logging.getLogger("vcf2parquet.headers")

    headers = ctx.obj["headers"]

    logger.debug(f"parameter: {output_path=}")

    logger.info(f"Start write headers in {output_path}")
    with open(output_path, "w") as fh_out:
        for line in headers:
            print(line, file=fh_out)
    logger.info(f"End write headers in {output_path}")

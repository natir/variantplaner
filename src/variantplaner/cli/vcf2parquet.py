"""Module contains vcf2parquet subcommand entry point function."""

# std import
from __future__ import annotations

import logging
import pathlib
import sys

# 3rd party import
import click
import polars

# project import
from variantplaner import Vcf, VcfParsingBehavior, cli, exception

logger = logging.getLogger("__name__")


@cli.main.group("vcf2parquet", chain=True)  # type: ignore[has-type]
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
@click.option(
    "-s",
    "--keep-star",
    help="Keep variant with * in alt.",
    type=bool,
    is_flag=False,
)
def vcf2parquet(
    ctx: click.Context,
    input_path: pathlib.Path,
    chrom2length_path: pathlib.Path | None,
    *,
    append: bool,
    keep_star: bool,
) -> None:
    """Convert a vcf in parquet."""
    logger = logging.getLogger("vcf2parquet")

    logger.debug(f"parameter: {input_path=} {chrom2length_path=} {append=} {keep_star=}")

    lf = Vcf()

    # Read vcf and manage structural variant
    logger.debug("Start read vcf")
    try:
        beahvior = VcfParsingBehavior.MANAGE_SV
        if keep_star:
            beahvior |= VcfParsingBehavior.KEEP_STAR

        lf.from_path(input_path, chrom2length_path, behavior=beahvior)
    except exception.NotVcfHeaderError:
        logging.error(f"Path {input_path} seems not contains Vcf.")  # noqa: TRY400  we are in cli exception isn't readable
        sys.exit(11)
    except exception.NotAVCFError:
        logging.error(f"Path {input_path} seems not contains Vcf.")  # noqa: TRY400  we are in cli exception isn't readable
        sys.exit(12)
    except exception.NoContigsLengthInformationError:
        logging.exception("Vcf didn't contains contigs length information you could use chrom2length-path argument.")
        sys.exit(13)
    logger.debug("End read vcf")

    ctx.obj["vcf_path"] = input_path
    ctx.obj["lazyframe"] = lf
    ctx.obj["append"] = append
    ctx.obj["headers"] = lf.header


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
    append = ctx.obj["append"]

    logger.debug(f"parameter: {output_path=}")

    logger.info(f"Start write variants in {output_path}")
    variants = lf.variants()

    if append:
        variants = __append(output_path, variants)

    try:
        variants.sink_parquet(output_path, maintain_order=False)
    except polars.exceptions.InvalidOperationError:
        variants.collect(streaming=True).write_parquet(output_path)
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
    append = ctx.obj["append"]

    logger.debug(f"parameter: {output_path=} {format_string=}")

    try:
        genotypes_data = lf.genotypes(format_string)
    except exception.NoGenotypeError:
        logger.error("It's seems vcf not contains genotypes information.")  # noqa: TRY400  we are in cli exception isn't readable
        sys.exit(12)

    if append:
        genotypes_data = __append(output_path, genotypes_data)

    logger.info(f"Start write genotypes in {output_path}")
    try:
        genotypes_data.lf.sink_parquet(output_path, maintain_order=False)
    except polars.exceptions.InvalidOperationError:
        genotypes_data.lf.collect(streaming=True).write_parquet(output_path)
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
    append = ctx.obj["append"]
    headers_obj = ctx.obj["headers"]

    logger.debug(f"parameter: {output_path=}")

    logger.info("Start extract annotations")
    annotations_data = lf.lf.with_columns(headers_obj.info_parser(info))
    annotations_data = annotations_data.drop(["chr", "pos", "ref", "alt", "filter", "qual", "info"])

    if rename_id:
        logger.info(f"Rename vcf variant id in {rename_id}")
        annotations_data = annotations_data.rename({"vid": rename_id})
    logger.info("End extract annotations")

    if append:
        annotations_data = __append(output_path, annotations_data)

    logger.info(f"Start write annotations in {output_path}")
    try:
        annotations_data.sink_parquet(output_path, maintain_order=False)
    except polars.exceptions.InvalidOperationError:
        annotations_data.collect(streaming=True).write_parquet(output_path)
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

    headers_obj = ctx.obj["headers"]

    logger.debug(f"parameter: {output_path=}")

    logger.info(f"Start write headers in {output_path}")
    with open(output_path, "w") as fh_out:
        for line in headers_obj:
            print(line, file=fh_out)
    logger.info(f"End write headers in {output_path}")


def __append(output_path: pathlib.Path, lf: polars.LazyFrame) -> polars.LazyFrame:
    """Concatenate contante of output_path with lf content if output_path exists.

    If parquet schema not match an error will be raise.
    """
    if not output_path.is_file():
        logger.error("Target file not exist append mode isn't apply.")
    else:
        lf = polars.concat([polars.scan_parquet(output_path), lf])

    return lf

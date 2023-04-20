"""Module that contains the command line application."""

# Why does this file exist, and why not put this in `__main__`?
#
# You might be tempted to import things from `__main__` later,
# but that will cause problems: the code will get executed twice:
#
# - When you run `python -m variantplanner` python will execute
#   `__main__.py` as a script. That means there won't be any
#   `variantplanner.__main__` in `sys.modules`.
# - When you import `__main__` it will get executed again (as a module) because
#   there's no `variantplanner.__main__` in `sys.modules`.

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
from variantplanner import exception, io, manipulation, structuration


@click.group(name="variantplanner")
@click.option("-t", "--threads", help="Number of threads usable", default=1, type=click.IntRange(0), show_default=True)
@click.option("-v", "--verbose", help="Verbosity level", count=True, type=click.IntRange(0, 4))
def main(threads: int, verbose: int) -> None:
    """Run VariantPlanner."""
    logging.basicConfig(
        style="{",
        format="{asctime} - {name}:{levelname}: {message}",
        encoding="utf-8",
        level=(4 - verbose) * 10,  # Python choose a bad logging levels order
        stream=sys.stderr,
    )

    logger = logging.getLogger("main")

    logger.debug(f"parameter: {threads=} {verbose=}")

    os.environ["POLARS_MAX_THREADS"] = str(threads)


###############
# vcf2parquet #
###############
@main.command()
@click.option(
    "-i",
    "--input-path",
    help="Path to vcf input file",
    type=click.Path(exists=True, dir_okay=False, readable=True, allow_dash=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-v",
    "--variants",
    help="Path where the variants will be written in parquet",
    type=click.Path(dir_okay=False, writable=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-g",
    "--genotypes",
    help="Path where the genotypes will be written in parquet",
    type=click.Path(dir_okay=False, writable=True, path_type=pathlib.Path),
)
def vcf2parquet(input_path: pathlib.Path, variants: pathlib.Path, genotypes: pathlib.Path) -> None:
    """Convert a vcf in parquet."""
    logger = logging.getLogger("vcf2parquet")

    logger.debug(f"parameter: {input_path=} {variants=} {genotypes=}")

    try:
        lf = io.vcf.into_lazyframe(input_path)
    except exception.NotAVCFError:
        logger.exception("")
        sys.exit(1)

    manipulation.extract_variants(lf).sink_parquet(variants)
    logger.info(f"finish write {variants}")

    if genotypes:
        try:
            manipulation.extract_genotypes(lf).sink_parquet(genotypes)
        except exception.NoGenotypeError:
            logger.exception("")
            sys.exit(2)

    logger.info(f"finish write {genotypes}")


#########
# Merge #
#########
@main.command()
@click.option(
    "-i",
    "--inputs-path",
    help="Paths of the variant files to be merged",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-m",
    "--merge-path",
    help="Path where merged variants will be written",
    type=click.Path(writable=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-b",
    "--bytes-memory-limit",
    help="Number of bytes used to build chunk of merge file",
    type=click.IntRange(0),
    show_default=True,
)
def merge(inputs_path: list[pathlib.Path], merge_path: pathlib.Path, bytes_memory_limit: int = 10_000_000_000) -> None:
    """Merge multiple variants parquet file in one.

    If you set TMPDIR, TEMP or TMP environment variable you can control where temp file is create.
    """
    logger = logging.getLogger("merge")

    logger.debug(f"parameter: {inputs_path=} {merge_path=} {bytes_memory_limit}")

    structuration.merge_variant.parquet(inputs_path, merge_path, bytes_memory_limit)


###############
# Annotations #
###############
@main.group("annotations")
@click.pass_context
@click.option(
    "-i",
    "--input-paths",
    help="Path to input file",
    type=click.Path(exists=True, dir_okay=False, readable=True, allow_dash=True, path_type=pathlib.Path),
    required=True,
    multiple=True,
)
@click.option(
    "-o",
    "--output-path",
    help="Path where variants annotations will be written in parquet",
    type=click.Path(writable=True, path_type=pathlib.Path),
    required=True,
)
def annotations_main(ctx: click.Context, input_paths: list[pathlib.Path], output_path: pathlib.Path) -> None:
    """Convert an annotation variation file in parquet."""
    logger = logging.getLogger("annotations")

    ctx.obj = {"input_paths": input_paths, "output_path": output_path}

    logger.debug(f"parameter: {input_paths=} {output_path=}")


@annotations_main.command("vcf")
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
def annotations_vcf(ctx: click.Context, info: set[str] | None = None, rename_id: str | None = None) -> None:
    """Convert an annotated vcf in parquet."""
    logger = logging.getLogger("annotations-vcf")

    ctx.ensure_object(dict)

    input_paths = ctx.obj["input_paths"]
    output_path = ctx.obj["output_path"]

    logger.debug(f"parameter: {input_paths=} {output_path=} {info=}")

    try:
        info_parser = io.vcf.info2expr(input_paths[0], info)
        lf = polars.concat(io.vcf.into_lazyframe(path) for path in input_paths)
    except exception.NotAVCFError:
        logger.exception("")
        sys.exit(1)

    lf = lf.with_columns(info_parser).drop(["chr", "pos", "ref", "alt", "filter", "qual", "info"])

    if rename_id:
        logger.info(f"Rename vcf variant id in {rename_id}")
        lf = lf.rename({"vid": rename_id})

    lf.collect().write_parquet(output_path)


@annotations_main.command("csv")
@click.pass_context
@click.option(
    "-c",
    "--chromosome",
    help="Name of chromosome column",
    type=str,
    required=True,
)
@click.option(
    "-p",
    "--position",
    help="Name of position column",
    type=str,
    required=True,
)
@click.option(
    "-r",
    "--reference",
    help="Name of reference column",
    type=str,
    required=True,
)
@click.option(
    "-a",
    "--alternative",
    help="Name of alternative column",
    type=str,
    required=True,
)
@click.option(
    "-i",
    "--info",
    multiple=True,
    help="List of columns that are kept if this list is empty all columns are kept",
    type=str,
)
def annotations_csv(
    ctx: click.Context,
    chromosome: str,
    position: str,
    reference: str,
    alternative: str,
    info: list[str],
) -> None:
    """Convert an annotated csv in parquet."""
    logger = logging.getLogger("annotations-vcf")

    ctx.ensure_object(dict)

    input_path = ctx.obj["input_path"]
    output_path = ctx.obj["output_path"]

    logger.debug(
        f"parameter: {input_path=} {output_path=} {chromosome=} {position=} {reference=} {alternative=} {info=}",
    )

    lf = io.csv.into_lazyframe(
        input_path,
        chromosome,
        position,
        reference,
        alternative,
        info,
        separator=",",
    )

    lf.sink_parquet(output_path)

    logging.info("end of annontations csv")


############
# Metadata #
############
@main.group()
@click.pass_context
@click.option(
    "-i",
    "--input-path",
    help="Path to input file",
    type=click.Path(exists=True, dir_okay=False, readable=True, allow_dash=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-o",
    "--output-path",
    help="Path where variants annotations will be written in parquet",
    type=click.Path(writable=True, path_type=pathlib.Path),
    required=True,
)
def metadata(ctx: click.Context, input_path: pathlib.Path, output_path: pathlib.Path) -> None:
    """Convert an metadata file in parquet."""
    logger = logging.getLogger("metadata")

    ctx.obj = {"input_path": input_path, "output_path": output_path}

    logger.debug(f"parameter: {input_path=} {output_path=}")


@metadata.command("json")
@click.pass_context
@click.option(
    "-f",
    "--fields",
    multiple=True,
    help="List of fields that are kept if this list is empty all fields are kept",
    type=str,
)
def metadata_json(ctx: click.Context, fields: list[str]) -> None:
    """Convert an metadata json in parquet."""
    logger = logging.getLogger("metadata-json")

    ctx.ensure_object(dict)

    input_path = ctx.obj["input_path"]
    output_path = ctx.obj["output_path"]

    logger.debug(f"parameter: {input_path=} {output_path=} {fields=}")


@metadata.command("csv")
@click.pass_context
@click.option(
    "-c",
    "--columns",
    multiple=True,
    help="List of columns that are kept if this list is empty all columns are kept",
    type=str,
)
def metadata_csv(ctx: click.Context, columns: list[str]) -> None:
    """Convert an metadata csv in parquet."""
    logger = logging.getLogger("metadata-csv")

    ctx.ensure_object(dict)

    input_path = ctx.obj["input_path"]
    output_path = ctx.obj["output_path"]

    logger.debug(f"parameter: {input_path=} {output_path=} {columns=}")

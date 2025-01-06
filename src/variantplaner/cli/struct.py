"""Module contains struct subcommand entry point function."""

# std import
from __future__ import annotations

import logging
import math
import os
import pathlib

# 3rd party import
import click

# project import
from variantplaner import cli
from variantplaner import struct as vp_struct


@cli.main.group("struct", chain=True)  # type: ignore[has-type]
@click.pass_context
@click.option(
    "-i",
    "--input-paths",
    help="Paths of the variant files to be merged.",
    cls=cli.MultipleValueOption,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-a",
    "--append",
    help="Switch in append mode.",
    type=bool,
    is_flag=True,
)
def struct(ctx: click.Context, input_paths: list[pathlib.Path], *, append: bool) -> None:
    """Subcommand to made struct operation on parquet file."""
    logger = logging.getLogger("struct")

    if not (isinstance(input_paths, (list, tuple))):
        input_paths = [input_paths]

    ctx.obj["input_paths"] = input_paths
    ctx.obj["append"] = append

    logger.debug(f"parameter: {input_paths=} {append=}")


@struct.command("variants")
@click.pass_context
@click.option(
    "-o",
    "--output-prefix",
    help="Prefix added before file where unique variants will be written.",
    type=click.Path(writable=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-c",
    "--chunk-size",
    help="Sizes of files chunks that will be grouped together to build the list of unique variants",
    type=click.IntRange(min=1),
    default=1_000_000,
    show_default=True,
)
@click.option(
    "-p",
    "--polars-threads",
    help="Number of threads use to merge each chunk of files.",
    type=click.IntRange(min=1),
    default=4,
    show_default=True,
)
def variants(
    ctx: click.Context,
    output_prefix: pathlib.Path,
    chunk_size: int,
    polars_threads: int,
) -> None:
    """Merge multiple variants parquet file in one.

    If you set TMPDIR, TEMP or TMP environment variable you can control where temp file is created.
    """
    logger = logging.getLogger("struct.variants")

    ctx.ensure_object(dict)

    input_paths = ctx.obj["input_paths"]
    append = ctx.obj["append"]

    logger.debug(f"parameter: {output_prefix=} {chunk_size=}")

    vp_struct.variants.merge(input_paths, output_prefix, chunk_size, polars_threads, append=append)


@struct.command("genotypes")
@click.pass_context
@click.option(
    "-p",
    "--prefix-path",
    help="Prefix add before genotype partitions",
    type=click.Path(file_okay=False, dir_okay=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-m",
    "--partition-mode",
    help="Partition mode",
    type=click.Choice(["random", "position"]),
    default="position",
    show_default=True,
)
@click.option(
    "-n",
    "--number-of-part",
    help="Maximal number of part, automaticly correct to next power of two",
    type=click.IntRange(min=1),
    default=256,
    show_default=True,
)
@click.option(
    "-f",
    "--file-per-thread",
    help="Number of file manage by on threads, reduce value to reduce memory usage",
    type=click.IntRange(min=1),
    default=15,
    show_default=True,
)
@click.option(
    "-P",
    "--polars-threads",
    help="Number of threads use by polars task",
    type=click.IntRange(1),
    default=4,
    show_default=True,
)
def genotypes(
    ctx: click.Context,
    prefix_path: pathlib.Path,
    partition_mode: str,
    number_of_part: int,
    file_per_thread: int,
    polars_threads: int,
) -> None:
    """Convert set of genotype parquet in hive like files structures."""
    logger = logging.getLogger("struct.genotypes")

    ctx.ensure_object(dict)

    input_paths = ctx.obj["input_paths"]
    threads = ctx.obj["threads"]
    append = ctx.obj["append"]

    os.environ["POLARS_MAX_THREADS"] = str(polars_threads)

    logger.debug(f"parameter: {prefix_path=} {partition_mode=} {number_of_part=} {file_per_thread=} {polars_threads=}")

    number_of_bits = math.ceil(math.log2(number_of_part))

    vp_struct.genotypes.hive(
        input_paths,
        prefix_path,
        threads,
        file_per_thread,
        append=append,
        number_of_bits=number_of_bits,
    )

"""Module contains transmission subcommand entry point function."""

# std import
from __future__ import annotations

import logging
import pathlib
import sys

# 3rd party import
import click
import polars

# project import
from variantplaner import cli, generate, io


@cli.main.command("transmission")  # type: ignore[has-type]
@click.option(
    "-g",
    "--genotypes-path",
    help="Path to genotypes parquet.",
    type=click.Path(exists=True, readable=True, allow_dash=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-p",
    "--pedigree-path",
    help="Path to pedigree file.",
    type=click.Path(exists=True, readable=True, path_type=pathlib.Path),
)
@click.option(
    "-i",
    "--index",
    help="Sample name of index.",
    type=str,
)
@click.option(
    "-m",
    "--mother",
    help="Sample name of mother.",
    type=str,
)
@click.option(
    "-f",
    "--father",
    help="Sample name of father.",
    type=str,
)
@click.option(
    "-o",
    "--output-path",
    help="Path where transmission will be write.",
    type=click.Path(writable=True, path_type=pathlib.Path),
)
def transmission(
    genotypes_path: pathlib.Path,
    output_path: pathlib.Path,
    pedigree_path: pathlib.Path | None,
    index: str | None,
    mother: str | None,
    father: str | None,
) -> None:
    """Generate transmission of a genotype set."""
    logger = logging.getLogger("vcf2parquet.genotypes")

    logger.debug(f"parameter: {genotypes_path=} {output_path=} {pedigree_path=} {index=} {mother=} {father=}")

    genotypes_lf = polars.scan_parquet(genotypes_path)

    if pedigree_path:
        pedigree_lf = io.ped.into_lazyframe(pedigree_path)
        transmission_lf = generate.transmission_ped(genotypes_lf, pedigree_lf)
    elif index:
        transmission_lf = generate.transmission(genotypes_lf, index, mother, father)
    else:
        logging.error("You must specify ped file or index almost sample name")
        sys.exit(31)

    transmission_lf.write_parquet(output_path)

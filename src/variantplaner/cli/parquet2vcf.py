"""Module contains parquet2vcf subcommand entry point function."""

# std import
from __future__ import annotations

import logging
import pathlib

# 3rd party import
import click
import polars

# project import
from variantplaner import Genotypes, Variants, Vcf, cli, io


@cli.main.command("parquet2vcf")  # type: ignore[has-type]
@click.option(
    "-v",
    "--variants-path",
    help="Path to variant parquet.",
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
    "-o",
    "--output-path",
    help="Path where the vcf is written.",
    type=click.Path(dir_okay=False, writable=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-g",
    "--genotypes-path",
    help="Path to genotype parquet.",
    type=click.Path(
        exists=True,
        dir_okay=False,
        readable=True,
        allow_dash=True,
        path_type=pathlib.Path,
    ),
)
@click.option(
    "-H",
    "--headers-path",
    help="Path to vcf header.",
    type=click.Path(
        exists=True,
        dir_okay=False,
        readable=True,
        allow_dash=True,
        path_type=pathlib.Path,
    ),
)
@click.option(
    "-s",
    "--select-chromosome",
    help="Output only record where chromosome match target",
    type=str,
)
@click.option(
    "-c",
    "--chromosome",
    help="Name of chromosome column.",
    type=str,
    default="chr",
    show_default=True,
)
@click.option(
    "-p",
    "--position",
    help="Name of position column.",
    type=str,
    default="pos",
    show_default=True,
)
@click.option(
    "-I",
    "--identifier",
    help="Name of identity column.",
    type=str,
    default="id",
    show_default=True,
)
@click.option(
    "-r",
    "--reference",
    help="Name of reference column.",
    default="ref",
    show_default=True,
)
@click.option(
    "-a",
    "--alternative",
    help="Name of alternative column.",
    default="alt",
    show_default=True,
)
@click.option(
    "-q",
    "--quality",
    help="Name of quality column.",
    type=str,
)
@click.option(
    "-f",
    "--filter",
    "filter_col",
    help="Name of filter column.",
    type=str,
)
@click.option(
    "-F",
    "--format",
    "format_str",
    help="Value of format column.",
    type=str,
)
def parquet2vcf(
    variants_path: pathlib.Path,
    output_path: pathlib.Path,
    genotypes_path: pathlib.Path | None = None,
    headers_path: pathlib.Path | None = None,
    select_chromosome: str | None = None,
    chromosome: str = "chr",
    position: str = "pos",
    identifier: str = "id",
    reference: str = "ref",
    alternative: str = "alt",
    quality: str | None = None,
    filter_col: str | None = None,
    format_str: str | None = None,
) -> None:
    """Convert variant parquet in vcf."""
    logger = logging.getLogger("vcf2parquet")

    logger.debug(
        f"parameter: {variants_path} {output_path} {genotypes_path} {headers_path} {chromosome} {position} {identifier} {reference} {alternative} {quality} {filter_col} {format_str}"
    )

    vcf = Vcf()

    lf = polars.scan_parquet(variants_path)
    if select_chromosome is not None:
        lf = lf.filter(polars.col(chromosome) == select_chromosome)

    vcf.set_variants(Variants(lf))

    if headers_path:
        vcf.header.from_files(headers_path)
    else:
        headers = None

    if genotypes_path and format_str:
        genotypes = Genotypes(polars.scan_parquet(genotypes_path))
        vcf.add_genotypes(genotypes)

        sample2vcf_col2polars_col: dict[str, dict[str, str]] = {}
        for sample in genotypes.samples_names():
            sample2vcf_col2polars_col[sample] = {}
            for format_col in format_str.split(":"):
                sample2vcf_col2polars_col[sample][format_col] = f"{sample}_{format_col.lower()}"
        rename_column = io.vcf.build_rename_column(
            chromosome,
            position,
            identifier,
            reference,
            alternative,
            quality,
            filter_col,
            [],
            format_str,
            sample2vcf_col2polars_col,
        )
    else:
        rename_column = io.vcf.build_rename_column(
            chromosome,
            position,
            identifier,
            reference,
            alternative,
            quality,
            filter_col,
        )

    io.vcf.lazyframe_in_vcf(
        vcf.lf,
        output_path,
        vcf_header=headers,
        renaming=rename_column,
    )

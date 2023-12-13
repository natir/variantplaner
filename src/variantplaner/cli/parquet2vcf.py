"""Module contains parquet2vcf subcommand entry point function."""

# std import
from __future__ import annotations

import logging
import pathlib

# 3rd party import
import click
import polars

# project import
from variantplaner import cli, extract, io


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

    variants_lf = polars.scan_parquet(variants_path)

    if genotypes_path and format_str:
        genotypes_lf = polars.scan_parquet(genotypes_path)
        sample_name = genotypes_lf.select("sample").collect().get_column("sample").to_list()
        merge_lf = extract.merge_variants_genotypes(variants_lf, genotypes_lf, sample_name)
        sample2vcf_col2polars_col: dict[str, dict[str, str]] = {}
        for sample in sample_name:
            sample2vcf_col2polars_col[sample] = {}
            for format_col in format_str.split(":"):
                sample2vcf_col2polars_col[sample][format_col] = f"{sample}_{format_col.lower()}"

        io.vcf.from_lazyframe(
            merge_lf,
            output_path,
            io.vcf.build_rename_column(
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
            ),
        )
    else:
        io.vcf.from_lazyframe(
            variants_lf,
            output_path,
            io.vcf.build_rename_column(
                chromosome,
                position,
                identifier,
                reference,
                alternative,
                quality,
                filter_col,
            ),
        )

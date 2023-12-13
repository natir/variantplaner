"""Module contains command line entry point function."""

# std import
from __future__ import annotations

import logging
import pathlib

# 3rd party import
import click
import polars

# project import
from variantplaner import cli


@cli.main.command(name="metadata")  # type: ignore[has-type]
@click.option(
    "-i",
    "--input-path",
    help="Input path.",
    type=click.Path(exists=True, readable=True, allow_dash=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-o",
    "--output-path",
    help="Output path.",
    type=click.Path(writable=True, allow_dash=True, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "-t",
    "--input-type",
    help="Type of input file.",
    type=click.Choice(["csv", "tsv", "ljson", "json"]),
    required=True,
)
def metadata(
    input_path: pathlib.Path,
    output_path: pathlib.Path,
    input_type: str,
) -> None:
    """Convert metadata file in parquet file."""
    logger = logging.getLogger("metadata")

    logger.debug(f"parameter: {input_path=} {output_path=} {input_type=}")

    if input_type == "csv":
        lf = polars.scan_csv(input_path)
        lf.sink_parquet(output_path, maintain_order=False)
    elif input_type == "tsv":
        lf = polars.scan_csv(input_path, separator="\t")
        lf.sink_parquet(output_path, maintain_order=False)
    elif input_type == "ljson":
        df = polars.read_ndjson(input_path)
        df.write_parquet(output_path)
    else:
        df = polars.read_json(input_path)

        df.write_parquet(output_path)

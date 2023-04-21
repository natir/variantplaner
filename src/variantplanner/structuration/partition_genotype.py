"""Function relate to partition of genotype."""

# std import
from __future__ import annotations

import logging
import typing

if typing.TYPE_CHECKING:  # pragma: no cover
    import pathlib

# 3rd party import
import polars

# project import

logger = logging.getLogger("structuration.partition_genotype")


def manage_group(
    prefix: pathlib.Path,
    group_columns: list[str],
) -> typing.Callable[[polars.DataFrame], polars.DataFrame]:
    """Function to generate function apply to group."""

    def get_value(column: str, columns: list[str], row: tuple[typing.Any, ...]) -> typing.Any:
        return row[columns.index(column)]

    def inner(group: polars.DataFrame) -> polars.DataFrame:
        row = group.row(0)
        columns = group.columns
        path = (
            prefix.joinpath(*[f"{column}={get_value(column, columns, row)}" for column in group_columns]) / "0.parquet"
        )

        write_or_add(group, path)

        return polars.DataFrame()

    return inner


def write_or_add(new_lf: polars.DataFrame, partition_path: pathlib.Path) -> None:
    """Create or add new data in parquet partition."""
    if partition_path.exists():
        old_lf = polars.read_parquet(partition_path)
        polars.concat([old_lf, new_lf]).write_parquet(partition_path)
    else:
        partition_path.parent.mkdir(parents=True, exist_ok=True)
        new_lf.write_parquet(partition_path)


def parquet(paths: list[pathlib.Path], output_prefix: pathlib.Path) -> None:
    """Reorganise genotypes struct in hive like struct.

    Args:
        paths: List of file you want reorganise
        output_prefix: prefix of hive

    Returns:
        None
    """
    # Iterate over each genotypes
    for path in paths:
        polars.scan_parquet(path).with_columns([polars.col("id").mod(50).alias("id_mod")]).groupby(
            "gt",
            "sample",
            "id_mod",
        ).apply(
            manage_group(
                output_prefix,
                [
                    "gt",
                    "sample",
                    "id_mod",
                ],
            ),
            schema=None,
        ).collect()

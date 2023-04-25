"""Function relate to partition of genotype."""

# std import
from __future__ import annotations

import logging
import multiprocessing
import typing

if typing.TYPE_CHECKING:  # pragma: no cover
    import pathlib

# 3rd party import
import polars

# project import

logger = logging.getLogger("struct.genotypes")


def manage_group(
    prefix: pathlib.Path,
    group_columns: list[str],
    basename: str,
) -> typing.Callable[[polars.DataFrame], polars.DataFrame]:
    """Function to generate function apply to group.

    Args:
        prefix: Prefix of hive
        group_columns: List of columns use to group
        basename: Basename of final parquet

    Returns:
        Function that perform operation on polars group by
    """

    def get_value(column: str, columns: list[str], row: tuple[typing.Any, ...]) -> typing.Any:
        return row[columns.index(column)]

    def inner(group: polars.DataFrame) -> polars.DataFrame:
        row = group.row(0)
        columns = group.columns
        path = (
            prefix.joinpath(*[f"{column}={get_value(column, columns, row)}" for column in group_columns])
            / f"{basename}.parquet"
        )

        write_or_add(group, path)

        return polars.DataFrame()

    return inner


def write_or_add(new_lf: polars.DataFrame, partition_path: pathlib.Path) -> None:
    """Create or add new data in parquet partition.

    Args:
        new_lf: Dataframe to add or write
        partition_path: Path where dataframe is write

    Returns:
        None
    """
    if partition_path.exists():
        old_lf = polars.read_parquet(partition_path)
        polars.concat([old_lf, new_lf]).write_parquet(partition_path)
    else:
        partition_path.parent.mkdir(parents=True, exist_ok=True)
        new_lf.write_parquet(partition_path)


def __hive_worker(path: pathlib.Path, output_prefix: pathlib.Path, basename: int) -> pathlib.Path:
    """Subprocess of hive function run in parallel.

    Args:
        path: List of file you want reorganise
        output_prefix: prefix of hive
        basename: name of parquet file

    Returns:
        None
    """
    polars.scan_parquet(path).with_columns(
        [
            polars.col("id").mod(50).alias("id_mod"),
            polars.col("sample").hash().mod(50).alias("sample_mod"),
        ],
    ).groupby(
        "id_mod",
        "sample_mod",
    ).apply(
        manage_group(
            output_prefix,
            [
                "id_mod",
                "sample_mod",
            ],
            str(basename),
        ),
        schema=None,
    ).collect()

    return path


def hive(paths: list[pathlib.Path], output_prefix: pathlib.Path, threads: int = 1) -> None:
    """Reorganise genotypes struct in hive like struct.

    Args:
        paths: List of file you want reorganise
        output_prefix: prefix of hive

    Returns:
        None
    """
    threads = min(threads, len(paths))

    with multiprocessing.get_context("spawn").Pool() as pool:
        pool.starmap(__hive_worker, [(path, output_prefix, i) for i, path in enumerate(paths)])

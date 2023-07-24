"""Function relate to genotype structuration."""

# std import
from __future__ import annotations

import itertools
import logging
import multiprocessing
import shutil
import tempfile
import typing

if typing.TYPE_CHECKING:  # pragma: no cover
    import pathlib

# 3rd party import
import polars

# project import

logger = logging.getLogger("struct.genotypes")


def __manage_group(
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

        __write_or_add(group, path)

        return polars.DataFrame()

    return inner


def __write_or_add(new_lf: polars.DataFrame, partition_path: pathlib.Path) -> None:
    """Create or add new data in parquet partition.

    Args:
        new_lf: Dataframe to add or write
        partition_path: Path where dataframe is write

    Returns:
        None
    """
    if partition_path.exists():
        tmp_file = tempfile.NamedTemporaryFile().name

        old_lf = polars.scan_parquet(partition_path)
        polars.concat([old_lf, new_lf.lazy()]).sink_parquet(tmp_file)

        shutil.move(tmp_file, partition_path)
    else:
        partition_path.parent.mkdir(parents=True, exist_ok=True)
        new_lf.write_parquet(partition_path)


def __hive_worker(lfs: tuple[polars.LazyFrame], output_prefix: pathlib.Path) -> None:
    """Subprocess of hive function run in parallel.

    Args:
        lfs: List of [polars.LazyFrame] you want reorganise
        output_prefix: prefix of hive

    Returns:
        None
    """
    basename = multiprocessing.current_process().name.split("-")[-1]

    polars.concat(lf for lf in lfs if lf is not None).with_columns(
        [
            polars.col("id").mod(256).alias("id_mod"),
        ],
    ).groupby(
        "id_mod",
    ).apply(
        __manage_group(
            output_prefix,
            [
                "id_mod",
            ],
            str(basename),
        ),
        schema={},
    ).collect()


def hive(paths: list[pathlib.Path], output_prefix: pathlib.Path, threads: int, file_per_thread: int) -> None:
    r"""Read all genotypes parquet file and use information to generate a hive like struct, based on $id\ \%\ 256$  with genotype information.

    Real number of threads use are equal to $min(threads, len(paths))$.

    Output format look like: `{output_prefix}/id_mod=[0..255]/[0..threads].parquet`.

    Args:
        paths: list of file you want reorganize
        output_prefix: prefix of hive
        threads: number of multiprocessing threads run
        file_per_thread: number of file manage per multiprocessing threads

    Returns:
        None
    """
    if len(paths) == 0:
        return

    lfs = (polars.scan_parquet(path) for path in paths)

    lf_groups = itertools.zip_longest(
        *[iter(lfs)] * file_per_thread,
        fillvalue=None,
    )

    with multiprocessing.get_context("spawn").Pool(threads) as pool:
        pool.starmap(
            __hive_worker,
            [(lf_group, output_prefix) for lf_group in lf_groups],
        )

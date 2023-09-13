"""Function relate to genotype structuration."""

# std import
from __future__ import annotations

import itertools
import logging
import multiprocessing
import typing

if typing.TYPE_CHECKING:  # pragma: no cover
    import pathlib

# 3rd party import
import polars

# project import

logger = logging.getLogger("struct.genotypes")


def __hive_worker(lfs: tuple[polars.LazyFrame], basename: str, output_prefix: pathlib.Path) -> None:
    """Concatenate multiple parquet file and group it by id % 256.

    Args:
        lfs: List of [polars.LazyFrame] you want reorganise
        basename: name of file
        output_prefix: prefix of hive

    Returns:
        None
    """
    lf = polars.concat(lf for lf in lfs if lf is not None).with_columns(
        [
            polars.col("id").mod(256).alias("id_mod"),
        ],
    )

    for id_mod in range(256):
        path = output_prefix / f"{id_mod}" / f"{basename}.parquet"
        path.parent.mkdir(parents=True, exist_ok=True)
        lf.filter(polars.col("id_mod") == id_mod).sink_parquet(path)


def __merge_file(prefix: pathlib.Path, basenames: list[str]) -> None:
    """Subprocess that merge file generate by __id_spliting.

    Args:
        prefix: prefix of hive struct
        basenames: list of all basenames

    Returns:
        None
    """
    lfs = [
        polars.scan_parquet(prefix / f"{basename}.parquet")
        for basename in basenames
        if (prefix / f"{basename}.parquet").is_file()
    ]

    if lfs:
        polars.concat(lfs).sink_parquet(prefix / "0.parquet")

    for basename in basenames:
        (prefix / f"{basename}.parquet").unlink(missing_ok=True)


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

    path_groups: typing.Iterable[typing.Iterable[pathlib.Path]] = (
        [[path] for path in paths] if file_per_thread <= 1 else itertools.zip_longest(*[iter(paths)] * file_per_thread)
    )

    basenames = ["_".join(p.stem for p in g_paths if p is not None) for g_paths in path_groups]

    lf_groups = [[polars.scan_parquet(p) for p in g_paths] for g_paths in path_groups]

    with multiprocessing.get_context("spawn").Pool(threads) as pool:
        pool.starmap(
            __hive_worker,
            [(lf_group, basename, output_prefix) for lf_group, basename in zip(lf_groups, basenames)],
        )
        pool.starmap(
            __merge_file,
            [(output_prefix / str(id_mod), basenames) for id_mod in range(256)],
        )

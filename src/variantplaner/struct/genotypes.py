"""Function relate to genotype structuration."""

# std import
from __future__ import annotations

import itertools
import logging
import multiprocessing
import shutil
import typing

if typing.TYPE_CHECKING:  # pragma: no cover
    import pathlib

# 3rd party import
import polars

# project import
from variantplaner import normalization

logger = logging.getLogger("struct.genotypes")

NUMBER_OF_BITS = 8
"""Number of high weight bit used in id to perform partitioning."""


def __hive_worker(lfs: tuple[polars.LazyFrame], basename: str, output_prefix: pathlib.Path) -> None:
    """Concatenate several parquet files and group them according to the bits between the 63rd and 55th bits included.

    Args:
        lfs: List of [polars.LazyFrame] you want reorganise
        basename: name of file
        output_prefix: prefix of hive

    Returns:
        None
    """
    logger.info(f"Call hive worker {lfs=}, {basename=}, {output_prefix=}")

    lf = normalization.add_id_part(polars.concat(lf for lf in lfs if lf is not None))

    for name, df in lf.collect().group_by(polars.col("id_part")):
        # required for python 3.11 and polars 0.20.4 wired
        fix_name = name[0] if isinstance(name, tuple) else name
        df.write_parquet(output_prefix / f"id_part={fix_name}" / f"{basename}.parquet")


def __merge_file(prefix: pathlib.Path, basenames: list[str], append: bool) -> None:  # noqa: FBT001
    """Subprocess that merge file generate by __id_spliting.

    Args:
        prefix: prefix of hive struct
        basenames: list of all basenames

    Returns:
        None
    """
    logger.info(f"Call merge file {prefix=}, {basenames=} {append=}")

    paths = [prefix / f"{basename}.parquet" for basename in basenames]
    if append and (prefix / "0.parquet").is_file():
        shutil.copyfile(prefix / "0.parquet", prefix / "0.bkp.parquet")
        paths.append(prefix / "0.bkp.parquet")

    logger.info(f"{paths=}")

    lfs = [polars.scan_parquet(path, hive_partitioning=False) for path in paths if path.is_file()]

    logger.info(f"{lfs=}")
    if lfs:
        logger.info(f"Merge multiple file in {prefix / '0.parquet'}")
        lf = polars.concat(lfs)
        lf = lf.unique(subset=["id", "sample", "gt"])
        lf.sink_parquet(prefix / "0.parquet", maintain_order=False)

    for path in paths:
        logger.info(f"Remove file {path}.parquet")
        path.unlink(missing_ok=True)


def hive(
    paths: list[pathlib.Path],
    output_prefix: pathlib.Path,
    threads: int,
    file_per_thread: int,
    *,
    append: bool,
) -> None:
    r"""Read all genotypes parquet file and use information to generate a hive like struct, based on 63rd and 55th bits included of variant id with genotype information.

    Real number of threads use are equal to $min(threads, len(paths))$.

    Output format look like: `{output_prefix}/id_part=[0..255]/0.parquet`.

    Args:
        paths: list of file you want reorganize
        output_prefix: prefix of hive
        threads: number of multiprocessing threads run
        file_per_thread: number of file manage per multiprocessing threads

    Returns:
        None
    """
    logger.info(f"{paths=} {output_prefix=}, {threads=}, {file_per_thread=}, {append=}")

    if len(paths) == 0:
        return

    for i in range(pow(2, NUMBER_OF_BITS)):
        (output_prefix / f"id_part={i}").mkdir(parents=True, exist_ok=True)

    path_groups: typing.Iterable[typing.Iterable[pathlib.Path]] = list(
        [[path] for path in paths]
        if file_per_thread < 2  # noqa: PLR2004 if number of file is lower than 2 file grouping isn't required
        else itertools.zip_longest(
            *[iter(paths)] * file_per_thread,
        ),
    )

    basenames = ["_".join(p.stem for p in g_paths if p is not None) for g_paths in path_groups]

    lf_groups = [[polars.scan_parquet(p) for p in g_paths if p is not None] for g_paths in path_groups]

    logger.info(f"{path_groups=}, {basenames=}")

    with multiprocessing.get_context("spawn").Pool(threads) as pool:
        pool.starmap(
            __hive_worker,
            [(lf_group, basename, output_prefix) for lf_group, basename in zip(lf_groups, basenames)],
        )

        pool.starmap(
            __merge_file,
            [(output_prefix / f"id_part={id_part}", basenames, append) for id_part in range(pow(2, NUMBER_OF_BITS))],
        )

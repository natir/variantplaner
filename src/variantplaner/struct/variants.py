"""Function relate to vcf structuration."""

# std import
from __future__ import annotations

import logging
import multiprocessing
import os
import pathlib
import random
import shutil
import string
import tempfile
import typing

# 3rd party import
import polars

# project import

logger = logging.getLogger("struct.variants")


def __random_string(size: int = 30) -> str:
    """Generate a random string of 30 ascii character.

    Args:
        size: length of string

    Returns:
        A random string ascii character
    """
    # We didn't use this random for cryptographics usage
    return "".join(random.choice(string.ascii_letters) for _ in range(size))  # noqa: S311


def __chunk_by_memory(
    paths: list[pathlib.Path],
    bytes_limit: int = 10_000_000_000,
) -> typing.Iterator[list[pathlib.Path]]:
    """A generator of chunk iterator on a list of path.

    Each chunk of file have bytes size just upper than bytes_limit.

    Args:
        paths: List of file you want chunked
        bytes_limit: Size limit of chunk file

    Returns:
        An iterator of chunk of path
    """
    total_bytes = 0
    ret = []
    for path in paths:
        total_bytes += int(path.stat().st_size)
        ret.append(path)

        if total_bytes > bytes_limit and len(ret) > 1:
            yield ret
            ret = []
            total_bytes = 0

    yield ret


def __concat_uniq(paths: list[pathlib.Path], output: pathlib.Path) -> None:
    """Merge multiple parquet file.

    Only one copy of each variants is kept, based on id.

    Args:
        paths: List of file you want chunked.
        output: Path where variants is write.

    Returns:
        None
    """
    lf = polars.concat([polars.scan_parquet(path) for path in paths])

    lf = lf.unique(subset=("chr", "pos", "ref", "alt"))

    lf.sink_parquet(output)


def merge(
    paths: list[pathlib.Path],
    output: pathlib.Path,
    memory_limit: int = 10_000_000_000,
    polars_threads: int = 4,
) -> None:
    """Perform merge of multiple parquet variants file in one file.

    These function generate temporary file, by default file are written in `/tmp` but you can control where these files are written by set TMPDIR, TEMP or TMP directory.

    Args:
        paths: List of file you want chunked.
        output: Path where variants is written.
        memory_limit: Size of each chunk in bytes.

    Returns:
        None
    """
    all_threads = int(os.environ["POLARS_MAX_THREADS"])
    multi_threads = max(all_threads // polars_threads, 1)
    os.environ["POLARS_MAX_THREADS"] = str(polars_threads)

    inputs = paths
    temp_directory = tempfile.TemporaryDirectory()
    temp_prefix = pathlib.Path(temp_directory.name)
    logger.debug(f"{temp_prefix=}")

    while len(inputs) != 1:
        new_inputs = []

        inputs_outputs = []
        for input_chunk in __chunk_by_memory(inputs, bytes_limit=memory_limit):
            logger.debug(f"{len(input_chunk)=}")
            if len(input_chunk) > 1:
                # general case
                temp_output = temp_prefix / __random_string()

                new_inputs.append(temp_output)
                inputs_outputs.append((input_chunk, temp_output))

            elif len(input_chunk) == 1:
                # if chunk containt only one file it's last file of inputs
                # we add it to new_inputs list
                new_inputs.append(input_chunk[0])

            logger.debug(f"{new_inputs=}")
            inputs = new_inputs

        with multiprocessing.get_context("spawn").Pool(multi_threads) as pool:
            pool.starmap(__concat_uniq, inputs_outputs)

    # When loop finish we have one file in inputs with all merge
    # We just have to rename it
    shutil.move(inputs[0], output)

    # Call cleanup to remove all tempfile generate durring merging
    temp_directory.cleanup()

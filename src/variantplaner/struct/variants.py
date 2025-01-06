"""Function relate to vcf structuration."""

# std import
from __future__ import annotations

import logging
import multiprocessing
import os
import pathlib
import shutil
import tempfile
import typing

# 3rd party import
import polars

# project import
import variantplaner

logger = logging.getLogger("struct.variants")


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

    if len(ret) > 0:
        yield ret


def __merge_split_unique(paths: list[pathlib.Path], out_prefix: pathlib.Path) -> set[str]:
    """Merge paths input, split chromosone, perform unique write result in out_prefix."""
    logger.info(f"{paths=} {out_prefix}")

    chr_names = set()

    lf = polars.concat([polars.scan_parquet(path) for path in paths])
    for (chr_name, *_), df_group in lf.collect().group_by(polars.col("chr")):
        chr_names.add(str(chr_name))
        df_group.unique(subset=("pos", "ref", "alt")).write_parquet(out_prefix / f"{chr_name}.parquet")

    return chr_names


def __merge_unique(paths: list[pathlib.Path], output: pathlib.Path) -> None:
    """Merge multiple parquet file.

    Only one copy of each variants is kept, based on id.

    Args:
        paths: List of file you want chunked.
        output: Path where variants is write.

    Returns:
        None
    """
    logger.info(f"{paths=} {output=}")

    lf = polars.concat([polars.scan_parquet(path) for path in paths])
    lf = lf.unique(subset=("pos", "ref", "alt"))

    lf.sink_parquet(output)


def merge(
    paths: list[pathlib.Path],
    output_prefix: pathlib.Path,
    memory_limit: int = 10_000_000_000,
    polars_threads: int = 4,
    *,
    append: bool,
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
    output_prefix.mkdir(parents=True, exist_ok=True)

    temp_prefix = pathlib.Path(tempfile.gettempdir()) / "variantplaner" / variantplaner.int2string(hash(output_prefix))
    temp_prefix.mkdir(parents=True, exist_ok=True)

    # merge file -> split by chromosome perform unique
    logger.debug("Start split first file")
    base_inputs_outputs: list[tuple[list[pathlib.Path], pathlib.Path]] = []
    for input_chunk in __chunk_by_memory(paths, bytes_limit=memory_limit):
        local_out_prefix = temp_prefix / variantplaner.int2string(hash(tuple(input_chunk)))
        local_out_prefix.mkdir(parents=True, exist_ok=True)

        base_inputs_outputs.append((input_chunk, local_out_prefix))

    with multiprocessing.get_context("spawn").Pool(multi_threads) as pool:
        chr_names = set().union(*pool.starmap(__merge_split_unique, base_inputs_outputs))
    logger.debug("End split first first file")

    if append and output_prefix.exists():
        chr_names |= {entry.name.split(".")[0] for entry in os.scandir(output_prefix) if entry.is_file()}

    # iterate over chromosme
    logger.debug("Start merge by chromosome")
    for chr_name in chr_names:
        logger.debug(f"start chromosome: {chr_name}")

        chr_temp_prefix = temp_prefix / chr_name
        chr_temp_prefix.mkdir(parents=True, exist_ok=True)

        inputs = [
            path / f"{chr_name}.parquet"
            for (_, path) in base_inputs_outputs
            if (path / f"{chr_name}.parquet").is_file()
        ]

        if append and (output_prefix / f"{chr_name}.parquet").exists():
            inputs.append(output_prefix / f"{chr_name}.parquet")

        if not inputs:
            continue

        while len(inputs) > 1:
            new_inputs = []

            inputs_outputs = []
            for input_chunk in __chunk_by_memory(inputs, bytes_limit=memory_limit):
                logger.debug(f"{input_chunk}")
                if len(input_chunk) == 1:
                    new_inputs.append(input_chunk[0])
                elif len(input_chunk) > 1:
                    temp_output = (
                        chr_temp_prefix / variantplaner.int2string(hash(tuple(input_chunk))) / f"{chr_name}.parquet"
                    )
                    temp_output.parent.mkdir(parents=True, exist_ok=True)

                    new_inputs.append(temp_output)
                    inputs_outputs.append((input_chunk, temp_output))

            inputs = new_inputs

            with multiprocessing.get_context("spawn").Pool(multi_threads) as pool:
                pool.starmap(__merge_unique, inputs_outputs)

        shutil.move(inputs[0], output_prefix / f"{chr_name}.parquet")
        logger.debug(f"end chromosome: {chr_name}")
    logger.debug("End merge by chromosome")

    # Call cleanup to remove all tempfile generate durring merging
    logger.debug("Star clean tmp file")
    shutil.rmtree(temp_prefix, ignore_errors=True)
    logger.debug("End clean tmp file")

"""Benchmark of struct."""

# std import
from __future__ import annotations

import os
import pathlib
import random
import shutil
import tempfile
import typing

# 3rd party import
import polars
import pytest

# project import
from variantplaner import normalization, struct

if typing.TYPE_CHECKING:  # pragma: no cover
    import pytest_benchmark


DATA_DIR = pathlib.Path(__file__).parent.parent / "tests" / "data"


def __generate_variant(common: int = 10_000, diff: int = 1000) -> polars.LazyFrame:
    """Generate a LazyFrame with many duplicate data and some not."""
    nuc = ["A", "C", "T", "G"]

    lf1 = polars.LazyFrame(
        {
            "chr": random.choices(range(1, 25), k=common),
            "pos": random.choices(range(1, 2**32), k=common),
            "ref": random.choices(nuc, k=common),
            "alt": random.choices(nuc, k=common),
        },
    )

    lf2 = polars.LazyFrame(
        {
            "chr": random.choices(range(1, 25), k=diff),
            "pos": random.choices(range(1, 2**32), k=diff),
            "ref": random.choices(nuc, k=diff),
            "alt": random.choices(nuc, k=diff),
        },
    )

    lf3 = polars.LazyFrame(
        {
            "chr": random.choices(range(1, 25), k=diff),
            "pos": random.choices(range(1, 2**32), k=diff),
            "ref": random.choices(nuc, k=diff),
            "alt": random.choices(nuc, k=diff),
        },
    )

    return normalization.add_variant_id(polars.concat([lf1, lf2, lf1, lf3, lf1]))


def __generate_variant_merge(
    common: int,
    merge_type: str,
) -> typing.Callable[[pytest_benchmark.BenchmarkSession], None]:
    @pytest.mark.benchmark(group=f"merge_{merge_type}")
    def inner(benchmark: pytest_benchmark.BenchmarkSession) -> None:
        """Merge variant."""
        lf = __generate_variant(common, common // 10)

        if merge_type == "id":
            benchmark(lambda: lf.unique(subset="id"))
        else:
            benchmark(lambda: lf.unique(subset=("chr", "pos", "ref", "alt")))

    if merge_type == "id":
        inner.__doc__ = f"""Merge variant on id with {common} variant and {common / 10 * 2} different."""
    else:
        inner.__doc__ = f"""Merge variant on variant with {common} variant and {common / 10 * 2} different."""

    return inner


def __generate_variant_merge_on_disk(
    threads: int,
    nb_file: int,
) -> typing.Callable[[pathlib.Path, pytest_benchmark.BenchmarkSession], None]:
    @pytest.mark.benchmark(group="merge_variant_on_disk")
    def inner(tmp_path: pathlib.Path, benchmark: pytest_benchmark.BenchmarkSession) -> None:
        """Merge variant on disk."""
        paths = []
        for i in range(nb_file):
            paths.append(tmp_path / f"{i}.parquet")
            __generate_variant().sink_parquet(tmp_path / f"{i}.parquet")

        os.environ["POLARS_MAX_THREADS"] = str(threads)
        benchmark(lambda: struct.variants.merge(paths, tmp_path / "output.parquet", memory_limit=5_000_000))

    return inner


def __generate_variant_merge_on_disk_old(
    threads: int,
    nb_file: int,
) -> typing.Callable[[pathlib.Path, pytest_benchmark.BenchmarkSession], None]:
    @pytest.mark.benchmark(group="merge_variant_on_disk_old")
    def inner(tmp_path: pathlib.Path, benchmark: pytest_benchmark.BenchmarkSession) -> None:
        """Merge variant on disk."""
        paths = []
        for i in range(nb_file):
            paths.append(tmp_path / f"{i}.parquet")
            __generate_variant().sink_parquet(tmp_path / f"{i}.parquet")

        os.environ["POLARS_MAX_THREADS"] = str(threads)
        benchmark(lambda: __old_merge(paths, tmp_path / "output.parquet", memory_limit=5_000_000))

    return inner


def __old_merge(paths: list[pathlib.Path], output: pathlib.Path, memory_limit: int = 10_000_000_000) -> None:
    """Perform merge of multiple parquet variants file in one file.

    These function generate temporary file, by default file are write in `/tmp` but you can control where these files are write by set TMPDIR, TEMP or TMP directory.

    Args:
        paths: List of file you want chunked.
        output: Path where variants is write.
        memory_limit: Size of each chunk in bytes.

    Returns:
        None
    """
    inputs = paths
    temp_directory = tempfile.TemporaryDirectory()
    temp_prefix = pathlib.Path(temp_directory.name)

    while len(inputs) != 1:
        new_inputs = []

        for input_chunk in struct.variants.__chunk_by_memory(inputs, bytes_limit=memory_limit):
            if len(input_chunk) > 1:
                # general case
                temp_output = temp_prefix / struct.variants.__random_string()

                new_inputs.append(temp_output)
                struct.variants.__concat_uniq(input_chunk, temp_output)

            elif len(input_chunk) == 1:
                # if chunk containt only one file it's last file of inputs
                # we add it to new_inputs list
                new_inputs.append(input_chunk[0])

            inputs = new_inputs

    # When loop finish we have one file in inputs with all merge
    # We just have to rename it
    shutil.move(inputs[0], output)

    # Call cleanup to remove all tempfile generate durring merging
    temp_directory.cleanup()


for i in range(15, 25):
    common = 2**i
    globals()[f"variant_merge_id_{common}"] = __generate_variant_merge(common, "id")
    globals()[f"variant_merge_variant_{common}"] = __generate_variant_merge(common, "variant")


for i in range(20, 41, 2):
    globals()[f"variant_merge_on_disk_{i}"] = __generate_variant_merge_on_disk(8, i)
    globals()[f"variant_merge_on_disk_old_{i}"] = __generate_variant_merge_on_disk(8, i)

"""Tests for help message of cli."""

# std import
from __future__ import annotations

# 3rd party import
from click.testing import CliRunner

# project import
from variantplanner import cli


def test_show_help() -> None:
    """Call cli '--help'."""
    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli.main, ["--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplanner [OPTIONS] COMMAND [ARGS]...

  Run VariantPlanner.

Options:
  -t, --threads INTEGER RANGE  Number of threads usable  [default: 1; x>=0]
  -v, --verbose                Verbosity level  [0<=x<=4]
  --help                       Show this message and exit.

Commands:
  annotations  Convert an annotation variation file in parquet.
  merge        Merge multiple variants parquet file in one.
  metadata     Convert an metadata file in parquet.
  vcf2parquet  Convert a vcf in parquet.
"""
    )


def test_show_help_annotations() -> None:
    """Call cli 'annotations --help'."""
    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli.main, ["annotations", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplanner annotations [OPTIONS] COMMAND [ARGS]...

  Convert an annotation variation file in parquet.

Options:
  -i, --input-paths FILE  Path to input file  [required]
  -o, --output-path PATH  Path where variants annotations will be written in
                          parquet  [required]
  --help                  Show this message and exit.

Commands:
  csv  Convert an annotated csv in parquet.
  vcf  Convert an annotated vcf in parquet.
"""
    )


def test_show_help_merge() -> None:
    """Call cli 'merge --help'."""
    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli.main, ["merge", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplanner merge [OPTIONS]

  Merge multiple variants parquet file in one.

  If you set TMPDIR, TEMP or TMP environment variable you can control where temp
  file is create.

Options:
  -i, --inputs-path FILE          Paths of the variant files to be merged
                                  [required]
  -m, --merge-path PATH           Path where merged variants will be written
                                  [required]
  -b, --bytes-memory-limit INTEGER RANGE
                                  Number of bytes used to build chunk of merge
                                  file  [default: 10000000000; x>=0]
  --help                          Show this message and exit.
"""
    )


def test_show_help_metadata() -> None:
    """Call cli 'metadata --help'."""
    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli.main, ["metadata", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplanner metadata [OPTIONS] COMMAND [ARGS]...

  Convert an metadata file in parquet.

Options:
  -i, --input-path FILE   Path to input file  [required]
  -o, --output-path PATH  Path where variants annotations will be written in
                          parquet  [required]
  --help                  Show this message and exit.

Commands:
  csv   Convert an metadata csv in parquet.
  json  Convert an metadata json in parquet.
"""
    )


def test_show_help_vcf2parquet() -> None:
    """Call cli 'vcf2parquet --help'."""
    runner = CliRunner()
    result = runner.invoke(cli.main, ["vcf2parquet", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplanner vcf2parquet [OPTIONS]

  Convert a vcf in parquet.

Options:
  -i, --input-path FILE  Path to vcf input file  [required]
  -v, --variants FILE    Path where the variants will be written in parquet
                         [required]
  -g, --genotypes FILE   Path where the genotypes will be written in parquet
  --help                 Show this message and exit.
"""
    )

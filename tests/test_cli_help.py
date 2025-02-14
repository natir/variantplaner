"""Tests for help message of cli."""

# std import
from __future__ import annotations

import os

# 3rd party import
import pytest
from click.testing import CliRunner

# project import
from variantplaner import cli


def test_show_help() -> None:
    """Call cli '--help'."""
    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli.main, ["--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplaner [OPTIONS] COMMAND [ARGS]...

  Run VariantPlanner.

Options:
  -t, --threads INTEGER RANGE  Number of threads usable  [default: 1; x>=0]
  -v, --verbose                Verbosity level  [0<=x<=4]
  --debug-info                 Get debug information
  -h, --help                   Show this message and exit.

Commands:
  metadata      Convert metadata file in parquet file.
  parquet2vcf   Convert variant parquet in vcf.
  struct        Subcommand to made struct operation on parquet file.
  transmission  Generate transmission of a genotype set.
  vcf2parquet   Convert a vcf in parquet.
"""
    )


@pytest.mark.skipif(
    os.environ.get("GITHUB_REPOSITORY", default="") == "SeqOIA-IT/variantplaner",
    reason="this test failled in github action",
)
def test_show_help_struct() -> None:
    """Call cli 'struct --help'."""
    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli.main, ["struct", "--help"])

    assert result.exit_code == 0

    # Second check are required due to inconsistency in cli module
    assert (
        result.stdout
        == """Usage: variantplaner struct [OPTIONS] COMMAND1 [ARGS]... [COMMAND2 [ARGS]...]...

  Subcommand to made struct operation on parquet file.

Options:
  -i, --input-paths FILE  Paths of the variant files to be merged.  [required]
  -a, --append            Switch in append mode.
  -h, --help              Show this message and exit.

Commands:
  genotypes  Convert set of genotype parquet in hive like files structures.
  variants   Merge multiple variants parquet file in one.
"""
    )


def test_show_help_metadata() -> None:
    """Call cli 'metadata --help'."""
    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli.main, ["metadata", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplaner metadata [OPTIONS]

  Convert metadata file in parquet file.

Options:
  -i, --input-path PATH           Input path.  [required]
  -o, --output-path PATH          Output path.  [required]
  -t, --input-type [csv|tsv|ljson|json]
                                  Type of input file.  [required]
  -h, --help                      Show this message and exit.
"""
    )


def test_show_help_vcf2parquet() -> None:
    """Call cli 'vcf2parquet --help'."""
    runner = CliRunner()
    result = runner.invoke(cli.main, ["vcf2parquet", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplaner vcf2parquet [OPTIONS] COMMAND1 [ARGS]... [COMMAND2
                                 [ARGS]...]...

  Convert a vcf in parquet.

Options:
  -i, --input-path FILE         Path to vcf input file.  [required]
  -c, --chrom2length-path FILE  CSV file that associates a chromosome name with
                                its size.
  -a, --append                  Switch in append mode.
  -s, --keep-star BOOLEAN       Keep variant with * in alt.
  -h, --help                    Show this message and exit.

Commands:
  annotations  Write annotations.
  genotypes    Write genotypes.
  headers      Write vcf headers.
  variants     Write variants.
"""
    )


def test_show_help_parquet2vcf() -> None:
    """Call cli 'parquet2vcf --help'."""
    runner = CliRunner()
    result = runner.invoke(cli.main, ["parquet2vcf", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplaner parquet2vcf [OPTIONS]

  Convert variant parquet in vcf.

Options:
  -v, --variants-path FILE      Path to variant parquet.  [required]
  -o, --output-path FILE        Path where the vcf is written.  [required]
  -g, --genotypes-path FILE     Path to genotype parquet.
  -H, --headers-path FILE       Path to vcf header.
  -s, --select-chromosome TEXT  Output only record where chromosome match target
  -c, --chromosome TEXT         Name of chromosome column.  [default: chr]
  -p, --position TEXT           Name of position column.  [default: pos]
  -I, --identifier TEXT         Name of identity column.  [default: id]
  -r, --reference TEXT          Name of reference column.  [default: ref]
  -a, --alternative TEXT        Name of alternative column.  [default: alt]
  -q, --quality TEXT            Name of quality column.
  -f, --filter TEXT             Name of filter column.
  -F, --format TEXT             Value of format column.
  -h, --help                    Show this message and exit.
"""
    )


def test_show_help_transmission() -> None:
    """Call cli 'transmission --help'."""
    runner = CliRunner()
    result = runner.invoke(cli.main, ["transmission", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplaner transmission [OPTIONS]

  Generate transmission of a genotype set.

Options:
  -g, --genotypes-path PATH  Path to genotypes parquet.  [required]
  -p, --pedigree-path PATH   Path to pedigree file.
  -i, --index TEXT           Sample name of index.
  -m, --mother TEXT          Sample name of mother.
  -f, --father TEXT          Sample name of father.
  -o, --output-path PATH     Path where transmission will be write.
  -h, --help                 Show this message and exit.
"""
    )

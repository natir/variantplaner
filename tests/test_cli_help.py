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
  metadata     Convert an metadata file in parquet.
  parquet2vcf  Convert parquet in vcf.
  struct       Struct operation on parquet file.
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


def test_show_help_struct() -> None:
    """Call cli 'struct --help'."""
    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli.main, ["struct", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplanner struct [OPTIONS] COMMAND [ARGS]...

  Struct operation on parquet file.

Options:
  -i, --input-paths FILE  Paths of the variant files to be merged  [required]
  --help                  Show this message and exit.

Commands:
  genotypes  Convert set of genotype parquet in hive file.
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


def test_show_help_parquet2vcf() -> None:
    """Call cli 'parquet2vcf --help'."""
    runner = CliRunner()
    result = runner.invoke(cli.main, ["parquet2vcf", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplanner parquet2vcf [OPTIONS]

  Convert parquet in vcf.

Options:
  -i, --input-path FILE   Path to variants in parquet format  [required]
  -o, --output FILE       Path where the vcf is write  [required]
  -c, --chromosome TEXT   Name of chromosome column  [default: chr]
  -p, --position TEXT     Name of position column  [default: pos]
  -I, --identifier TEXT   Name of identity column  [default: id]
  -r, --reference TEXT    Name of reference column  [default: ref]
  -a, --alternative TEXT  Name of alternative column  [default: alt]
  -q, --quality TEXT      Name of quality column
  -f, --filter TEXT       Name of filter column
  --help                  Show this message and exit.
"""
    )

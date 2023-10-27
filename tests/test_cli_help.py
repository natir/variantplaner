"""Tests for help message of cli."""

# std import
from __future__ import annotations

# 3rd party import
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
  --help                       Show this message and exit.

Commands:
  annotations  Convert an annotation variation file in a compatible parquet.
  generate     Generate information.
  metadata     Convert metadata file in parquet file.
  parquet2vcf  Convert some parquet file in vcf.
  struct       Subcommand to made struct operation on parquet file.
  vcf2parquet  Convert a vcf in multiple parquet file.
"""
    )


def test_show_help_annotations() -> None:
    """Call cli 'annotations --help'."""
    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli.main, ["annotations", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplaner annotations [OPTIONS] COMMAND [ARGS]...

  Convert an annotation variation file in a compatible parquet.

Options:
  -i, --input-path FILE         Path to input file  [required]
  -o, --output-path PATH        Path where variants annotations will be written
                                in parquet  [required]
  -c, --chrom2length-path FILE  CSV file that associates a chromosome name with
                                its size  [required]
  --help                        Show this message and exit.

Commands:
  csv  Convert annotations store in csv file to compatible parquet file.
  vcf  Convert a vcf file with INFO field in compatible parquet file.
"""
    )


def test_show_help_struct() -> None:
    """Call cli 'struct --help'."""
    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli.main, ["struct", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplaner struct [OPTIONS] COMMAND [ARGS]...

  Subcommand to made struct operation on parquet file.

Options:
  -i, --input-paths FILE  Paths of the variant files to be merged  [required]
  --help                  Show this message and exit.

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
        == """Usage: variantplaner metadata [OPTIONS] COMMAND [ARGS]...

  Convert metadata file in parquet file.

Options:
  -i, --input-path FILE   Path to input file  [required]
  -o, --output-path PATH  Path where variants annotations will be written in
                          parquet  [required]
  --help                  Show this message and exit.

Commands:
  csv   Convert metadata csv file in parquet file.
  json  Convert metadata json file in parquet file.
"""
    )


def test_show_help_vcf2parquet() -> None:
    """Call cli 'vcf2parquet --help'."""
    runner = CliRunner()
    result = runner.invoke(cli.main, ["vcf2parquet", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplaner vcf2parquet [OPTIONS]

  Convert a vcf in multiple parquet file.

Options:
  -i, --input-path FILE         Path to vcf input file  [required]
  -c, --chrom2length-path FILE  CSV file that associates a chromosome name with
                                its size  [required]
  -v, --variants FILE           Path where the variants will be written in
                                parquet  [required]
  -g, --genotypes FILE          Path where the genotypes will be written in
                                parquet
  -a, --annotations FILE        Path where the annotations will be written in
                                parquet (if no info file is empty)
  -f, --format-string TEXT      Value of FORMAT column, line not match with this
                                are ignored  [default: GT:AD:DP:GQ]
  --help                        Show this message and exit.
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

  Convert some parquet file in vcf.

Options:
  -i, --input-path FILE      Path to variants in parquet format  [required]
  -g, --genotypes-path FILE  Path to genotypes in parquet format
  -o, --output FILE          Path where the vcf is written  [required]
  -c, --chromosome TEXT      Name of chromosome column  [default: chr]
  -p, --position TEXT        Name of position column  [default: pos]
  -I, --identifier TEXT      Name of identity column  [default: id]
  -r, --reference TEXT       Name of reference column  [default: ref]
  -a, --alternative TEXT     Name of alternative column  [default: alt]
  -q, --quality TEXT         Name of quality column
  -f, --filter TEXT          Name of filter column
  -F, --format TEXT          Value of format column
  --help                     Show this message and exit.
"""
    )


def test_show_help_generate() -> None:
    """Call cli 'generate --help'."""
    runner = CliRunner()
    result = runner.invoke(cli.main, ["generate", "--help"])

    assert result.exit_code == 0
    assert (
        result.stdout
        == """Usage: variantplaner generate [OPTIONS] COMMAND [ARGS]...

  Generate information.

Options:
  --help  Show this message and exit.

Commands:
  transmission  Generate transmission file from genotypes and pedigree.
"""
    )

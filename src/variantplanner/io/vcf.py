"""Function to manage input and output of vcf file."""

# std import
import pathlib

# 3rd party import
import polars

# project import


def info2expr(input_path: pathlib.Path) -> list[polars.Expr]:
    """Read vcf header to generate a list of polars.Expr to extract info.

    Args:
        input_path: Path to vcf file.

    Returns:
        List of polars expr to parse info columns.
    """
    return []


def format2expr(input_path: pathlib.Path) -> list[polars.Expr]:
    """Read vcf header to generate a list of polars.Expr to extract genotype.

    Args:
        input_path: Path to vcf file.

    Returns:
        List of polars expr to parse format columns.
    """
    return []


def sample_index(input_path: pathlib.Path) -> dict[str, int]:
    """Read vcf header to generate an association map between sample name and index.

    Args:
        input_path: Path to vcf file.

    Returns:
        Map that associate a sample name to is sample index.
    """
    return {}


def into_lazyframe(input_path: pathlib.Path) -> polars.LazyFrame:
    """Read a vcf file and convert it in lazyframe.

    Args:
        input_path: Path to vcf file.

    Returns:
        A lazyframe that containt vcf information.
    """
    return polars.LazyFrame()


def from_lazyframe(lf: polars.LazyFrame, output_path: pathlib.Path) -> None:
    """Write lazyframe in vcf format.

    Args:
        lf: LazyFrame contains information.
        output_path: Path to where vcf to write.

    Returns:
        None
    """
    return

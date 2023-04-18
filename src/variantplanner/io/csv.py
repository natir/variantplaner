"""Function to manage input and output of csv file."""

# std import
import pathlib

# 3rd party import
import polars

# project import


def into_lazyframe(
    input_path: pathlib.Path,
    chromosome_col: str,
    position_col: str,
    reference_col: str,
    alternative_col: str,
    info_cols: list[str],
) -> polars.LazyFrame:
    """Read a csv file and convert it in lazyframe.

    Args:
        input_path: Path to csv.
        chromosome_col: Name of the column that holds the chromosomes.
        position_col: Name of the column that holds the positions.
        reference_col: Name of the column that holds the reference sequence.
        alternative_col: Name of the column that hold the alternative sequence.

    Returns:
        A lazyframe that containt csv information
    """
    return polars.LazyFrame()


def from_lazyframe(lf: polars.LazyFrame, output_path: pathlib.Path) -> None:
    """Write lazy frame in csv format.

    Args:
        lf: LazyFrame must be write.
        output_path: Path where csv to write.

    Returns:
        None
    """
    return

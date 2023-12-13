"""Module contains command line entry point function."""

# std import
from __future__ import annotations

import logging
import os
import sys
import typing

# 3rd party import
import click
import polars

# project import
from variantplaner import debug


class MultipleValueOption(click.Option):
    """Class use to manage option with variadic value."""

    def __init__(self, *args: list[typing.Any], **kwargs: dict[typing.Any, typing.Any]):
        """Intialise click option parser."""
        super(MultipleValueOption, self).__init__(*args,**kwargs) # type: ignore[arg-type] # noqa: UP008  false positive and complexe type
        self._previous_parser_process = None
        self._eat_all_parser = None

    def add_to_parser(self, parser: typing.Any, ctx: typing.Any) -> typing.Any:
        """Add parser."""

        def parser_process(value: typing.Any, state: typing.Any) -> typing.Any:
            # method to hook to the parser.process
            done = False
            value = [value]

            # grab everything up to the next option
            while state.rargs and not done:
                if self._eat_all_parser:
                    for prefix in self._eat_all_parser.prefixes:
                        if not state.rargs or state.rargs[0].startswith(prefix):
                            done = True
                        if not done:
                            value.append(state.rargs.pop(0))

            if len(value) == 1:
                value = value[0]
            else:
                value = tuple(value)
                self.nargs = len(value)

            # call the actual process
            if self._previous_parser_process:
                self._previous_parser_process(value, state)

        retval = super(MultipleValueOption, self).add_to_parser(parser, ctx)  # noqa: UP008 false positive
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break

        return retval


@click.group(name="variantplaner", context_settings={"help_option_names": ["-h", "--help"]})
@click.pass_context
@click.option(
    "-t",
    "--threads",
    help="Number of threads usable",
    default=1,
    type=click.IntRange(0),
    show_default=True,
)
@click.option("-v", "--verbose", help="Verbosity level", count=True, type=click.IntRange(0, 4))
@click.option(
    "--debug-info",
    help="Get debug information",
    is_flag=True,
    show_default=True,
    default=False,
)
def main(ctx: click.Context, *, threads: int = 1, verbose: int = 0, debug_info: bool = False) -> None:
    """Run VariantPlanner."""
    logging.basicConfig(
        style="{",
        format="{asctime} - {name}:{levelname}: {message}",
        encoding="utf-8",
        level=(4 - verbose) * 10,  # Python choose a strange logging levels order
        stream=sys.stderr,
    )

    logger = logging.getLogger("main")

    logger.debug(f"parameter: {threads=} {verbose=} {debug_info=}")

    if debug_info:
        debug.print_info()
        sys.exit(0)

    ctx.obj = {
        "threads": threads,
    }

    polars.set_random_seed(42)
    os.environ["POLARS_MAX_THREADS"] = str(threads)


# module import required after main definition
from variantplaner.cli import metadata # noqa: E402 F401 I001 these import should be here
from variantplaner.cli import parquet2vcf # noqa: E402 F401  these import should be here
from variantplaner.cli import struct # noqa: E402 F401  these import should be here
from variantplaner.cli import transmission # noqa: E402 F401  these import should be here
from variantplaner.cli import vcf2parquet # noqa: E402 F401  these import should be here

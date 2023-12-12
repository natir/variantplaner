"""Module contains command line entry point function."""

# Why does this file exist, and why not put this in `__main__`?
#
# You might be tempted to import things from `__main__` later,
# but that will cause problems: the code will get executed twice:
#
# - When you run `python -m variantplaner` python will execute
#   `__main__.py` as a script. That means there won't be any
#   `variantplaner.__main__` in `sys.modules`.
# - When you import `__main__` it will get executed again (as a module) because
#   there's no `variantplaner.__main__` in `sys.modules`.

# std import
from __future__ import annotations

import logging
import os
import pathlib
import sys

# 3rd party import
import click
import polars

# project import
from variantplaner import debug, exception, extract, generate, io, struct

class MultipleValueOption(click.Option):
    """Click class to have multiple value option without repeat option name."""

    def __init__(self, *args, **kwargs):
        nargs = kwargs.pop('nargs', -1)
        assert nargs == -1, 'nargs, if set, must be -1 not {}'.format(nargs)
        super(OptionEatAll, self).__init__(*args, **kwargs)
        self._previous_parser_process = None
        self._eat_all_parser = None

    def add_to_parser(self, parser, ctx):
        def parser_process(value, state):
            # method to hook to the parser.process
            done = False
            value = [value]

            # grab everything up to the next option
            while state.rargs and not done:
                for prefix in self._eat_all_parser.prefixes:
                    if state.rargs[0].startswith(prefix):
                        done = True
                    if not done:
                        value.append(state.rargs.pop(0))

            value = tuple(value)

            # call the actual process
            self._previous_parser_process(value, state)

        retval = super(OptionEatAll, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break

        return retval


@click.group(name="variantplaner", context_settings={"help_option_names": ["-h", "--help"]})
@click.option("-t", "--threads", help="Number of threads usable", default=1, type=click.IntRange(0), show_default=True)
@click.option("-v", "--verbose", help="Verbosity level", count=True, type=click.IntRange(0, 4))
@click.option("--debug-info", help="Get debug information", is_flag=True, show_default=True, default=False)
def main(*, threads: int = 1, verbose: int = 0, debug_info: bool = False) -> None:
    """Run VariantPlanner."""
    logging.basicConfig(
        style="{",
        format="{asctime} - {name}:{levelname}: {message}",
        encoding="utf-8",
        level=(4 - verbose) * 10,  # Python choose a bad logging levels order
        stream=sys.stderr,
    )

    logger = logging.getLogger("main")

    logger.debug(f"parameter: {threads=} {verbose=} {debug_info=}")

    if debug_info:
        debug.print_info()
        sys.exit(0)

    polars.set_random_seed(42)
    os.environ["POLARS_MAX_THREADS"] = str(threads)


# module import required after main definition
from .vcf2parquet import *

"""
Main CLI entry point for airsspy.

Provides the ``airss`` command with subcommands for deploying searches,
querying results, checking the environment, and utility tools.
"""

import logging

import click

from airsspy.log import setup_logging

_COMMANDS = {
    "check": (".cmd_check", "check"),
    "convert": (".cmd_convert", "convert"),
    "db": (".cmd_db", "db"),
    "deploy": (".cmd_deploy", "deploy"),
    "pack": (".cmd_pack", "pack"),
    "rank": (".cmd_rank", "rank"),
    "run": (".cmd_run", "run"),
    "tools": (".cmd_tools", "tools"),
    "unpack": (".cmd_pack", "unpack"),
}


class LazyGroup(click.Group):
    """A Click group that imports command modules only when needed."""

    def list_commands(self, ctx):
        return sorted(_COMMANDS)

    def get_command(self, ctx, cmd_name):
        import importlib

        try:
            module_name, attr = _COMMANDS[cmd_name]
        except KeyError:
            return None
        module = importlib.import_module(module_name, __package__)
        return getattr(module, attr)


@click.group("airss", cls=LazyGroup)
@click.version_option(version="0.1.4", prog_name="airsspy")
@click.pass_context
@click.option(
    "-v", "--verbose", count=True, help="Increase verbosity (-v debug, -vv trace)."
)
@click.option(
    "-q",
    "--quiet",
    count=True,
    help="Decrease verbosity (-q warnings, -qq errors, -qqq silent).",
)
@click.option(
    "--db-host",
    default="localhost",
    help="MongoDB host for the jobflow store.",
    show_default=True,
)
@click.option(
    "--db-port",
    default=27017,
    type=int,
    help="MongoDB port for the jobflow store.",
    show_default=True,
)
@click.option(
    "--db-name",
    default="airss",
    help="MongoDB database name.",
    show_default=True,
)
def cli(ctx, verbose, quiet, db_host, db_port, db_name):
    """Command-line interface for AIRSS structure searches."""
    if verbose and quiet:
        raise click.UsageError("--verbose and --quiet are mutually exclusive.")
    level = logging.INFO - (verbose * 10) + (quiet * 10)
    level = max(level, logging.CRITICAL + 10)
    setup_logging(level=level)
    ctx.ensure_object(dict)
    ctx.obj["db_host"] = db_host
    ctx.obj["db_port"] = db_port
    ctx.obj["db_name"] = db_name

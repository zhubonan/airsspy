"""
Main CLI entry point for airsspy.

Provides the ``airss`` command with subcommands for deploying searches,
querying results, checking the environment, and utility tools.
"""

import click

from .cmd_check import check
from .cmd_db import db
from .cmd_deploy import deploy
from .cmd_tools import tools


@click.group("airss")
@click.version_option(version="0.1.4", prog_name="airsspy")
@click.pass_context
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
def cli(ctx, db_host, db_port, db_name):
    """Command-line interface for AIRSS structure searches."""
    ctx.ensure_object(dict)
    ctx.obj["db_host"] = db_host
    ctx.obj["db_port"] = db_port
    ctx.obj["db_name"] = db_name


cli.add_command(deploy)
cli.add_command(db)
cli.add_command(check)
cli.add_command(tools)

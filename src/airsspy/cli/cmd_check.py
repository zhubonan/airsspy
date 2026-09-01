"""
CLI commands for checking the AIRSS environment.
"""

import subprocess

import click

from airsspy.scheduler import Dummy, Scheduler


@click.group("check")
def check():
    """Check environment and installation."""


@check.command("scheduler")
@click.option(
    "--allow-dummy", is_flag=True, default=False, help="Allow dummy scheduler"
)
def check_scheduler(allow_dummy):
    """Check the status of the job scheduler."""
    obj = Scheduler.get_scheduler()
    if isinstance(obj, Dummy) and not allow_dummy:
        raise click.ClickException(
            "Not in a scheduler environment (use --allow-dummy to override)"
        )
    click.echo(f"Scheduler:     {obj}")
    click.echo(f"JOB ID:        {obj.job_id}")
    click.echo(f"USERNAME:      {obj.user_name}")
    click.echo(f"SECONDS LEFT:  {obj.get_remaining_seconds()}")


@check.command("airss")
def check_airss():
    """Check the AIRSS package installation."""
    click.echo("Checking essential AIRSS scripts:")
    is_ok = True
    for name in ["buildcell", "cabal", "castep_relax"]:
        try:
            loc = subprocess.check_output(
                ["which", name], universal_newlines=True
            ).strip()
        except subprocess.CalledProcessError:
            click.echo(f"  {name}: NOT FOUND")
            is_ok = False
        else:
            click.echo(f"  {name}: {loc}")

    if is_ok:
        click.echo("All GOOD!")
    else:
        raise click.ClickException("Some AIRSS components are missing")


@check.command("database")
@click.pass_context
def check_database(ctx):
    """Check the connection to the jobflow store."""
    from airsspy.jf.store import SearchStore

    host = ctx.obj["db_host"]
    port = ctx.obj["db_port"]
    database = ctx.obj["db_name"]

    click.echo(f"Connecting to {host}:{port}/{database}...")
    try:
        store = SearchStore(database=database, host=host, port=port)
        store.connect()
        projects = store.list_projects()
        store.close()
    except Exception as exc:
        raise click.ClickException(f"Cannot connect to the database: {exc}") from exc

    click.echo("Connection successful: OK")
    click.echo(f"Projects found: {len(projects)}")
    for proj in projects:
        click.echo(f"  {proj}")

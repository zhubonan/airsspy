"""
CLI commands for querying the AIRSS jobflow store.
"""

from pathlib import Path

import click


@click.group("db")
@click.pass_context
def db(ctx):
    """Query and manage AIRSS search results."""
    click.echo(
        f"Using database: {ctx.obj['db_host']}:{ctx.obj['db_port']}/{ctx.obj['db_name']}",
        err=True,
    )


@db.command("list-projects")
@click.pass_context
def list_projects(ctx):
    """List all projects in the store."""
    from airsspy.jf.store import SearchStore

    store = SearchStore(
        database=ctx.obj["db_name"],
        host=ctx.obj["db_host"],
        port=ctx.obj["db_port"],
    )
    store.connect()
    projects = store.list_projects()
    store.close()

    for proj in projects:
        click.echo(proj)


@db.command("list-seeds")
@click.option("--project", "-p", help="Filter by project name")
@click.pass_context
def list_seeds(ctx, project):
    """List seed names, optionally filtered by project."""
    from airsspy.jf.store import SearchStore

    store = SearchStore(
        database=ctx.obj["db_name"],
        host=ctx.obj["db_host"],
        port=ctx.obj["db_port"],
    )
    store.connect()
    seeds = store.list_seeds(project_name=project)
    store.close()

    for seed in seeds:
        click.echo(seed)


@db.command("summary")
@click.option("--project", "-p", help="Filter by project name")
@click.pass_context
def summary(ctx, project):
    """Display a summary of structures per project/seed."""
    from tabulate import tabulate

    from airsspy.jf.store import SearchStore

    store = SearchStore(
        database=ctx.obj["db_name"],
        host=ctx.obj["db_host"],
        port=ctx.obj["db_port"],
    )
    store.connect()
    df = store.show_struct_counts(project_name=project)
    store.close()

    if df.empty:
        click.echo("No data available")
        return
    click.echo(tabulate(df, headers="keys", tablefmt="simple", showindex=False))


@db.command("throughput")
@click.option("--past-days", "-d", default=2, type=float, help="Look back N days")
@click.option("--project", "-p", help="Filter by project name")
@click.pass_context
def throughput(ctx, past_days, project):
    """Summarise search throughput over recent days."""
    from tabulate import tabulate

    from airsspy.jf.store import SearchStore

    store = SearchStore(
        database=ctx.obj["db_name"],
        host=ctx.obj["db_host"],
        port=ctx.obj["db_port"],
    )
    store.connect()
    df = store.throughput_summary(past_days=past_days, project_name=project)
    store.close()

    if df.empty:
        click.echo(f"No results found for the past {24 * past_days:.0f} hours.")
        return
    click.echo(tabulate(df, headers="keys", tablefmt="simple", showindex=False))


@db.command("retrieve-project")
@click.option("--project", "-p", required=True, help="Project name")
@click.option("--output", "-o", default=".", help="Output directory for RES files")
@click.pass_context
def retrieve_project(ctx, project, output):
    """Retrieve RES files for a particular project."""
    from tqdm import tqdm

    from airsspy.jf.store import SearchStore

    out_dir = Path(output)
    out_dir.mkdir(parents=True, exist_ok=True)

    store = SearchStore(
        database=ctx.obj["db_name"],
        host=ctx.obj["db_host"],
        port=ctx.obj["db_port"],
    )
    store.connect()
    results = store.retrieve_project(project)
    store.close()

    click.echo(f"Retrieving {len(results)} results for project: {project}")
    for result in tqdm(results):
        if result.res_content:
            (out_dir / f"{result.struct_name}.res").write_text(result.res_content)

    click.echo(f"Written to {out_dir}/")

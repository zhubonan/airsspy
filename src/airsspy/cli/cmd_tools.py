"""
CLI commands for miscellaneous tools.
"""

import click


@click.group("tools")
@click.pass_context
def tools(ctx):
    """Collection of utility tools."""
    _ = ctx


@tools.command("modcell")
@click.argument("base_cell")
@click.argument("other_cell")
def modcell(base_cell, other_cell):
    """Modify the structure of a CELL file using another.

    BASE_CELL is the template cell file. OTHER_CELL is the file
    providing the new structure (any ASE-supported format).
    """
    from ase.io import read

    from airsspy.tools.modcell import modify_cell

    atoms = read(other_cell)
    lines = modify_cell(base_cell, atoms)
    click.echo("\n".join(lines))


def _parse_alpha_grid(raw: str) -> list[float]:
    return [float(token.strip()) for token in raw.split(",") if token.strip()]


@tools.command("volume-minsep-build-dataset")
@click.option(
    "--mp-docs",
    required=True,
    type=click.Path(exists=True),
    help="Materials Project document dump.",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(),
    help="Output raw minsep/volume dataset JSON.",
)
@click.option(
    "--progress/--no-progress",
    default=True,
    help="Show progress while converting MP docs.",
)
def volume_minsep_build_dataset(mp_docs, output, progress):
    """Build a raw minsep/volume dataset from MP documents."""
    import json

    from airsspy.volume_minsep_data import build_minsep_volume_dataset

    summary = build_minsep_volume_dataset(
        mp_docs_path=mp_docs,
        output_path=output,
        show_progress=progress,
    )
    click.echo(json.dumps(summary, indent=2, sort_keys=True))


@tools.command("volume-minsep-curate-dataset")
@click.option(
    "--dataset",
    required=True,
    type=click.Path(exists=True),
    help="Raw minsep/volume dataset JSON.",
)
@click.option(
    "--mp-docs",
    required=True,
    type=click.Path(exists=True),
    help="Materials Project document dump with energy_above_hull.",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(),
    help="Curated output JSON path.",
)
def volume_minsep_curate_dataset(dataset, mp_docs, output):
    """Keep one minsep/volume row per reduced formula."""
    import json

    from airsspy.volume_minsep_data import curate_minsep_volume_dataset

    summary = curate_minsep_volume_dataset(
        minsep_dataset_path=dataset,
        mp_docs_path=mp_docs,
        output_path=output,
    )
    click.echo(json.dumps(summary, indent=2, sort_keys=True))


@tools.command("volume-minsep-train-baseline")
@click.option(
    "--dataset",
    required=True,
    type=click.Path(exists=True),
    help="Raw or curated minsep/volume dataset JSON.",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(),
    help="Directory for trained baseline artifacts.",
)
@click.option("--seed", type=int, default=17, show_default=True)
@click.option("--pair-min-count", type=int, default=1, show_default=True)
@click.option("--volume-alpha-grid", default="1e-6,1e-4,1e-2,1e0,1e2")
@click.option("--minsep-alpha-grid", default="1e-6,1e-4,1e-2,1e0,1e2")
@click.option("--max-volume-rows", type=int, default=None)
@click.option("--max-pair-rows", type=int, default=None)
def volume_minsep_train_baseline(
    dataset,
    output_dir,
    seed,
    pair_min_count,
    volume_alpha_grid,
    minsep_alpha_grid,
    max_volume_rows,
    max_pair_rows,
):
    """Train the lightweight baseline volume/minsep predictor."""
    import json

    from airsspy.volume_minsep_model import train_baseline_models

    summary = train_baseline_models(
        dataset_path=dataset,
        output_dir=output_dir,
        seed=seed,
        pair_min_count=pair_min_count,
        volume_alpha_grid=_parse_alpha_grid(volume_alpha_grid),
        minsep_alpha_grid=_parse_alpha_grid(minsep_alpha_grid),
        max_volume_rows=max_volume_rows,
        max_pair_rows=max_pair_rows,
    )
    click.echo(json.dumps(summary, indent=2, sort_keys=True))

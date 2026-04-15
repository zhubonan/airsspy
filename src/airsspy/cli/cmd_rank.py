"""CLI command for ranking AIRSS structures by energy."""

import sys

import click


@click.command("rank")
@click.argument("files", nargs=-1, type=click.Path(exists=True))
@click.option(
    "-t", "--top", "top_n", type=int, default=None, help="Show only top N structures"
)
@click.option(
    "-de",
    "--delta-e",
    type=float,
    default=None,
    help="Filter structures with relative energy per atom above this threshold (eV)",
)
@click.option(
    "-f", "--formula", default=None, help="Filter by reduced formula (e.g. SiO2)"
)
@click.option(
    "-nr",
    "--absolute",
    is_flag=True,
    help="Show absolute enthalpy for all entries",
)
@click.option("-l", "--long-labels", is_flag=True, help="Show full structure labels")
@click.option(
    "-s",
    "--summary",
    is_flag=True,
    help="Show only the most stable per composition",
)
@click.option(
    "-u",
    "--unite",
    type=float,
    default=None,
    help="Merge similar structures with given distance threshold",
)
@click.option(
    "--input-format",
    default=None,
    type=click.Choice(["res", "extxyz", "auto"]),
    help="Force input format (default: auto-detect from extension)",
)
@click.option(
    "--fingerprint-cutoff",
    type=float,
    default=4.0,
    help="Distance cutoff (Å) for fingerprint computation (default: 4.0)",
)
def rank(
    files,
    top_n,
    delta_e,
    formula,
    absolute,
    long_labels,
    summary,
    unite,
    input_format,
    fingerprint_cutoff,
):
    """Rank structures by enthalpy per formula unit.

    Reads SHELX .res files from stdin and/or FILE arguments.
    For extxyz files, energy is read from atoms.info['energy'].

    Example usage:

    \b
      cat *.res | airss rank
      airss rank *.res -t 20
      airss rank structures.xyz --input-format extxyz
      cat packed.res | airss rank -de 0.05 -f SiO2
      cat *.res | airss rank -u 0.1 -s
    """
    from airsspy.ranking import (
        eliminate_similar,
        format_header,
        format_rank_line,
        rank_structures,
        read_extxyz_file,
        read_res_file,
        read_res_stream,
        summary_structures,
    )

    records = []

    # Read from file arguments
    for fpath in files:
        fmt = input_format
        if fmt is None or fmt == "auto":
            if fpath.endswith((".xyz", ".extxyz")):
                fmt = "extxyz"
            else:
                fmt = "res"

        if fmt == "extxyz":
            records.extend(read_extxyz_file(fpath))
        else:
            records.extend(read_res_file(fpath))

    # Read from stdin (if not a terminal)
    if not sys.stdin.isatty():
        records.extend(read_res_stream(sys.stdin))

    if not records:
        click.echo("No structures found.", err=True)
        return

    click.echo(f"Read {len(records)} structures", err=True)

    # Merge similar structures if requested
    if unite is not None:
        click.echo(f"Merging similar structures (threshold={unite})", err=True)
        records = eliminate_similar(records, unite, cutoff=fingerprint_cutoff)
        click.echo(f"After merging: {len(records)} structures", err=True)

    # Determine if any structure has non-zero spin
    has_spin = any(r.spin != 0.0 or r.spin_abs != 0.0 for r in records)

    if summary:
        ranked, total = summary_structures(records, delta_e=delta_e)
        click.echo(format_header(show_spin=has_spin, summary_mode=True, long_labels=long_labels), err=True)
        for rec in ranked:
            line = format_rank_line(
                rec, long_labels=long_labels, show_spin=has_spin, summary_mode=True
            )
            click.echo(line)
        click.echo(f"Number of structures   : {total}", err=True)
        click.echo(f"Number of compositions : {len(ranked)}", err=True)
    else:
        ranked = rank_structures(
            records,
            delta_e=delta_e,
            formula_filter=formula,
            top_n=top_n,
            absolute=absolute,
        )

        if not ranked:
            click.echo("No structures remaining after filtering.", err=True)
            return

        click.echo(format_header(show_spin=has_spin, long_labels=long_labels), err=True)
        for rec in ranked:
            line = format_rank_line(rec, long_labels=long_labels, show_spin=has_spin)
            click.echo(line)

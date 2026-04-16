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
    "-m",
    "--maxwell",
    is_flag=True,
    help="Compute convex hull (Maxwell construction) using pymatgen PhaseDiagram",
)
@click.option(
    "-el",
    "--element-list",
    default=None,
    help="Comma-separated element list for phase diagram (e.g. Si,O). Auto-detected if not given.",
)
@click.option(
    "--plot",
    "plot_path",
    default=None,
    type=click.Path(),
    help="Save phase diagram plot as HTML to PATH (requires -m)",
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
    maxwell,
    element_list,
    plot_path,
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
      cat *.res | airss rank -m -el Si,O
      cat *.res | airss rank -m -el Si,O --plot hull.html
    """
    from airsspy.ranking import (
        eliminate_similar,
        format_header,
        format_maxwell_header,
        format_maxwell_line,
        format_rank_line,
        maxwell_construction,
        rank_structures,
        read_extxyz_file,
        read_res_file,
        read_res_stream,
        summary_structures,
    )

    if maxwell and formula:
        click.echo("Error: -f/--formula cannot be used with -m/--maxwell", err=True)
        sys.exit(1)

    if plot_path and not maxwell:
        click.echo("Error: --plot requires -m/--maxwell", err=True)
        sys.exit(1)

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

    if maxwell:
        # Parse element list
        elements = None
        if element_list:
            elements = [el.strip() for el in element_list.split(",")]

        try:
            ranked, pd, inferred_elements = maxwell_construction(
                records, elements=elements, delta_e=delta_e
            )
        except ValueError as exc:
            click.echo(f"Error: {exc}", err=True)
            sys.exit(1)

        # Use inferred element names (pd.elements may show DummySpecies for non-standard species)
        elem_names = ",".join(inferred_elements)
        n_stable = sum(1 for r in ranked if r["on_hull"])
        click.echo(
            f"Maxwell construction for {elem_names}: "
            f"{len(ranked)} structures, {n_stable} on hull",
            err=True,
        )

        click.echo(
            format_maxwell_header(show_spin=has_spin, long_labels=long_labels), err=True
        )

        if summary:
            # Summary mode: show only lowest e_above_hull per composition
            best_per_formula: dict[str, dict] = {}
            for rec in ranked:
                f = rec["formula"]
                if f not in best_per_formula or rec["e_above_hull"] < best_per_formula[f]["e_above_hull"]:
                    best_per_formula[f] = rec
            summary_recs = sorted(
                best_per_formula.values(), key=lambda r: r["e_above_hull"]
            )
            for rec in summary_recs:
                line = format_maxwell_line(rec, long_labels=long_labels, show_spin=has_spin)
                click.echo(line)
        else:
            for rec in ranked:
                line = format_maxwell_line(rec, long_labels=long_labels, show_spin=has_spin)
                click.echo(line)

        # Generate plot if requested
        if plot_path:
            try:
                from airsspy.ranking import plot_maxwell

                fig = plot_maxwell(ranked, inferred_elements)
                fig.write_html(str(plot_path))
                click.echo(f"Plot saved to {plot_path}", err=True)
            except Exception as exc:
                click.echo(f"Warning: could not generate plot: {exc}", err=True)

    elif summary:
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

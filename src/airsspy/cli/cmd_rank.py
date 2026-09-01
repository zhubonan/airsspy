"""CLI command for ranking AIRSS structures by energy."""

import sys
from pathlib import Path

import click


def _write_record_file(rec, path: Path, fmt: str) -> None:
    """Write a StructureRecord to a .res or .xyz file."""
    path.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "res":
        if rec._raw_lines:
            lines = [ln.rstrip("\n") for ln in rec._raw_lines]
            if lines and not lines[-1].strip().startswith("END"):
                lines.append("END")
            path.write_text("\n".join(lines) + "\n")
        elif rec._atoms is not None:
            from airsspy.convert import _atoms_to_res_lines

            lines = _atoms_to_res_lines(rec._atoms)
            path.write_text("\n".join(lines) + "\n")
        else:
            click.echo(f"Warning: no data for {rec.label}, skipping", err=True)
    else:
        if rec._atoms is not None:
            from ase.io import write as ase_write

            ase_write(str(path), rec._atoms, format="extxyz")
        elif rec._raw_lines:
            from ase.io import write as ase_write

            from airsspy.restools import RESFile

            res = RESFile.from_lines(rec._raw_lines, include_structure=True)
            if res.atoms is not None:
                atoms = res.atoms
                atoms.info["label"] = rec.label
                atoms.info["pressure"] = rec.pressure
                atoms.info["spin"] = rec.spin
                atoms.info["spin_abs"] = rec.spin_abs
                atoms.info["symm"] = rec.symm
                ase_write(str(path), atoms, format="extxyz")
        else:
            click.echo(f"Warning: no data for {rec.label}, skipping", err=True)


def _echo_pathology_diagnostics(kept_count, rejected, diagnostics) -> None:
    """Print post-search pathological pruning diagnostics."""
    click.echo(
        f"Pathological prune: kept {kept_count}, rejected {len(rejected)}",
        err=True,
    )

    skipped = 0
    for diag in diagnostics:
        if diag.get("status") != "applied":
            skipped += 1
            continue

        message = (
            "  {formula}: cutoff={cutoff:.6f} eV/atom "
            "(median={median:.6f}, robust_sigma={robust_sigma:.6f}, "
            "tail={tail_size}, baseline={baseline_size}, rejected={rejected_count})"
        ).format(**diag)
        click.echo(message, err=True)

    if skipped:
        click.echo(
            f"  skipped {skipped} formula groups with insufficient/zero-MAD statistics",
            err=True,
        )

    if rejected:
        labels = ", ".join(rec.label for rec in rejected)
        click.echo(f"Rejected pathological structures: {labels}", err=True)


@click.command(
    "rank",
    context_settings={"help_option_names": ["-h", "--help", "-?"]},
)
@click.argument("files", nargs=-1, type=click.Path(exists=True))
@click.option(
    "-r",
    "--rank",
    "rank_compat",
    is_flag=True,
    expose_value=False,
    help="Compatibility flag; ranking is the default action.",
)
@click.option(
    "-t", "--top", "top_n", type=int, default=None, help="Show only top N structures"
)
@click.option(
    "-de",
    "--delta-e",
    "--denergy",
    "--delta_e",
    "delta_e",
    type=float,
    default=None,
    help="Filter structures with relative energy per atom above this threshold (eV)",
)
@click.option(
    "-f",
    "--formula",
    default=None,
    help="Filter by formula: exact (SiO2), elements (Si,O), or glob (Si*)",
)
@click.option(
    "--filter-name",
    default=None,
    help="Filter structures by label using glob pattern (e.g. 'BiSI-0*')",
)
@click.option(
    "-fu",
    "--formula-unit",
    "--formula_unit",
    "formula_units",
    type=int,
    default=None,
    help="Filter by exact number of formula units.",
)
@click.option(
    "-sn",
    "--speciesnumber",
    "--species-number",
    "species_number",
    type=int,
    default=None,
    help="Filter by exact number of species.",
)
@click.option(
    "-in",
    "--ionsnumber",
    "--ions-number",
    "ions_number",
    type=int,
    default=None,
    help="Filter by ion count; negative values mean up to abs(N).",
)
@click.option(
    "--prune-pathological",
    is_flag=True,
    help="Remove suspicious low-energy pathological structures before ranking.",
)
@click.option(
    "--pathology-tail-fraction",
    type=float,
    default=0.10,
    show_default=True,
    help="Lowest-energy fraction used for pathological pruning statistics.",
)
@click.option(
    "--pathology-sigma-factor",
    type=float,
    default=3.0,
    show_default=True,
    help="Robust-sigma multiplier for pathological pruning cutoff.",
)
@click.option(
    "--pathology-trim-count",
    type=int,
    default=1,
    show_default=True,
    help="Number of lowest tail structures excluded from baseline statistics.",
)
@click.option(
    "--pathology-min-tail-size",
    type=int,
    default=5,
    show_default=True,
    help="Minimum baseline tail size required for pathological pruning.",
)
@click.option(
    "-nr",
    "--absolute",
    "--not-relative",
    "--not_relative",
    "absolute",
    is_flag=True,
    help="Show absolute enthalpy for all entries",
)
@click.option(
    "-l",
    "--long-labels",
    "--long",
    "long_labels",
    is_flag=True,
    help="Show full structure labels",
)
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
    "--unite-output",
    default=None,
    type=click.Path(),
    help="Write merged structure groups to this directory (requires -u)",
)
@click.option(
    "--unite-output-format",
    type=click.Choice(["res", "extxyz"]),
    default="extxyz",
    help="Output format for united structures (default: extxyz)",
)
@click.option(
    "--unite-ethresh",
    type=float,
    default=0.1,
    help="Energy threshold (eV/atom) for merge candidates above the minimum (default: 0.1)",
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
    "--elementlist",
    "--element_list",
    "element_list",
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
    "-dr",
    "--distance",
    "--fingerprint-cutoff",
    "fingerprint_cutoff",
    type=float,
    default=10.0,
    help="Distance cutoff (Å) for fingerprint computation (default: 10.0)",
)
@click.option(
    "--unite-zweight",
    is_flag=True,
    help="Enable Z-weighting of fingerprint distances when merging (default: unweighted)",
)
@click.option(
    "-p",
    "--pressure",
    "ext_pressure",
    type=float,
    default=None,
    help="Apply external pressure (GPa). Enthalpy becomes H = E + PV.",
)
@click.option(
    "--symprec",
    type=float,
    default=0.01,
    help="Symmetry tolerance for spacegroup detection from extxyz (default: 0.01 Å)",
)
@click.option(
    "--energy-field",
    default=None,
    help="atoms.info key for energy (extxyz only; default: auto-detect)",
)
@click.option(
    "--label-field",
    default=None,
    help="atoms.info key for label (extxyz only; default: auto-detect)",
)
@click.option(
    "--pressure-field",
    default=None,
    help="atoms.info key for pressure in GPa (extxyz only; default: auto-detect)",
)
def rank(
    files,
    top_n,
    delta_e,
    formula,
    filter_name,
    formula_units,
    species_number,
    ions_number,
    prune_pathological,
    pathology_tail_fraction,
    pathology_sigma_factor,
    pathology_trim_count,
    pathology_min_tail_size,
    absolute,
    long_labels,
    summary,
    unite,
    unite_output,
    unite_output_format,
    unite_ethresh,
    maxwell,
    element_list,
    plot_path,
    input_format,
    fingerprint_cutoff,
    unite_zweight,
    ext_pressure,
    symprec,
    energy_field,
    label_field,
    pressure_field,
):
    """Rank structures by enthalpy per formula unit.

    Reads SHELX .res files from stdin and/or FILE arguments.
    For extxyz files, energy, label and pressure are auto-detected
    from atoms.info / calculator.  Override with --energy-field,
    --label-field, --pressure-field.

    Example usage:

    \b
      cat *.res | ap rank
      ap rank *.res -t 20
      ap rank structures.xyz --input-format extxyz
      ap rank structures.xyz --label-field structure_id
      cat packed.res | ap rank -de 0.05 -f SiO2
      cat *.res | ap rank -u 0.1 -s
      cat *.res | ap rank -m -el Si,O
      cat *.res | ap rank -m -el Si,O --plot hull.html
    """
    from airsspy.ranking import (
        apply_external_pressure,
        eliminate_similar,
        fill_dict_symm,
        filter_by_formula_units,
        filter_by_ions_number,
        filter_by_species_number,
        format_header,
        format_maxwell_header,
        format_maxwell_line,
        format_rank_line,
        maxwell_construction,
        prefilter_records,
        prune_pathological_records,
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
    keep_res_raw = unite is not None or bool(unite_output)

    # Read from file arguments
    for fpath in files:
        fmt = input_format
        if fmt is None or fmt == "auto":
            if fpath.endswith((".xyz", ".extxyz")):
                fmt = "extxyz"
            else:
                fmt = "res"

        if fmt == "extxyz":
            records.extend(
                read_extxyz_file(
                    fpath,
                    energy_field=energy_field,
                    label_field=label_field,
                    pressure_field=pressure_field,
                )
            )
        else:
            records.extend(read_res_file(fpath, keep_raw=keep_res_raw))

    # Read from stdin (if not a terminal)
    if not sys.stdin.isatty():
        records.extend(read_res_stream(sys.stdin, keep_raw=keep_res_raw))

    if not records:
        click.echo("No structures found.", err=True)
        return

    click.echo(f"Read {len(records)} structures", err=True)

    # Step 2: Apply external pressure
    if ext_pressure is not None:
        apply_external_pressure(records, ext_pressure)
        click.echo(f"Applied external pressure: {ext_pressure} GPa", err=True)

    # Step 3: Filter by name and formula
    if filter_name is not None:
        from airsspy.ranking import filter_by_name as _filter_by_name

        before = len(records)
        records = _filter_by_name(records, filter_name)
        click.echo(
            f"Filter --filter-name '{filter_name}': {before} → {len(records)}", err=True
        )

    if formula is not None and not maxwell:
        from airsspy.ranking import filter_by_formula as _filter_by_formula

        before = len(records)
        records = _filter_by_formula(records, formula)
        click.echo(f"Filter -f '{formula}': {before} → {len(records)}", err=True)

    if formula_units is not None:
        before = len(records)
        records = filter_by_formula_units(records, formula_units)
        click.echo(
            f"Filter -fu {formula_units}: {before} → {len(records)}",
            err=True,
        )

    if species_number is not None:
        before = len(records)
        records = filter_by_species_number(records, species_number)
        click.echo(
            f"Filter -sn {species_number}: {before} → {len(records)}",
            err=True,
        )

    if ions_number is not None:
        before = len(records)
        records = filter_by_ions_number(records, ions_number)
        click.echo(
            f"Filter -in {ions_number}: {before} → {len(records)}",
            err=True,
        )

    if not records:
        click.echo("No structures remaining after filtering.", err=True)
        return

    if prune_pathological:
        try:
            records, rejected, diagnostics = prune_pathological_records(
                records,
                tail_fraction=pathology_tail_fraction,
                sigma_factor=pathology_sigma_factor,
                trim_count=pathology_trim_count,
                min_tail_size=pathology_min_tail_size,
            )
        except ValueError as exc:
            click.echo(f"Error: {exc}", err=True)
            sys.exit(1)

        _echo_pathology_diagnostics(len(records), rejected, diagnostics)
        if not records:
            click.echo(
                "No structures remaining after pathological pruning.",
                err=True,
            )
            return

    # Step 4: Pre-rank filter (reduce set before merge)
    if unite is not None:
        before = len(records)
        records = prefilter_records(records, ethresh=unite_ethresh)
        click.echo(
            f"Pre-rank filter (ethresh={unite_ethresh}): {before} → {len(records)} structures",
            err=True,
        )

    # Step 5: Merge similar structures
    if unite is not None:
        click.echo(f"Merging similar structures (threshold={unite})", err=True)
        records = eliminate_similar(
            records, unite, cutoff=fingerprint_cutoff, zweight=unite_zweight
        )
        click.echo(f"After merging: {len(records)} structures", err=True)

    # Step 5b: Write united output
    if unite_output and unite is not None:
        out_dir = Path(unite_output)
        out_dir.mkdir(parents=True, exist_ok=True)
        ext = "xyz" if unite_output_format == "extxyz" else "res"
        n_groups = 0
        n_files = 0
        for rec in records:
            group_dir = out_dir / f"{rec.label}.{ext}"
            group_dir.mkdir(exist_ok=True)
            _write_record_file(
                rec, group_dir / f"{rec.label}.{ext}", unite_output_format
            )
            n_files += 1
            for peer in rec._merged_peers:
                _write_record_file(
                    peer, group_dir / f"{peer.label}.{ext}", unite_output_format
                )
                n_files += 1
            n_groups += 1
        click.echo(
            f"Wrote {n_files} structures in {n_groups} groups to {out_dir}/", err=True
        )

    # Step 6: Post-rank and display
    has_spin = any(r.spin != 0.0 or r.spin_abs != 0.0 for r in records)

    if maxwell:
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

        elem_names = ",".join(inferred_elements)
        n_stable = sum(1 for r in ranked if r["on_hull"])
        click.echo(
            f"Maxwell construction for {elem_names}: "
            f"{len(ranked)} structures, {n_stable} on hull",
            err=True,
        )

        if summary:
            best_per_formula: dict[str, dict] = {}
            for rec in ranked:
                f = rec["formula"]
                if (
                    f not in best_per_formula
                    or rec["e_above_hull"] < best_per_formula[f]["e_above_hull"]
                ):
                    best_per_formula[f] = rec
            summary_recs = sorted(
                best_per_formula.values(), key=lambda r: r["e_above_hull"]
            )
            fill_dict_symm(summary_recs, symprec=symprec)
        else:
            fill_dict_symm(ranked, symprec=symprec)

        click.echo(
            format_maxwell_header(show_spin=has_spin, long_labels=long_labels), err=True
        )

        if summary:
            for rec in summary_recs:
                line = format_maxwell_line(
                    rec, long_labels=long_labels, show_spin=has_spin
                )
                click.echo(line)
        else:
            for rec in ranked:
                line = format_maxwell_line(
                    rec, long_labels=long_labels, show_spin=has_spin
                )
                click.echo(line)

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
        fill_dict_symm(ranked, symprec=symprec)
        click.echo(
            format_header(
                show_spin=has_spin, summary_mode=True, long_labels=long_labels
            ),
            err=True,
        )
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
            top_n=top_n,
            absolute=absolute,
        )

        if not ranked:
            click.echo("No structures remaining after filtering.", err=True)
            return

        fill_dict_symm(ranked, symprec=symprec)
        click.echo(format_header(show_spin=has_spin, long_labels=long_labels), err=True)
        for rec in ranked:
            line = format_rank_line(rec, long_labels=long_labels, show_spin=has_spin)
            click.echo(line)

"""CLI command for converting between RES and extxyz formats."""

import sys
from pathlib import Path

import click


@click.command("convert")
@click.argument("input_path", type=click.Path(exists=True))
@click.argument("output_path", type=click.Path())
@click.option(
    "-l",
    "--label",
    default=None,
    help="Extract a single structure by label instead of bulk conversion",
)
def convert(input_path, output_path, label):
    """Convert between packed RES and extxyz formats.

    \b
    Bulk conversion (all structures):
      ap convert packed.res output.xyz       # .res → extxyz
      ap convert input.xyz output_dir/       # extxyz → unpacked .res files

    \b
    Extract a single structure by label:
      ap convert -l Si-002 packed.res Si-002.res
      ap convert -l Si-001 input.xyz Si-001.res
    """
    from airsspy.convert import extract_structure, extxyz_to_res, res_to_extxyz

    inp = Path(input_path)
    out = Path(output_path)

    if label is not None:
        # Extract single structure
        found = extract_structure(inp, label, out)
        if not found:
            click.echo(f"Structure '{label}' not found in {inp}", err=True)
            sys.exit(1)
        click.echo(f"Extracted '{label}' → {out}", err=True)
        return

    # Bulk conversion
    in_fmt = "extxyz" if inp.suffix in (".xyz", ".extxyz") else "res"

    if in_fmt == "res":
        if out.suffix in (".xyz", ".extxyz"):
            # packed .res → extxyz
            n = res_to_extxyz(inp, out)
            click.echo(f"Converted {n} structures: {inp} → {out}", err=True)
        else:
            # packed .res → unpacked .res files (just split)
            out.mkdir(parents=True, exist_ok=True)
            from airsspy.restools import RESFile

            # Reload with structures for full line preservation
            res_full = RESFile.from_packed(str(inp), include_structure=True)
            for res in res_full:
                lbl = res.label or "Unknown"
                safe = lbl.replace("/", "_").replace(" ", "_")
                raw = res.lines
                lines = [ln.rstrip("\n") for ln in raw] if raw else res.to_res_lines()
                if lines and not lines[-1].strip().startswith("END"):
                    lines.append("END")
                (out / f"{safe}.res").write_text("\n".join(lines) + "\n")
            click.echo(f"Unpacked {len(res_full)} structures to {out}/", err=True)
    else:
        # extxyz → unpacked .res files
        out.mkdir(parents=True, exist_ok=True)
        n = extxyz_to_res(inp, out)
        click.echo(f"Converted {n} structures: {inp} → {out}/", err=True)

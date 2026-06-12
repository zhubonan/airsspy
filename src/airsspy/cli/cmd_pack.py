"""CLI commands for packing and unpacking structure files."""

from pathlib import Path

import click


@click.command("pack")
@click.argument("inputs", nargs=-1, type=click.Path(exists=True))
@click.argument("output", type=click.Path())
@click.option(
    "--from-dir",
    default=None,
    type=click.Path(exists=True),
    help="Pack all .res or .xyz files from this directory",
)
def pack(inputs, output, from_dir):
    """Concatenate individual structure files into one packed file.

    \b
      ap pack *.res output.res
      ap pack *.xyz output.xyz
      ap pack --from-dir res_files/ output.res
      ap pack --from-dir xyz_files/ output.xyz
    """
    out = Path(output)

    files: list[Path] = []
    if from_dir:
        d = Path(from_dir)
        ext = out.suffix
        if ext not in (".res", ".xyz", ".extxyz"):
            click.echo("Error: output must be .res, .xyz, or .extxyz", err=True)
            raise SystemExit(1)
        if ext in (".xyz", ".extxyz"):
            patterns = [".xyz", ".extxyz"]
        else:
            patterns = [ext]
        for p in patterns:
            files.extend(sorted(d.glob(f"*{p}")))
    else:
        files = [Path(f) for f in inputs]

    if not files:
        click.echo("No input files found.", err=True)
        return

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fout:
        for fpath in files:
            content = fpath.read_text()
            fout.write(content)
            if not content.endswith("\n"):
                fout.write("\n")

    click.echo(f"Packed {len(files)} files → {out}", err=True)


@click.command("unpack")
@click.argument("input", type=click.Path(exists=True))
@click.argument("output_dir", type=click.Path())
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["res", "extxyz"]),
    default=None,
    help="Output format (default: same as input)",
)
def unpack(input, output_dir, fmt):
    """Split a packed file into individual structure files.

    \b
      ap unpack packed.res output_dir/
      ap unpack packed.xyz output_dir/

    Use ``ap convert`` to change between RES and extxyz formats.
    """
    inp = Path(input)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    if fmt is None:
        if inp.suffix == ".res":
            fmt = "res"
        else:
            fmt = "extxyz"

    if fmt == "res":
        _unpack_res(inp, out)
    else:
        out_ext = inp.suffix if inp.suffix in (".xyz", ".extxyz") else ".xyz"
        _unpack_extxyz(inp, out, out_ext)


def _unpack_res(inp: Path, out: Path) -> None:
    """Split packed .res on TITL...END blocks."""
    text = inp.read_text()
    lines = text.splitlines(True)

    current: list[str] = []
    label = None
    count = 0

    def flush():
        nonlocal label, count
        if not current or label is None:
            return
        safe = label.replace("/", "_").replace(" ", "_")
        (out / f"{safe}.res").write_text("".join(current))
        count += 1
        current.clear()
        label = None

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("TITL"):
            if current:
                flush()
            tokens = stripped.split()
            label = tokens[1] if len(tokens) > 1 else f"struct_{count:04d}"
        current.append(line)
        if stripped == "END":
            flush()

    if current:
        flush()

    click.echo(f"Unpacked {count} structures → {out}/", err=True)


def _unpack_extxyz(inp: Path, out: Path, out_ext: str = ".xyz") -> None:
    """Split packed extxyz on structure boundaries.

    Each structure starts with a line containing an integer (atom count)
    followed by a properties line, then that many atom lines.
    """
    text = inp.read_text()
    lines = text.splitlines()

    structures: list[list[str]] = []
    current: list[str] = []
    i = 0

    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        try:
            natoms = int(line.split()[0])
            if natoms > 0 and len(line.split()) == 1:
                if current:
                    structures.append(current)
                    current = []
                current.append(lines[i])
                i += 1
                if i < len(lines):
                    current.append(lines[i])
                    i += 1
                    for _ in range(natoms):
                        if i < len(lines):
                            current.append(lines[i])
                            i += 1
                continue
        except (ValueError, IndexError):
            pass

        current.append(lines[i])
        i += 1

    if current:
        structures.append(current)

    count = 0
    for struct_lines in structures:
        label = _extract_label_from_extxyz(struct_lines, count)
        safe = label.replace("/", "_").replace(" ", "_")
        (out / f"{safe}{out_ext}").write_text("\n".join(struct_lines) + "\n")
        count += 1

    click.echo(f"Unpacked {count} structures → {out}/", err=True)


def _extract_label_from_extxyz(lines: list[str], fallback_idx: int) -> str:
    """Extract label from extxyz info line (second line of structure)."""
    if len(lines) >= 2:
        info_line = lines[1]
        for part in info_line.split():
            if part.startswith("label="):
                return part.split("=", 1)[1].strip("'\"")
            if part.startswith("name="):
                return part.split("=", 1)[1].strip("'\"")
            if part.startswith("structure_id="):
                return part.split("=", 1)[1].strip("'\"")
    return f"struct_{fallback_idx:04d}"

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
            ends_with_newline = True
            with open(fpath) as fin:
                while True:
                    chunk = fin.read(1024 * 1024)
                    if not chunk:
                        break
                    ends_with_newline = chunk.endswith("\n")
                    fout.write(chunk)
            if not ends_with_newline:
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
    from airsspy.restools import iter_res_blocks

    count = 0
    with open(inp) as stream:
        for current in iter_res_blocks(stream):
            label = None
            for line in current:
                stripped = line.strip()
                if stripped.startswith("TITL"):
                    tokens = stripped.split()
                    label = tokens[1] if len(tokens) > 1 else None
                    break
            if label is None:
                label = f"struct_{count:04d}"
            safe = label.replace("/", "_").replace(" ", "_")
            if current and not current[-1].strip().startswith("END"):
                current = [*current, "END\n"]
            (out / f"{safe}.res").write_text("".join(current))
            count += 1

    click.echo(f"Unpacked {count} structures → {out}/", err=True)


def _unpack_extxyz(inp: Path, out: Path, out_ext: str = ".xyz") -> None:
    """Split packed extxyz on structure boundaries.

    Each structure starts with a line containing an integer (atom count)
    followed by a properties line, then that many atom lines.
    """
    current: list[str] = []
    count = 0
    with open(inp) as stream:
        lines_iter = iter(stream)
        for first in lines_iter:
            line = first.strip()
            if not line:
                continue
            try:
                natoms = int(line.split()[0])
            except (ValueError, IndexError):
                current.append(first.rstrip("\n"))
                continue
            if natoms <= 0 or len(line.split()) != 1:
                current.append(first.rstrip("\n"))
                continue
            current = [first.rstrip("\n")]
            try:
                current.append(next(lines_iter).rstrip("\n"))
            except StopIteration:
                break
            for _ in range(natoms):
                try:
                    current.append(next(lines_iter).rstrip("\n"))
                except StopIteration:
                    break
            label = _extract_label_from_extxyz(current, count)
            safe = label.replace("/", "_").replace(" ", "_")
            (out / f"{safe}{out_ext}").write_text("\n".join(current) + "\n")
            count += 1
            continue

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

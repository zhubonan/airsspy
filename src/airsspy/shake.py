"""Helpers for shaking candidate structures before re-relaxation."""

from __future__ import annotations

import csv
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read

from airsspy.restools import RESFile, _get_res_lines


@dataclass(frozen=True)
class ShakeRecord:
    """Metadata for one shaken structure written to disk."""

    source_path: Path
    source_label: str
    shake_label: str
    shake_index: int
    random_seed: int
    input_path: Path
    source_enthalpy: float | None
    source_enthalpy_per_atom: float | None


def collect_shake_input_paths(
    files: Sequence[str | Path],
    *,
    filelist: str | Path | None = None,
    base_dir: str | Path | None = None,
) -> list[Path]:
    """Collect positional and filelist inputs, preserving order and uniqueness."""
    root = Path(base_dir or ".").resolve()
    entries: list[tuple[str | Path, Path]] = [(entry, root) for entry in files]
    if filelist is not None:
        filelist_path = _resolve_existing_path(filelist, root)
        for line in filelist_path.read_text().splitlines():
            item = line.strip()
            if not item or item.startswith("#"):
                continue
            entries.append((item, filelist_path.parent))

    paths: list[Path] = []
    seen: set[Path] = set()
    for entry, entry_root in entries:
        path = _resolve_existing_path(entry, entry_root)
        if path in seen:
            continue
        seen.add(path)
        paths.append(path)
    return paths


def prepare_shaken_inputs(
    input_paths: Sequence[Path],
    output_dir: str | Path,
    *,
    rattles: int = 1,
    rattle_sigma: float = 0.05,
    strain_sigma: float = 0.0,
    seed: int = 1,
    clean: bool = False,
) -> list[ShakeRecord]:
    """Write shaken RES inputs and return their manifest records."""
    if rattles < 1:
        raise ValueError("rattles must be at least 1")
    if rattle_sigma < 0:
        raise ValueError("rattle_sigma must be non-negative")
    if strain_sigma < 0:
        raise ValueError("strain_sigma must be non-negative")

    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if clean:
        for old_file in output.glob("*.res"):
            old_file.unlink()

    records: list[ShakeRecord] = []
    used_labels: set[str] = set()
    for input_number, source_path in enumerate(input_paths):
        atoms, source_label, enthalpy, spins = read_shake_source(source_path)
        enthalpy_per_atom = enthalpy / len(atoms) if enthalpy is not None else None
        base_label = _sanitize_label(source_label or source_path.stem)

        for shake_index in range(1, rattles + 1):
            random_seed = seed + input_number * rattles + shake_index - 1
            shake_label = _unique_label(
                f"{base_label}-shake{shake_index:02d}", used_labels
            )
            shaken = shake_atoms(
                atoms,
                rattle_sigma=rattle_sigma,
                strain_sigma=strain_sigma,
                seed=random_seed,
            )
            target = output / f"{shake_label}.res"
            write_shaken_res(
                shaken,
                target,
                label=shake_label,
                source_path=source_path,
                source_label=source_label,
                random_seed=random_seed,
                rattle_sigma=rattle_sigma,
                strain_sigma=strain_sigma,
                spins=spins,
            )
            records.append(
                ShakeRecord(
                    source_path=source_path,
                    source_label=source_label,
                    shake_label=shake_label,
                    shake_index=shake_index,
                    random_seed=random_seed,
                    input_path=target,
                    source_enthalpy=enthalpy,
                    source_enthalpy_per_atom=enthalpy_per_atom,
                )
            )
    return records


def read_shake_source(
    path: str | Path,
) -> tuple[Atoms, str, float | None, list[float]]:
    """Read a structure that can be shaken."""
    source = Path(path).resolve()
    if source.suffix.lower() == ".res":
        res = RESFile.from_file(str(source))
        atoms = res.atoms
        if atoms is None:
            raise ValueError(f"RES file has no readable structure: {source}")
        enthalpy = float(res.enthalpy) if res.enthalpy is not None else None
        return atoms, str(res.label or source.stem), enthalpy, list(res.spins)

    atoms = read(str(source))
    label = str(atoms.info.get("label") or source.stem)
    energy = atoms.info.get("enthalpy", atoms.info.get("energy"))
    return atoms, label, float(energy) if energy is not None else None, []


def shake_atoms(
    atoms: Atoms,
    *,
    rattle_sigma: float = 0.05,
    strain_sigma: float = 0.0,
    seed: int = 1,
) -> Atoms:
    """Return a copy of *atoms* with random positional and optional cell noise."""
    shaken = atoms.copy()
    if strain_sigma:
        rng = np.random.default_rng(seed)
        strain = rng.normal(scale=strain_sigma, size=(3, 3))
        strain = 0.5 * (strain + strain.T)
        shaken.set_cell((np.eye(3) + strain) @ shaken.cell.array, scale_atoms=True)
    if rattle_sigma:
        shaken.rattle(stdev=rattle_sigma, seed=seed)
    shaken.wrap()
    return shaken


def write_shaken_res(
    atoms: Atoms,
    path: str | Path,
    *,
    label: str,
    source_path: str | Path,
    source_label: str,
    random_seed: int,
    rattle_sigma: float,
    strain_sigma: float,
    spins: list[float] | None = None,
) -> None:
    """Write a shaken structure as a single AIRSS RES file."""
    species = atoms.get_chemical_symbols()
    scaled_positions = atoms.get_scaled_positions(wrap=True).tolist()
    cellpar = atoms.cell.cellpar().tolist()
    titl = [
        label,
        0.0,
        atoms.get_volume(),
        0.0,
        0.0,
        0.0,
        len(atoms),
        "P1",
        "n",
        "-",
        "1",
    ]
    rems = [
        f"SHAKE_SOURCE {Path(source_path).resolve()}",
        f"SHAKE_SOURCE_LABEL {source_label}",
        f"SHAKE_SEED {random_seed}",
        f"SHAKE_RATTLE_SIGMA {rattle_sigma}",
        f"SHAKE_STRAIN_SIGMA {strain_sigma}",
    ]
    lines = _get_res_lines(
        titl,
        species,
        scaled_positions,
        cellpar,
        rem_lines=rems,
        spins=spins,
    )
    Path(path).write_text("\n".join(lines) + "\n")


def write_shake_manifest(records: Sequence[ShakeRecord], path: str | Path) -> None:
    """Write a CSV manifest for prepared shaken inputs."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "source_path",
                "source_label",
                "shake_label",
                "shake_index",
                "random_seed",
                "input_path",
                "source_enthalpy",
                "source_enthalpy_per_atom",
            ],
        )
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "source_path": str(record.source_path),
                    "source_label": record.source_label,
                    "shake_label": record.shake_label,
                    "shake_index": record.shake_index,
                    "random_seed": record.random_seed,
                    "input_path": str(record.input_path),
                    "source_enthalpy": _optional_float(record.source_enthalpy),
                    "source_enthalpy_per_atom": _optional_float(
                        record.source_enthalpy_per_atom
                    ),
                }
            )


def _resolve_existing_path(entry: str | Path, root: Path) -> Path:
    path = Path(entry).expanduser()
    if not path.is_absolute():
        path = root / path
    path = path.resolve()
    if not path.exists():
        raise FileNotFoundError(str(path))
    if not path.is_file():
        raise IsADirectoryError(str(path))
    return path


def _sanitize_label(label: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.+-]+", "-", label.strip())
    return safe.strip("-") or "shake"


def _unique_label(label: str, used: set[str]) -> str:
    if label not in used:
        used.add(label)
        return label
    counter = 2
    while f"{label}-{counter}" in used:
        counter += 1
    unique = f"{label}-{counter}"
    used.add(unique)
    return unique


def _optional_float(value: float | None) -> str:
    return "" if value is None else f"{value:.12g}"

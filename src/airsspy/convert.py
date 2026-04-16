"""Lossless conversion between AIRSS .res and extended XYZ formats.

Supports round-tripping all data including forces (stored as extra
columns 8-10 on atom lines in .res), per-atom spins, REM metadata,
and all TITL fields. Output .res files are fully compatible with
cryan and other AIRSS tools.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.geometry import cell_to_cellpar
from ase.io import read as ase_read
from ase.io import write as ase_write

from .restools import RESFile, _get_res_lines, parse_titl

# ---------------------------------------------------------------------------
# Force parsing helpers
# ---------------------------------------------------------------------------


def _parse_res_forces(lines: list[str]) -> list[list[float]] | None:
    """Parse force columns from RES atom lines.

    Atom line format: Symbol index x y z occ [spin] [fx fy fz]
    Base columns (6): Symbol(0) index(1) x(2) y(3) z(4) occ(5)
    Optional spin(6), then forces(7,8,9) — or if no spin, forces(6,7,8).

    We detect forces by checking if there are 3+ extra columns beyond base
    (for no-spin case: 9 tokens → forces at 6,7,8) or 4+ extra (with spin:
    10 tokens → forces at 7,8,9).

    Returns None if no force columns are present.
    """
    # First pass: determine if spin column is present
    has_spin = False
    for line in lines:
        tokens = line.split()
        if not tokens or tokens[0] != "SFAC":
            continue
        # Look at the next atom line to determine column count
        break

    for line in lines:
        tokens = line.split()
        if not tokens:
            continue
        if tokens[0] == "SFAC":
            # Check first atom line after SFAC for column count
            continue
        if tokens[0] in ("TITL", "CELL", "LATT", "REM", "END"):
            continue
        if tokens[0] and tokens[0][0].isalpha() and len(tokens) > 6:
            # Check if column 6 looks like a spin value (small number)
            # vs a force value. If we have 10 tokens, col 6 is spin, 7-9 are forces.
            # If 9 tokens, col 6 is first force component.
            if len(tokens) >= 10:
                has_spin = True
            break

    # Second pass: parse forces
    forces: list[list[float]] = []
    in_sfac = False
    found_forces = False
    force_start = 7 if has_spin else 6  # index of first force component

    for line in lines:
        tokens = line.split()
        if not tokens:
            continue
        if tokens[0] == "SFAC":
            in_sfac = True
            continue
        if tokens[0] == "END":
            in_sfac = False
            continue
        if tokens[0] in ("TITL", "CELL", "LATT", "REM"):
            continue
        if in_sfac and tokens[0] and tokens[0][0].isalpha():
            n_extra = len(tokens) - 6  # columns beyond base (symbol..occ)
            # With spin: n_extra=1 (spin only), 4 (spin+forces)
            # Without spin: n_extra=0 (base only), 3 (forces only)
            expected_for_force = 4 if has_spin else 3
            if n_extra >= expected_for_force:
                try:
                    fx = float(tokens[force_start])
                    fy = float(tokens[force_start + 1])
                    fz = float(tokens[force_start + 2])
                    forces.append([fx, fy, fz])
                    found_forces = True
                except (ValueError, IndexError):
                    forces.append([0.0, 0.0, 0.0])
            else:
                forces.append([0.0, 0.0, 0.0])

    return forces if found_forces else None


# ---------------------------------------------------------------------------
# RES → extxyz
# ---------------------------------------------------------------------------


def _resfile_to_atoms(res: RESFile) -> Atoms:
    """Convert a RESFile to an ASE Atoms object with full metadata."""
    if res.structure is None:
        raise ValueError(f"RESFile '{res.label}' has no structure loaded")

    atoms = res.atoms
    if atoms is None:
        raise ValueError(f"RESFile '{res.label}' could not produce Atoms")

    # Structure-level info
    atoms.info["label"] = res.label or "Unknown"
    atoms.info["pressure"] = res.pressure
    atoms.info["spin"] = res.spin
    atoms.info["spin_abs"] = res.spin_abs
    atoms.info["symm"] = res.symm or "-"

    # Copies from raw TITL line
    if res.lines:
        for line in res.lines:
            if "TITL" in line:
                try:
                    ti = parse_titl(line)
                    atoms.info["copies"] = int(ti.flag3)
                except (ValueError, IndexError):
                    atoms.info["copies"] = 1
                break
    else:
        atoms.info["copies"] = 1

    # REM lines
    rem = res.rem
    if rem:
        atoms.info["rem"] = rem

    # Per-atom spins
    spins = res.spins
    if spins:
        atoms.set_initial_magnetic_moments(spins)

    # Forces from raw lines
    if res.lines:
        forces = _parse_res_forces(res.lines)
        if forces is not None:
            forces_arr = np.array(forces)
            calc = SinglePointCalculator(atoms, energy=res.enthalpy, forces=forces_arr)
            atoms.calc = calc
            return atoms

    # Energy only (no forces)
    if res.enthalpy is not None:
        calc = SinglePointCalculator(atoms, energy=res.enthalpy)
        atoms.calc = calc

    return atoms


def res_to_extxyz(res_path: str | Path, extxyz_path: str | Path) -> int:
    """Convert a (packed) .res file to extxyz.

    Returns the number of structures converted.
    """
    res_path = Path(res_path)
    extxyz_path = Path(extxyz_path)

    res_objs = RESFile.from_packed(str(res_path), include_structure=True)
    atoms_list = [_resfile_to_atoms(r) for r in res_objs]

    ase_write(str(extxyz_path), atoms_list, format="extxyz")
    return len(atoms_list)


# ---------------------------------------------------------------------------
# extxyz → RES
# ---------------------------------------------------------------------------


def _extract_energy(atoms: Atoms) -> float:
    """Extract energy from calculator or info dict."""
    if atoms.calc is not None:
        try:
            return float(atoms.get_potential_energy())
        except Exception:
            pass
    return float(atoms.info.get("energy", 0.0))


def _extract_forces(atoms: Atoms) -> np.ndarray | None:
    """Extract forces from calculator or arrays."""
    if atoms.calc is not None:
        try:
            return np.array(atoms.get_forces())
        except Exception:
            pass
    if "forces" in atoms.arrays:
        return np.array(atoms.arrays["forces"])
    return None


def _atoms_to_res_lines(atoms: Atoms) -> list[str]:
    """Convert an ASE Atoms object to RES lines (with optional force columns)."""
    natoms = len(atoms)
    species = atoms.get_chemical_symbols()
    scaled = atoms.get_scaled_positions()
    cellpar = cell_to_cellpar(atoms.cell)

    # Extract metadata
    label = str(atoms.info.get("label", atoms.info.get("name", "Unknown")))
    pressure = float(atoms.info.get("pressure", 0.0))
    energy = _extract_energy(atoms)
    volume = atoms.get_volume()
    spin = float(atoms.info.get("spin", 0.0))
    spin_abs = float(atoms.info.get("spin_abs", 0.0))
    symm = str(atoms.info.get("symm", atoms.info.get("spacegroup", "-")))
    copies = int(atoms.info.get("copies", 1))

    # Per-atom spins
    spins = None
    if atoms.has("initial_magmoms"):
        magmoms = atoms.get_initial_magnetic_moments()
        if np.any(magmoms != 0):
            spins = [float(m) for m in magmoms]

    # REM lines
    rem_lines = atoms.info.get("rem")
    if isinstance(rem_lines, str):
        rem_lines = [rem_lines]

    # Build TITL list matching _get_res_lines format
    # Format: "{} {:.3f} {:.3f} {:.4f} {:.2f} {:.2f} {} ({}) {} {} {}"
    # symm goes into "({})" so pass bare name
    symm_bare = symm.strip("()")
    titl_list = [
        label,
        pressure,
        volume,
        energy,
        spin,
        spin_abs,
        natoms,
        symm_bare,
        "n",
        "-",
        str(copies),
    ]

    # Get base RES lines (without forces)
    lines = _get_res_lines(
        titl_list, species, scaled.tolist(), cellpar, rem_lines, spins
    )

    # Append force columns if present
    forces = _extract_forces(atoms)
    if forces is not None:
        # Find atom lines (between SFAC and END) and append forces
        new_lines: list[str] = []
        atom_idx = 0
        in_atoms = False
        for line in lines:
            if line.startswith("SFAC"):
                in_atoms = True
                new_lines.append(line)
                continue
            if line.strip() == "END":
                in_atoms = False
                new_lines.append(line)
                continue
            if in_atoms and line and line[0].isalpha():
                # Atom line — append forces
                if atom_idx < len(forces):
                    fx, fy, fz = forces[atom_idx]
                    line = f"{line} {fx:>12.6f} {fy:>12.6f} {fz:>12.6f}"
                    atom_idx += 1
            new_lines.append(line)
        lines = new_lines

    return lines


def extxyz_to_res(extxyz_path: str | Path, output_dir: str | Path) -> int:
    """Convert an extxyz file to individual .res files (one per structure).

    Each structure is written to ``<output_dir>/<label>.res``.
    If output_dir doesn't exist, it is created.

    Returns the number of structures converted.
    """
    extxyz_path = Path(extxyz_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        atoms_list = ase_read(str(extxyz_path), index=":")
    except Exception:
        atoms_list = [ase_read(str(extxyz_path))]

    for i, atoms in enumerate(atoms_list):
        res_lines = _atoms_to_res_lines(atoms)
        label = str(atoms.info.get("label", atoms.info.get("name", f"struct_{i:04d}")))
        # Sanitize label for filename
        safe_label = label.replace("/", "_").replace(" ", "_")
        res_file = output_dir / f"{safe_label}.res"
        res_file.write_text("\n".join(res_lines) + "\n")

    return len(atoms_list)


# ---------------------------------------------------------------------------
# Extract a single structure by label
# ---------------------------------------------------------------------------


def extract_structure(
    source_path: str | Path,
    label: str,
    output_path: str | Path,
    input_format: str | None = None,
) -> bool:
    """Extract a single structure by label from a packed .res or extxyz file.

    Args:
        source_path: Path to packed .res or extxyz file.
        label: Structure label to match (from TITL field 1 or atoms.info['label']).
        output_path: Path to write the extracted structure.
            The format is determined by the file extension (.res or .xyz/.extxyz).
        input_format: Force input format ('res', 'extxyz', or None for auto-detect).

    Returns:
        True if the structure was found and written, False otherwise.
    """
    source_path = Path(source_path)
    output_path = Path(output_path)

    fmt = input_format
    if fmt is None:
        if source_path.suffix in (".xyz", ".extxyz"):
            fmt = "extxyz"
        else:
            fmt = "res"

    out_fmt = None
    if output_path.suffix in (".xyz", ".extxyz"):
        out_fmt = "extxyz"
    else:
        out_fmt = "res"

    # Find the matching structure
    if fmt == "extxyz":
        atoms = _find_in_extxyz(source_path, label)
        if atoms is None:
            return False
        if out_fmt == "res":
            res_lines = _atoms_to_res_lines(atoms)
            output_path.write_text("\n".join(res_lines) + "\n")
        else:
            ase_write(str(output_path), atoms, format="extxyz")
    else:
        res = _find_in_packed_res(source_path, label)
        if res is None:
            return False
        if out_fmt == "extxyz":
            atoms = _resfile_to_atoms(res)
            ase_write(str(output_path), atoms, format="extxyz")
        else:
            # RES → RES: write raw lines (from_packed strips END, so add it back)
            if res.lines:
                lines = [ln.rstrip("\n") for ln in res.lines]
                if not lines[-1].strip().startswith("END"):
                    lines.append("END")
                output_path.write_text("\n".join(lines) + "\n")
            else:
                res_lines = res.to_res_lines()
                output_path.write_text("\n".join(res_lines))

    return True


def _find_in_extxyz(source_path: Path, label: str) -> Atoms | None:
    """Find a structure by label in an extxyz file."""
    try:
        atoms_list = ase_read(str(source_path), index=":")
    except Exception:
        atoms_list = [ase_read(str(source_path))]

    for atoms in atoms_list:
        al = str(atoms.info.get("label", atoms.info.get("name", "")))
        if al == label:
            return atoms
    return None


def _find_in_packed_res(source_path: Path, label: str) -> RESFile | None:
    """Find a structure by label in a packed .res file.

    Uses fast TITL-only parsing first to avoid loading full structures
    for non-matching entries.
    """
    res_objs = RESFile.from_packed(
        str(source_path), include_structure=False, only_titl=True
    )

    target_idx = None
    for i, res in enumerate(res_objs):
        if res.label == label:
            target_idx = i
            break

    if target_idx is None:
        return None

    # Reload only the matching structure with full data
    return RESFile.from_packed(str(source_path), include_structure=True)[target_idx]

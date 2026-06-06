"""Structure ranking utilities for AIRSS search results.

Provides fast parsing of SHELX .res files, ranking by enthalpy per
formula unit, and optional merging of similar structures using
distance fingerprint comparison (equivalent to ``cryan -u``).
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field
from functools import reduce
from math import ceil, gcd, isfinite
from pathlib import Path
from typing import TextIO

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class StructureRecord:
    """Lightweight record for a ranked structure."""

    label: str
    pressure: float
    volume: float
    enthalpy: float
    spin: float = 0.0
    spin_abs: float = 0.0
    natoms: int = 0
    symm: str = ""
    species_counts: dict[str, int] = field(default_factory=dict)
    copies: int = 1
    source: str = ""
    # Raw lines kept for lazy full-structure loading (needed by eliminate_similar)
    _raw_lines: list[str] = field(default_factory=list, repr=False)
    _atoms: object | None = field(default=None, repr=False)
    _merged_peers: list = field(default_factory=list, repr=False)

    @property
    def reduced_formula(self) -> str:
        """Hill-system reduced formula, e.g. 'SiO2'."""
        return _reduce_formula(self.species_counts)

    @property
    def n_formula_units(self) -> int:
        """Number of formula units (GCD of species counts)."""
        counts = list(self.species_counts.values())
        if not counts:
            return 1
        if len(counts) == 1:
            return counts[0]
        return reduce(gcd, counts)

    @property
    def enthalpy_per_fu(self) -> float:
        """Enthalpy per formula unit."""
        nfu = self.n_formula_units
        return self.enthalpy / nfu if nfu > 0 else self.enthalpy

    @property
    def volume_per_fu(self) -> float:
        """Volume per formula unit."""
        nfu = self.n_formula_units
        return self.volume / nfu if nfu > 0 else self.volume


# ---------------------------------------------------------------------------
# Formula helpers
# ---------------------------------------------------------------------------


# Element symbols ordered by atomic number for formula generation.
# Matches cryan's ordering: C first, H second (if C present),
# others by atomic number, O always last.
_ELEMENT_ORDER: dict[str, int] = {}
_PERIODIC_TABLE = [
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]
for _i, _el in enumerate(_PERIODIC_TABLE):
    _ELEMENT_ORDER[_el] = _i + 1  # Z starts at 1


def _element_sort_key(el: str) -> float:
    """Sort key matching cryan: C first, H second (if C present), O last, rest by Z."""
    if el == "O":
        return float("inf")  # Oxygen always last
    if el == "C":
        return -2.0
    if el == "H":
        return -1.0
    return float(_ELEMENT_ORDER.get(el, 200))


def _reduce_formula(species_counts: dict[str, int]) -> str:
    """Compute reduced formula matching cryan's ordering convention.

    Ordering: C first, H second (only if C present), then by atomic
    number, O always last.
    """
    if not species_counts:
        return ""
    counts = list(species_counts.values())
    if len(counts) == 1:
        el = next(iter(species_counts))
        ct = counts[0]
        return el if ct == 1 else f"{el}{ct}"
    divisor = reduce(gcd, counts)

    elements = sorted(species_counts.keys(), key=_element_sort_key)

    parts: list[str] = []
    for el in elements:
        ct = species_counts[el] // divisor
        parts.append(el if ct == 1 else f"{el}{ct}")
    return "".join(parts)


def _parse_formula_counts(formula: str) -> dict[str, int] | None:
    """Parse a simple chemical formula into integer element counts."""
    import re

    if not formula:
        return None

    counts: dict[str, int] = {}
    pos = 0
    for match in re.finditer(r"([A-Z][a-z]?)(\d*)", formula):
        if match.start() != pos:
            return None
        element, count_text = match.groups()
        count = int(count_text or "1")
        if count <= 0:
            return None
        counts[element] = counts.get(element, 0) + count
        pos = match.end()

    if pos != len(formula):
        return None
    return counts


# ---------------------------------------------------------------------------
# Label truncation
# ---------------------------------------------------------------------------


def _truncate_label(label: str, width: int = 20) -> str:
    """Truncate *label* to *width* chars, appending ``...`` if truncated."""
    if len(label) <= width:
        return label
    return label[: width - 3] + "..."


# ---------------------------------------------------------------------------
# Extxyz field extraction helpers
# ---------------------------------------------------------------------------

_EV_A3_TO_GPA = 160.21766208


def apply_external_pressure(
    records: list[StructureRecord], pressure_gpa: float
) -> None:
    """Apply external pressure correction in-place.

    Adds ``P * V / _EV_A3_TO_GPA`` to each record's enthalpy, where
    *P* is in GPa and *V* in Å³.  Positive *pressure_gpa* favours
    denser (smaller-volume) structures.
    """
    for rec in records:
        pv_eV = pressure_gpa * rec.volume / _EV_A3_TO_GPA
        rec.enthalpy += pv_eV
        rec.pressure += pressure_gpa


def filter_by_name(
    records: list[StructureRecord], pattern: str
) -> list[StructureRecord]:
    """Filter records by label using a glob pattern.

    Supports ``*``, ``?``, and ``[seq]`` wildcards (fnmatch).
    """
    from fnmatch import fnmatch

    return [r for r in records if fnmatch(r.label, pattern)]


def filter_by_formula(
    records: list[StructureRecord], formula: str
) -> list[StructureRecord]:
    """Filter records by chemical formula.

    Three modes:

    * Exact reduced formula: ``-f SiO2``
    * Comma-separated elements: ``-f Si,O`` — matches any composition
      containing *all* listed elements
    * Glob on reduced formula: ``-f "Si*"`` — fnmatch on the reduced formula
    """
    if "," in formula:
        elements = {el.strip() for el in formula.split(",")}
        return [
            r for r in records if elements.issubset(r.species_counts.keys())
        ]

    from fnmatch import fnmatch

    if any(char in formula for char in "*?["):
        return [r for r in records if fnmatch(r.reduced_formula, formula)]

    parsed_counts = _parse_formula_counts(formula)
    if parsed_counts is not None:
        formula = _reduce_formula(parsed_counts)

    return [r for r in records if r.reduced_formula == formula]


def filter_by_formula_units(
    records: list[StructureRecord], n_formula_units: int
) -> list[StructureRecord]:
    """Filter records by exact number of formula units."""
    return [r for r in records if r.n_formula_units == n_formula_units]


def filter_by_species_number(
    records: list[StructureRecord], species_number: int
) -> list[StructureRecord]:
    """Filter records by exact number of distinct species."""
    return [r for r in records if len(r.species_counts) == species_number]


def filter_by_ions_number(
    records: list[StructureRecord], ions_number: int
) -> list[StructureRecord]:
    """Filter records by ion count.

    Positive values require an exact ``natoms`` match.  Negative values match
    records with ``natoms <= abs(ions_number)``, following cryan's range form.
    """
    if ions_number < 0:
        limit = abs(ions_number)
        return [r for r in records if r.natoms <= limit]
    return [r for r in records if r.natoms == ions_number]


def _extract_energy(atoms, field: str | None = None) -> float:
    """Extract energy from an ASE Atoms object.

    If *field* is given, read only ``atoms.info[field]``.
    Otherwise try common keys in order, then the attached calculator.
    Falls back to 0.0 with a warning.
    """
    if field is not None:
        val = atoms.info.get(field)
        if val is None:
            logger.warning("energy field %r not found in atoms.info", field)
            return 0.0
        return float(val)

    for key in ("energy", "enthalpy", "free_energy"):
        if key in atoms.info:
            return float(atoms.info[key])

    if atoms.calc is not None:
        try:
            return float(atoms.get_potential_energy())
        except Exception:
            pass

    logger.warning("no energy field found, defaulting to 0.0")
    return 0.0


def _extract_label(atoms, path: str, index: int, field: str | None = None) -> str:
    """Extract a structure label from an ASE Atoms object.

    If *field* is given, read only ``atoms.info[field]``.
    Otherwise try common keys in order.
    Falls back to ``"<basename>:<index>"``.
    """
    if field is not None:
        val = atoms.info.get(field)
        if val is not None:
            return str(val)
        logger.warning("label field %r not found in atoms.info", field)

    for key in ("label", "name", "structure_id", "source_label"):
        val = atoms.info.get(key)
        if val is not None and str(val).strip():
            return str(val)

    return f"{Path(path).name}:{index}"


def _stress_to_pressure_gpa(atoms) -> float:
    """Compute pressure (GPa) from the stress tensor on *atoms*.

    Pressure = -trace(stress) / 3, converted from eV/Å³ to GPa.
    Tries the calculator first, then ``atoms.info["stress"]``.
    """
    stress = None
    if atoms.calc is not None:
        try:
            stress = atoms.get_stress()
        except Exception:
            pass
    if stress is None:
        raw = atoms.info.get("stress")
        if raw is not None:
            stress = np.asarray(raw, dtype=float).ravel()

    if stress is None:
        return 0.0

    stress = np.asarray(stress, dtype=float).ravel()
    if stress.shape == (6,):
        trace = stress[0] + stress[1] + stress[2]
    elif stress.shape == (9,):
        trace = stress[0] + stress[4] + stress[8]
    elif stress.shape == (3, 3):
        trace = float(stress[0, 0] + stress[1, 1] + stress[2, 2])
    else:
        return 0.0

    return float(-trace / 3.0 * _EV_A3_TO_GPA)


def _extract_pressure(atoms, field: str | None = None) -> float:
    """Extract pressure (GPa) from an ASE Atoms object.

    If *field* is given, read only ``atoms.info[field]`` as a scalar.
    Otherwise try scalar keys in order, then the stress tensor.
    Falls back to 0.0.
    """
    if field is not None:
        val = atoms.info.get(field)
        if val is None:
            logger.warning("pressure field %r not found in atoms.info", field)
            return 0.0
        return float(val)

    for key in ("pressure", "extern_pressure"):
        if key in atoms.info:
            return float(atoms.info[key])

    return _stress_to_pressure_gpa(atoms)


def _extract_symm(atoms) -> str:
    """Extract spacegroup symbol from ``atoms.info`` only (no spglib)."""
    for key in ("symm", "spacegroup"):
        val = atoms.info.get(key)
        if val is not None and str(val).strip():
            return str(val)
    return ""


def fill_missing_spacegroups(
    records: list[StructureRecord], symprec: float = 0.01
) -> None:
    """Detect spacegroups via spglib for records with empty ``symm``.

    Only processes records that have an ``_atoms`` reference (extxyz).
    Modifies records in place.
    """
    for rec in records:
        if rec.symm:
            continue
        if rec._atoms is None:
            continue
        try:
            import spglib

            atoms = rec._atoms
            dataset = spglib.get_symmetry_dataset(
                (
                    atoms.cell,
                    atoms.get_scaled_positions(),
                    atoms.get_atomic_numbers(),
                ),
                symprec=symprec,
            )
            if dataset is not None and dataset.international:
                rec.symm = f"({dataset.international})"
        except Exception:
            pass


def fill_dict_symm(dicts: list[dict], symprec: float = 0.01) -> None:
    """Fill missing ``symm`` in ranking output dicts via spglib.

    Each dict must have a ``_record`` key referencing the source
    :class:`StructureRecord`.  Modifies dicts in place.
    """
    recs = [d["_record"] for d in dicts if not d.get("symm") and d.get("_record")]
    fill_missing_spacegroups(recs, symprec=symprec)
    for d in dicts:
        rec = d.get("_record")
        if rec is not None and not d.get("symm"):
            d["symm"] = rec.symm


# ---------------------------------------------------------------------------
# Fast RES parsing (no full structure construction)
# ---------------------------------------------------------------------------

_RES_KEYWORDS = frozenset({"TITL", "CELL", "LATT", "SFAC", "REM", "END", ""})


def _parse_res_fast(lines: list[str]) -> StructureRecord | None:
    """Parse a single RES structure from lines (TITL + species counts only).

    No pymatgen/ASE structure is constructed -- fast path for ranking.
    """
    label: str | None = None
    pressure = 0.0
    volume = 0.0
    enthalpy = 0.0
    spin = 0.0
    spin_abs = 0.0
    natoms = 0
    symm = ""
    copies = 1
    species_counts: Counter = Counter()
    in_sfac = False

    for line in lines:
        tokens = line.split()
        if not tokens:
            continue

        if tokens[0] == "TITL":
            # TITL label P V H spin spin_abs nat (symm) n - copies
            # But TITL may have variable-length fields
            # tokens[0]="TITL", [1]=label, [2]=P, [3]=V, [4]=H ...
            ntok = len(tokens) - 1  # exclude "TITL" itself
            if ntok >= 5:
                label = tokens[1]
                try:
                    pressure = float(tokens[2])
                except ValueError:
                    pressure = 0.0
                try:
                    volume = float(tokens[3])
                except ValueError:
                    volume = 0.0
                try:
                    enthalpy = float(tokens[4])
                except ValueError:
                    enthalpy = 0.0
            if ntok >= 7:
                try:
                    spin = float(tokens[5])
                    spin_abs = float(tokens[6])
                except ValueError:
                    pass
            if ntok >= 8:
                try:
                    natoms = int(tokens[7])
                except ValueError:
                    natoms = 0
            if ntok >= 9:
                symm = tokens[8]
            # Parse copies from "n - <N>" at end of line
            # Look for "n" followed by "-" followed by a number
            for i in range(len(tokens) - 2):
                if tokens[i] == "n" and tokens[i + 1] == "-":
                    try:
                        copies = int(tokens[i + 2])
                    except (ValueError, IndexError):
                        copies = 1
                    break

        elif tokens[0] == "SFAC":
            in_sfac = True

        elif tokens[0] == "END":
            in_sfac = False

        elif in_sfac and tokens[0] not in _RES_KEYWORDS:
            # Atom line: Symbol index x y z occ [spin]
            sp = tokens[0]
            if sp and sp[0].isalpha():
                species_counts[sp] += 1

        # Skip CELL, LATT, REM lines

    if label is None:
        return None

    actual_nat = sum(species_counts.values())
    if natoms == 0:
        natoms = actual_nat

    return StructureRecord(
        label=label,
        pressure=pressure,
        volume=volume,
        enthalpy=enthalpy,
        spin=spin,
        spin_abs=spin_abs,
        natoms=natoms,
        symm=symm,
        species_counts=dict(species_counts),
        copies=copies,
        _raw_lines=list(lines),
    )


# ---------------------------------------------------------------------------
# Input readers
# ---------------------------------------------------------------------------


def _read_res_from_lines_iter(line_iter) -> list[StructureRecord]:
    """Read concatenated RES structures from an iterable of lines."""
    records: list[StructureRecord] = []
    current: list[str] = []

    for line in line_iter:
        if isinstance(line, str):
            line = line.rstrip("\n")
        if line.startswith("END"):
            if current:
                rec = _parse_res_fast(current)
                if rec is not None:
                    records.append(rec)
            current = []
        else:
            current.append(line)

    # Handle trailing structure without END
    if current:
        rec = _parse_res_fast(current)
        if rec is not None:
            records.append(rec)

    return records


def read_res_stream(stream: TextIO) -> list[StructureRecord]:
    """Read concatenated RES structures from a text stream (stdin or file)."""
    records = _read_res_from_lines_iter(stream)
    for rec in records:
        rec.source = "stdin"
    return records


def read_res_file(path: str) -> list[StructureRecord]:
    """Read RES structures from a file (may be packed)."""
    with open(path) as fh:
        records = _read_res_from_lines_iter(fh)
    for rec in records:
        rec.source = path
    return records


def read_extxyz_file(
    path: str,
    energy_field: str | None = None,
    label_field: str | None = None,
    pressure_field: str | None = None,
) -> list[StructureRecord]:
    """Read structures from an extxyz file using ASE.

    Field names for energy, label and pressure are auto-detected from
    ``atoms.info`` / the attached calculator.  Override detection by
    passing *energy_field*, *label_field*, or *pressure_field*.

    Spacegroup is read from ``atoms.info`` if present; otherwise it
    is left empty and can be filled later via
    :func:`fill_missing_spacegroups`.
    """
    from ase.io import read as ase_read

    records: list[StructureRecord] = []
    try:
        atoms_list = ase_read(path, index=":")
    except Exception:
        atoms_list = [ase_read(path)]

    for i, atoms in enumerate(atoms_list):
        species_counts = dict(Counter(atoms.get_chemical_symbols()))
        energy = _extract_energy(atoms, field=energy_field)
        label = _extract_label(atoms, path, i, field=label_field)
        pressure = _extract_pressure(atoms, field=pressure_field)
        volume = atoms.get_volume()
        natoms = len(atoms)
        spin = float(atoms.info.get("spin", 0.0))
        spin_abs = float(atoms.info.get("spin_abs", 0.0))
        symm = _extract_symm(atoms)

        records.append(
            StructureRecord(
                label=label,
                pressure=pressure,
                volume=volume,
                enthalpy=energy,
                spin=spin,
                spin_abs=spin_abs,
                natoms=natoms,
                symm=symm,
                species_counts=species_counts,
                source=path,
                _atoms=atoms,
            )
        )
    return records


# ---------------------------------------------------------------------------
# Similarity / merging (cryan -u equivalent)
# ---------------------------------------------------------------------------


def _compute_distance_fingerprint(
    record: StructureRecord,
    cutoff: float = 4.0,
    zweight: bool = False,
) -> np.ndarray | None:
    """Compute a sorted distance fingerprint for a structure.

    Uses pymatgen's ``get_all_neighbors(cutoff)`` to find all distances
    to periodic images within *cutoff*, matching cryan's distance
    fingerprint algorithm.  When *zweight* is True, each distance *d* is
    weighted as ``d * (1 + log10(zmax² / (Z_i * Z_j)))`` to distinguish different
    atom-type pairs.

    Works for records loaded from either RES (via ``_raw_lines``) or
    extxyz (via ``_atoms``).  Returns None if the structure cannot be
    loaded.
    """
    structure = None

    if record._raw_lines:
        from .restools import RESFile

        try:
            res = RESFile.from_lines(record._raw_lines, include_structure=True)
            if res.structure is None:
                return None
            structure = res.structure
        except Exception:
            return None
    elif record._atoms is not None:
        try:
            from pymatgen.io.ase import AseAtomsAdaptor

            structure = AseAtomsAdaptor.get_structure(record._atoms)
        except Exception:
            return None
    else:
        return None

    neighbors = structure.get_all_neighbors(cutoff)

    if not zweight:
        all_dists = [
            float(n.nn_distance)
            for nlist in neighbors
            for n in nlist
        ]
    else:
        zmax = max(site.specie.Z for site in structure)
        all_dists = [
            float(n.nn_distance * (1.0 + np.log10(zmax * zmax / (structure[i].specie.Z * n.specie.Z))))
            for i, nlist in enumerate(neighbors)
            for n in nlist
        ]

    if not all_dists:
        return None

    return np.sort(np.array(all_dists, dtype=np.float64))


def eliminate_similar(
    records: list[StructureRecord],
    threshold: float,
    cutoff: float = 4.0,
    zweight: bool = False,
) -> list[StructureRecord]:
    """Merge similar structures by comparing distance fingerprints.

    Matches cryan's ``-u`` behaviour:
    1. Sort by enthalpy_per_fu ascending (most stable first)
    2. For each pair of same-formula records, compare scaled distance fingerprints
    3. If max difference < threshold * mean_min_distance, merge copies

    *cutoff* controls the neighbour search radius (Å) for fingerprint
    computation (default 4.0, matching cryan's ``rmax / 1.75``).

    When *zweight* is True, distances are weighted by ``d * zmax² / (Z_i·Z_j)``
    to distinguish different atom-type pairs.

    Merged peers are tracked in each surviving record's ``_merged_peers``
    list for later output.

    Returns the deduplicated list with accumulated copies.
    """
    from tqdm import tqdm

    groups: dict[str, list[StructureRecord]] = {}
    for rec in records:
        key = rec.reduced_formula
        groups.setdefault(key, []).append(rec)

    all_fps: dict[int, np.ndarray | None] = {}
    for rec in tqdm(records, desc="Computing fingerprints", unit="struct"):
        all_fps[id(rec)] = _compute_distance_fingerprint(
            rec, cutoff=cutoff, zweight=zweight
        )

    result: list[StructureRecord] = []
    total = len(records)

    with tqdm(total=total, desc="Comparing", unit="struct") as pbar:
        for _formula, group in groups.items():
            group.sort(key=lambda r: r.enthalpy_per_fu)

            n = len(group)
            fps = [all_fps[id(rec)] for rec in group]
            vpf = np.array([rec.volume_per_fu for rec in group])
            nfu = np.array([rec.n_formula_units for rec in group])

            merged_into = [-1] * n

            for i in range(n):
                pbar.update(1)
                if merged_into[i] >= 0:
                    continue
                if fps[i] is None:
                    continue

                fi = fps[i]
                nfi = len(fi)
                min_dist_i = fi[0] if nfi > 0 else 1.0
                vpf_i = vpf[i]
                nfi_rec = nfu[i]
                inv_vpf_i = 1.0 / vpf_i

                for j in range(i + 1, n):
                    if merged_into[j] >= 0:
                        continue
                    if fps[j] is None:
                        continue

                    fj = fps[j]
                    vpf_j = vpf[j]

                    if vpf_i <= 0 or vpf_j <= 0:
                        continue

                    rel_diff = abs(vpf_i - vpf_j) / max(vpf_i, vpf_j)
                    if rel_diff > 0.5:
                        continue

                    scale_a = ((vpf_j * inv_vpf_i + 1.0) * 0.5) ** (1.0 / 3.0)
                    scale_b = ((vpf_j + vpf_i) / (2.0 * vpf_j)) ** (1.0 / 3.0)

                    nfj_rec = nfu[j]
                    n_compare = min(nfi * nfj_rec, len(fj) * nfi_rec)

                    min_dist_j = fj[0] if len(fj) > 0 else 1.0
                    mean_min = (min_dist_i * scale_a + min_dist_j * scale_b) * 0.5
                    thresh = threshold * mean_min

                    ks = np.arange(n_compare)
                    iis = np.minimum(ks // nfj_rec, nfi - 1)
                    jjs = np.minimum(ks // nfi_rec, len(fj) - 1)

                    diffs = np.abs(fi[iis] * scale_a - fj[jjs] * scale_b)
                    nonzero = fi[iis] > 1e-10

                    if nonzero.any() and (diffs[nonzero] > thresh).any():
                        continue

                    merged_into[j] = i

            for i in range(n):
                if merged_into[i] >= 0:
                    target = merged_into[i]
                    while merged_into[target] >= 0:
                        target = merged_into[target]
                    group[target].copies += group[i].copies
                    group[target]._merged_peers.append(group[i])

            for i in range(n):
                if merged_into[i] < 0:
                    result.append(group[i])

    return result


# ---------------------------------------------------------------------------
# Phase diagram / Maxwell construction
# ---------------------------------------------------------------------------


def infer_elements(records: list[StructureRecord]) -> list[str]:
    """Infer the element list from all structures' species_counts.

    Returns elements sorted by atomic number.
    """
    all_elements: set[str] = set()
    for rec in records:
        all_elements.update(rec.species_counts.keys())
    return sorted(all_elements, key=lambda el: _ELEMENT_ORDER.get(el, 200))


def check_elemental_references(
    records: list[StructureRecord], elements: list[str]
) -> list[str]:
    """Check which elemental references are present.

    Returns a list of elements that have no pure-element structure.
    """
    # Build a set of compositions that are pure elements
    has_pure: set[str] = set()
    for rec in records:
        if len(rec.species_counts) == 1:
            el = next(iter(rec.species_counts))
            has_pure.add(el)
    return [el for el in elements if el not in has_pure]


def records_to_pd_entries(records: list[StructureRecord]) -> list:
    """Convert StructureRecords to pymatgen PDEntry objects.

    Uses Composition(species_counts) and total enthalpy as energy.
    Sets entry.name = label for display.
    """
    from pymatgen.analysis.phase_diagram import PDEntry
    from pymatgen.core.composition import Composition

    entries = []
    for rec in records:
        comp = Composition(rec.species_counts)
        entry = PDEntry(comp, energy=rec.enthalpy, name=rec.label)
        entries.append(entry)
    return entries


def maxwell_construction(
    records: list[StructureRecord],
    elements: list[str] | None = None,
    delta_e: float | None = None,
    verbose: bool = True,
) -> tuple[list[dict], object]:
    """Compute convex hull using pymatgen PhaseDiagram.

    Args:
        records: Structure records to analyse.
        elements: Element list for the chemical system. If None, inferred.
        delta_e: Filter structures with e_above_hull above this (eV/atom).
        verbose: If True, print warnings to stderr.

    Returns:
        Tuple of (output_records, PhaseDiagram, elements).
        Each output dict has all rank fields plus:
        - e_above_hull: energy above hull (eV/atom)
        - hull_energy_per_atom: hull energy at this composition (eV/atom)
        - formation_energy_per_atom: formation energy (eV/atom)
        - on_hull: True if structure is on the convex hull

    Raises:
        ValueError: If fewer than 2 elements in the system.
    """
    import warnings

    from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
    from pymatgen.core import Element
    from pymatgen.core.composition import Composition

    if elements is None:
        elements = infer_elements(records)

    if len(elements) < 2:
        raise ValueError(
            f"Need at least 2 elements for a phase diagram, got {len(elements)}: "
            f"{','.join(elements)}"
        )

    # Filter out records with no species data (e.g., unfinished runs)
    valid_records = [r for r in records if r.species_counts]
    if len(valid_records) < len(records):
        skipped = len(records) - len(valid_records)
        if verbose:
            logger.warning("skipping %d structures with no atom data", skipped)

    # Maxwell output is composition-level: use the lowest-enthalpy
    # representative for each reduced formula and accumulate copies.
    representative_records: list[StructureRecord] = []
    representative_copies: dict[str, int] = {}
    representative_by_formula: dict[str, StructureRecord] = {}
    representative_index: dict[str, int] = {}
    for rec in valid_records:
        formula = rec.reduced_formula
        representative_copies[formula] = (
            representative_copies.get(formula, 0) + rec.copies
        )

        if formula not in representative_by_formula:
            representative_by_formula[formula] = rec
            representative_index[formula] = len(representative_records)
            representative_records.append(rec)
            continue

        current = representative_by_formula[formula]
        if rec.enthalpy_per_fu < current.enthalpy_per_fu:
            representative_by_formula[formula] = rec
            representative_records[representative_index[formula]] = rec

    # Convert records to PDEntry
    entries = records_to_pd_entries(representative_records)

    # Add fake elemental references for any missing pure elements
    missing = check_elemental_references(representative_records, elements)
    fake_entries = []
    if missing:
        msg = f"Warning: no structures for pure elements: {', '.join(missing)}. Using E=0 references."
        if verbose:
            logger.warning(msg)
        for el in missing:
            fake_entry = PDEntry(
                Composition(el), energy=0.0, name=f"{el} (ref)"
            )
            fake_entries.append(fake_entry)
        entries.extend(fake_entries)

    # Build phase diagram — try with real Element objects first,
    # fall back to letting pymatgen auto-detect (handles dummy species like A, B)
    try:
        pd_elements = [Element(el) for el in elements]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pd = PhaseDiagram(entries, elements=pd_elements)
    except ValueError:
        # Non-standard species names (e.g., A, B) — let pymatgen auto-detect
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pd = PhaseDiagram(entries)

    # Compute hull data for each representative composition
    output_records: list[dict] = []
    for i, rec in enumerate(representative_records):
        entry = entries[i]
        comp = entry.composition
        formula = rec.reduced_formula

        # Energy per atom (enthalpy / natoms)
        h_per_atom = rec.enthalpy / rec.natoms if rec.natoms > 0 else rec.enthalpy

        # Hull energy and energy above hull
        hull_e = pd.get_hull_energy_per_atom(comp)
        e_above = pd.get_e_above_hull(entry)
        if e_above is None:
            e_above = float("inf")

        # Formation energy per atom
        form_e = pd.get_form_energy_per_atom(entry)

        on_hull = e_above < 1e-4

        # Filter by delta_e
        if delta_e is not None and e_above > delta_e:
            continue

        output_records.append(
            {
                "label": rec.label,
                "pressure": rec.pressure,
                "volume_per_fu": rec.volume_per_fu,
                "enthalpy_per_fu": rec.enthalpy_per_fu,
                "enthalpy_per_atom": h_per_atom,
                "hull_energy_per_atom": hull_e,
                "e_above_hull": e_above,
                "formation_energy_per_atom": form_e,
                "on_hull": on_hull,
                "spin_per_fu": rec.spin / rec.n_formula_units
                if rec.n_formula_units > 0
                else 0.0,
                "spin_abs_per_fu": rec.spin_abs / rec.n_formula_units
                if rec.n_formula_units > 0
                else 0.0,
                "nfu": rec.n_formula_units,
                "formula": formula,
                "symm": rec.symm,
                "copies": representative_copies[formula],
                "source": rec.source,
                "species_counts": rec.species_counts,
                "_record": rec,
            }
        )

    # Sort by e_above_hull (stable first), then by enthalpy_per_atom
    output_records.sort(key=lambda r: (r["e_above_hull"], r["enthalpy_per_atom"]))

    return output_records, pd, elements


# ---------------------------------------------------------------------------
# Ranking and output
# ---------------------------------------------------------------------------


def prefilter_records(
    records: list[StructureRecord],
    ethresh: float = 0.1,
) -> list[StructureRecord]:
    """Filter records by energy threshold per atom relative to the minimum.

    Groups by formula, computes relative enthalpy per atom within each
    group, removes structures above *ethresh* eV/atom.  Used to reduce
    the candidate set before merging.  Returns surviving
    ``StructureRecord`` objects.
    """
    groups: dict[str, list[StructureRecord]] = {}
    for rec in records:
        groups.setdefault(rec.reduced_formula, []).append(rec)

    surviving: list[StructureRecord] = []
    for group in groups.values():
        min_h_per_fu = min(r.enthalpy_per_fu for r in group)

        for rec in group:
            rel_h = rec.enthalpy_per_fu - min_h_per_fu
            rel_h_per_atom = (
                rel_h * rec.n_formula_units / rec.natoms
                if rec.natoms > 0
                else 0.0
            )
            if rel_h_per_atom > ethresh:
                continue
            surviving.append(rec)

    surviving.sort(key=lambda r: r.enthalpy_per_fu)

    return surviving


def _enthalpy_per_atom(rec: StructureRecord) -> float:
    return rec.enthalpy / rec.natoms if rec.natoms > 0 else rec.enthalpy


def prune_pathological_records(
    records: list[StructureRecord],
    tail_fraction: float = 0.10,
    sigma_factor: float = 3.0,
    trim_count: int = 1,
    min_tail_size: int = 5,
) -> tuple[list[StructureRecord], list[StructureRecord], list[dict]]:
    """Remove suspiciously low-energy records using a trimmed MAD cutoff.

    The filter is applied independently for each reduced formula.  Energies are
    compared as enthalpy per atom.  For each formula group, the lowest
    ``tail_fraction`` of records is used as the candidate tail, the lowest
    ``trim_count`` of those records are excluded from the baseline statistics,
    and the cutoff is ``median - sigma_factor * 1.4826 * MAD``.

    Returns ``(kept, rejected, diagnostics)``.  Diagnostics are dictionaries so
    callers can report skipped groups and per-formula cutoffs without redoing
    the statistics.
    """
    if not isfinite(tail_fraction) or not 0.0 < tail_fraction <= 1.0:
        raise ValueError("tail_fraction must be in (0, 1]")
    if not isfinite(sigma_factor) or sigma_factor < 0.0:
        raise ValueError("sigma_factor must be finite and >= 0")
    if trim_count < 0:
        raise ValueError("trim_count must be >= 0")
    if min_tail_size < 1:
        raise ValueError("min_tail_size must be >= 1")

    groups: dict[str, list[StructureRecord]] = {}
    for rec in records:
        groups.setdefault(rec.reduced_formula, []).append(rec)

    rejected_ids: set[int] = set()
    diagnostics: list[dict] = []

    for formula, group in groups.items():
        ranked = sorted(group, key=_enthalpy_per_atom)
        tail_count = min(len(ranked), max(1, ceil(tail_fraction * len(ranked))))
        tail = ranked[:tail_count]
        baseline = tail[min(trim_count, len(tail)) :]

        diagnostic = {
            "formula": formula,
            "group_size": len(group),
            "tail_size": len(tail),
            "trim_count": min(trim_count, len(tail)),
            "baseline_size": len(baseline),
        }

        if len(baseline) < min_tail_size:
            diagnostic.update({"status": "skipped", "reason": "insufficient_tail"})
            diagnostics.append(diagnostic)
            continue

        values = np.array([_enthalpy_per_atom(rec) for rec in baseline], dtype=float)
        median = float(np.median(values))
        mad = float(np.median(np.abs(values - median)))
        robust_sigma = 1.4826 * mad

        diagnostic.update(
            {
                "median": median,
                "mad": mad,
                "robust_sigma": robust_sigma,
            }
        )

        if robust_sigma <= 0.0:
            diagnostic.update({"status": "skipped", "reason": "zero_mad"})
            diagnostics.append(diagnostic)
            continue

        cutoff = median - sigma_factor * robust_sigma
        formula_rejected = [
            rec for rec in group if _enthalpy_per_atom(rec) < cutoff
        ]
        rejected_ids.update(id(rec) for rec in formula_rejected)

        diagnostic.update(
            {
                "status": "applied",
                "cutoff": cutoff,
                "rejected_count": len(formula_rejected),
                "rejected_labels": [rec.label for rec in formula_rejected],
            }
        )
        diagnostics.append(diagnostic)

    kept = [rec for rec in records if id(rec) not in rejected_ids]
    rejected = [rec for rec in records if id(rec) in rejected_ids]
    return kept, rejected, diagnostics


def rank_structures(
    records: list[StructureRecord],
    delta_e: float | None = None,
    top_n: int | None = None,
    absolute: bool = False,
) -> list[dict]:
    """Rank structures by enthalpy per formula unit.

    Returns a list of output dicts with keys needed for formatting.
    """

    # Group by composition
    groups: dict[str, list[StructureRecord]] = {}
    for rec in records:
        key = rec.reduced_formula
        groups.setdefault(key, []).append(rec)

    output_records: list[dict] = []
    for formula, group in groups.items():
        min_h_per_fu = min(r.enthalpy_per_fu for r in group)

        for rec in group:
            rel_h = rec.enthalpy_per_fu - min_h_per_fu
            rel_h_per_atom = (
                rel_h * rec.n_formula_units / rec.natoms if rec.natoms > 0 else 0.0
            )

            # Filter by delta_e (per atom)
            if delta_e is not None and rel_h_per_atom > delta_e:
                continue

            output_records.append(
                {
                    "label": rec.label,
                    "pressure": rec.pressure,
                    "volume_per_fu": rec.volume_per_fu,
                    "enthalpy_per_fu": rec.enthalpy_per_fu,
                    "relative_enthalpy": rel_h,
                    "relative_enthalpy_per_atom": rel_h_per_atom,
                    "spin_per_fu": rec.spin / rec.n_formula_units
                    if rec.n_formula_units > 0
                    else 0.0,
                    "spin_abs_per_fu": rec.spin_abs / rec.n_formula_units
                    if rec.n_formula_units > 0
                    else 0.0,
                    "nfu": rec.n_formula_units,
                    "formula": formula,
                    "symm": rec.symm,
                    "copies": rec.copies,
                    "source": rec.source,
                    "_record": rec,
                }
            )

    # Sort by enthalpy per formula unit (most stable first)
    output_records.sort(key=lambda r: r["enthalpy_per_fu"])

    # Set display enthalpy: first entry absolute, rest relative (unless -nr)
    if not absolute:
        for rec in output_records:
            rec["display_enthalpy"] = rec["relative_enthalpy"]
        if output_records:
            output_records[0]["display_enthalpy"] = output_records[0]["enthalpy_per_fu"]
    else:
        for rec in output_records:
            rec["display_enthalpy"] = rec["enthalpy_per_fu"]

    if top_n is not None:
        output_records = output_records[:top_n]

    return output_records


def summary_structures(
    records: list[StructureRecord],
    delta_e: float | None = None,
) -> list[dict]:
    """Return only the most stable structure per composition.

    Output matches cryan's ``-s`` flag.
    """
    # Group by composition
    groups: dict[str, list[StructureRecord]] = {}
    for rec in records:
        key = rec.reduced_formula
        groups.setdefault(key, []).append(rec)

    output_records: list[dict] = []
    total_copies = 0

    for formula, group in groups.items():
        # Find most stable
        best = min(group, key=lambda r: r.enthalpy_per_fu)
        total_copies += sum(r.copies for r in group)
        group_copies = sum(r.copies for r in group)

        # Count how many distinct structures in this composition
        n_structures = len(group)

        output_records.append(
            {
                "label": best.label,
                "pressure": best.pressure,
                "volume_per_fu": best.volume_per_fu,
                "enthalpy_per_fu": best.enthalpy_per_fu,
                "spin_per_fu": best.spin / best.n_formula_units
                if best.n_formula_units > 0
                else 0.0,
                "spin_abs_per_fu": best.spin_abs / best.n_formula_units
                if best.n_formula_units > 0
                else 0.0,
                "nfu": best.n_formula_units,
                "formula": formula,
                "symm": best.symm,
                "copies": best.copies,
                "group_copies": group_copies,
                "n_structures": n_structures,
                "source": best.source,
                "_record": best,
            }
        )

    # Sort by enthalpy descending (most negative = most stable last in cryan)
    # Actually cryan sorts by -energy (descending energy = most stable first)
    output_records.sort(key=lambda r: -r["enthalpy_per_fu"])

    # Set display enthalpy (always absolute for summary)
    for rec in output_records:
        rec["display_enthalpy"] = rec["enthalpy_per_fu"]

    return output_records, total_copies


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------


def format_header(
    show_spin: bool = False,
    summary_mode: bool = False,
    long_labels: bool = False,
) -> str:
    """Format the header line (goes to stderr)."""
    struct_fmt = f"{'structure':>40s}" if long_labels else f"{'structure':<20s}"

    if summary_mode:
        parts = [
            struct_fmt,
            f"{'P/GPa':>9s}",
            f"{'V/A^3':>10s}",
            f"{'H/eV':>12s}",
        ]
        if show_spin:
            parts.append(f"{'S':>6s}")
            parts.append(f"{'|S|':>6s}")
        parts.extend(
            [
                f"{'nfu':>4s}",
                f"{'formula':>18s}",
                f"{'space_group':>11s}",
                f"{'#':>6s}",
                f"{'tot#':>6s}",
            ]
        )
        return " ".join(parts)

    parts = [
        struct_fmt,
        f"{'P/GPa':>9s}",
        f"{'V/A^3':>10s}",
        f"{'H/eV':>12s}",
    ]
    if show_spin:
        parts.append(f"{'S':>6s}")
        parts.append(f"{'|S|':>6s}")
    parts.extend(
        [
            f"{'nfu':>4s}",
            f"{'formula':>18s}",
            f"{'space_group':>11s}",
            f"{'#':>5s}",
        ]
    )
    return " ".join(parts)


def format_rank_line(
    rec: dict,
    long_labels: bool = False,
    show_spin: bool = False,
    summary_mode: bool = False,
) -> str:
    """Format a single ranked record as a cryan-compatible output line."""
    label = rec["label"] if long_labels else _truncate_label(rec["label"])

    parts = [
        f"{label:>40s}" if long_labels else f"{label:<20s}",
        f"{rec['pressure']:>9.2f}",
        f"{rec['volume_per_fu']:>10.3f}",
        f"{rec['display_enthalpy']:>12.6f}",
    ]
    if show_spin:
        parts.append(f"{rec.get('spin_per_fu', 0.0):>6.2f}")
        parts.append(f"{rec.get('spin_abs_per_fu', 0.0):>6.2f}")
    parts.extend(
        [
            f"{rec['nfu']:>4d}",
            f"{rec['formula']:>18s}",
            f"{rec['symm'].strip('()'):>11s}",
            f"{rec['copies']:>5d}",
        ]
    )
    if summary_mode:
        parts.append(f"{rec.get('group_copies', rec['copies']):>6d}")

    return " ".join(parts)


# ---------------------------------------------------------------------------
# Maxwell plotting
# ---------------------------------------------------------------------------


def plot_maxwell(
    ranked: list[dict],
    elements: list[str],
) -> object:
    """Build a cryan-style convex hull plot using plotly.

    Returns a plotly ``go.Figure``.  The layout has a main panel showing
    formation enthalpy vs composition with all structures as scatter
    points and the convex hull as lines, plus a right-hand panel listing
    the stable phases (formula, nfu, space group, composition).

    Args:
        ranked: Output dicts from ``maxwell_construction()``.
        elements: Element list (2 for binary).
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    if len(elements) != 2:
        raise ValueError("Custom plot only supported for binary systems")

    el_b, el_a = elements  # sorted by atomic number
    x_label = f"x in {el_b}<sub>1-x</sub>{el_a}<sub>x</sub>"

    # --- Collect scatter data ---
    x_all, y_all, colors_all, labels_all = [], [], [], []
    stable_seen: dict[str, dict] = {}

    for rec in ranked:
        sc = rec["species_counts"]
        total = sum(sc.values())
        x = sc.get(el_a, 0) / total if total > 0 else 0.0
        y = rec["formation_energy_per_atom"]

        x_all.append(x)
        y_all.append(y)
        labels_all.append(rec["label"])
        colors_all.append("black" if rec["on_hull"] else "red")

        if rec["on_hull"]:
            formula = rec["formula"]
            symm = rec["symm"].strip("()")
            if formula not in stable_seen or y < stable_seen[formula]["y"]:
                stable_seen[formula] = {
                    "formula": formula,
                    "nfu": rec["nfu"],
                    "symm": symm,
                    "x": x,
                    "y": y,
                }

    stable_entries = sorted(stable_seen.values(), key=lambda s: s["x"])

    # --- Build figure with two columns ---
    fig = make_subplots(
        rows=1,
        cols=2,
        column_widths=[0.72, 0.28],
        horizontal_spacing=0.03,
        specs=[[{"type": "xy"}, {"type": "domain"}]],
        print_grid=False,
    )

    # Main panel: all structures
    fig.add_trace(
        go.Scatter(
            x=x_all,
            y=y_all,
            mode="markers",
            marker={"size": 5, "color": colors_all},
            text=labels_all,
            hovertemplate="x=%{x:.3f}<br>E<sub>f</sub>=%{y:.4f}<br>%{text}<extra></extra>",
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    # Hull line over stable entries
    if len(stable_entries) >= 2:
        fig.add_trace(
            go.Scatter(
                x=[s["x"] for s in stable_entries],
                y=[s["y"] for s in stable_entries],
                mode="lines+markers",
                line={"color": "black", "width": 2},
                marker={"size": 7, "color": "black", "symbol": "circle"},
                text=[s["formula"] for s in stable_entries],
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            ),
            row=1,
            col=1,
        )

    # Right panel: stable phases table
    fig.add_trace(
        go.Table(
            header={
                "values": ["Phase", "nfu", "Space Group", "x"],
                "fill_color": "lightgrey",
                "align": "left",
                "font": {"size": 12},
                "line_color": "black",
            },
            cells={
                "values": [
                    [s["formula"] for s in stable_entries],
                    [str(s["nfu"]) for s in stable_entries],
                    [s["symm"] for s in stable_entries],
                    [f"{s['x']:.3f}" for s in stable_entries],
                ],
                "fill_color": "white",
                "align": "left",
                "font": {"size": 11},
                "line_color": "black",
            },
        ),
        row=1,
        col=2,
    )

    # Styling
    fig.update_xaxes(title_text=x_label, range=[-0.02, 1.02], row=1, col=1)
    fig.update_yaxes(title_text="Formation Enthalpy (eV/atom)", row=1, col=1)
    fig.update_xaxes(showticklabels=False, showgrid=False, row=1, col=2)
    fig.update_yaxes(showticklabels=False, showgrid=False, row=1, col=2)

    fig.update_layout(
        title=f"Convex Hull: {el_b}-{el_a}",
        height=600,
        width=900,
        margin={"l": 60, "r": 20, "t": 50, "b": 60},
    )

    return fig


# ---------------------------------------------------------------------------
# Maxwell formatting
# ---------------------------------------------------------------------------


def format_maxwell_header(
    show_spin: bool = False,
    long_labels: bool = False,
) -> str:
    """Format the header line for Maxwell construction output."""
    struct_fmt = f"{'structure':>40s}" if long_labels else f"{'structure':<20s}"

    parts = [
        struct_fmt,
        f"{'P/GPa':>9s}",
        f"{'V/A^3':>10s}",
        f"{'H/eV/atom':>12s}",
        f"{'hull(eV/at)':>12s}",
        f"{'e_hull(eV)':>11s}",
        f"{'st':>3s}",
    ]
    if show_spin:
        parts.append(f"{'S':>6s}")
        parts.append(f"{'|S|':>6s}")
    parts.extend(
        [
            f"{'nfu':>4s}",
            f"{'formula':>18s}",
            f"{'space_group':>11s}",
            f"{'#':>5s}",
        ]
    )
    return " ".join(parts)


def format_maxwell_line(
    rec: dict,
    long_labels: bool = False,
    show_spin: bool = False,
) -> str:
    """Format a single Maxwell construction record as a cryan-compatible output line."""
    label = rec["label"] if long_labels else _truncate_label(rec["label"])
    status = "+" if rec["on_hull"] else "-"

    parts = [
        f"{label:>40s}" if long_labels else f"{label:<20s}",
        f"{rec['pressure']:>9.2f}",
        f"{rec['volume_per_fu']:>10.3f}",
        f"{rec['enthalpy_per_atom']:>12.6f}",
        f"{rec['hull_energy_per_atom']:>12.6f}",
        f"{rec['e_above_hull']:>11.6f}",
        f"{status:>3s}",
    ]
    if show_spin:
        parts.append(f"{rec.get('spin_per_fu', 0.0):>6.2f}")
        parts.append(f"{rec.get('spin_abs_per_fu', 0.0):>6.2f}")
    parts.extend(
        [
            f"{rec['nfu']:>4d}",
            f"{rec['formula']:>18s}",
            f"{rec['symm'].strip('()'):>11s}",
            f"{rec['copies']:>5d}",
        ]
    )
    return " ".join(parts)

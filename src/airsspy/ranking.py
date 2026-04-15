"""Structure ranking utilities for AIRSS search results.

Provides fast parsing of SHELX .res files, ranking by enthalpy per
formula unit, and optional merging of similar structures using
distance fingerprint comparison (equivalent to ``cryan -u``).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from functools import reduce
from math import gcd
from typing import TextIO


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


def read_extxyz_file(path: str) -> list[StructureRecord]:
    """Read structures from an extxyz file using ASE.

    Energy is expected in ``atoms.info["energy"]``.
    """
    from ase.io import read as ase_read

    records: list[StructureRecord] = []
    try:
        atoms_list = ase_read(path, index=":")
    except Exception:
        atoms_list = [ase_read(path)]

    for i, atoms in enumerate(atoms_list):
        species_counts = dict(Counter(atoms.get_chemical_symbols()))
        energy = atoms.info.get("energy", 0.0)
        label = atoms.info.get("label", atoms.info.get("name", f"{path}:{i}"))
        pressure = atoms.info.get("pressure", atoms.info.get("extern_pressure", 0.0))
        volume = atoms.get_volume()
        natoms = len(atoms)
        spin = atoms.info.get("spin", 0.0)
        spin_abs = atoms.info.get("spin_abs", 0.0)
        symm = atoms.info.get("symm", atoms.info.get("spacegroup", ""))

        records.append(
            StructureRecord(
                label=str(label),
                pressure=float(pressure),
                volume=volume,
                enthalpy=float(energy),
                spin=float(spin),
                spin_abs=float(spin_abs),
                natoms=natoms,
                symm=str(symm),
                species_counts=species_counts,
                source=path,
            )
        )
    return records


# ---------------------------------------------------------------------------
# Similarity / merging (cryan -u equivalent)
# ---------------------------------------------------------------------------


def _compute_distance_fingerprint(
    record: StructureRecord,
    cutoff: float = 4.0,
    zweight: bool = True,
) -> list[float] | None:
    """Compute a sorted distance fingerprint for a structure.

    Uses pymatgen's ``get_all_neighbors(cutoff)`` to find all distances
    to periodic images within *cutoff*, matching cryan's distance
    fingerprint algorithm.  When *zweight* is True, each distance *d* is
    weighted as ``d * zmax**2 / (Z_i * Z_j)`` to distinguish different
    atom-type pairs.

    Returns None if the structure cannot be loaded.
    """
    from .restools import RESFile

    if record._raw_lines:
        try:
            res = RESFile.from_lines(record._raw_lines, include_structure=True)
            if res.structure is None:
                return None
        except Exception:
            return None
    else:
        return None

    structure = res.structure
    neighbors = structure.get_all_neighbors(cutoff)

    if not zweight:
        all_dists: list[float] = []
        for nlist in neighbors:
            for n in nlist:
                all_dists.append(float(n.nn_distance))
    else:
        zmax = max(site.specie.Z for site in structure)
        zmax2 = zmax * zmax
        all_dists = []
        for i, nlist in enumerate(neighbors):
            zi = structure[i].specie.Z
            for n in nlist:
                zj = n.specie.Z
                all_dists.append(float(n.nn_distance * zmax2 / (zi * zj)))

    all_dists.sort()
    return all_dists


def eliminate_similar(
    records: list[StructureRecord],
    threshold: float,
    cutoff: float = 4.0,
) -> list[StructureRecord]:
    """Merge similar structures by comparing distance fingerprints.

    Matches cryan's ``-u`` behaviour:
    1. Sort by enthalpy_per_fu ascending (most stable first)
    2. For each pair of same-formula records, compare scaled distance fingerprints
    3. If max difference < threshold * mean_min_distance, merge copies

    *cutoff* controls the neighbour search radius (Å) for fingerprint
    computation (default 4.0, matching cryan's ``rmax / 1.75``).

    Returns the deduplicated list with accumulated copies.
    """
    # Group by formula
    groups: dict[str, list[StructureRecord]] = {}
    for rec in records:
        key = rec.reduced_formula
        groups.setdefault(key, []).append(rec)

    result: list[StructureRecord] = []

    for _formula, group in groups.items():
        # Sort by energy (most stable first)
        group.sort(key=lambda r: r.enthalpy_per_fu)

        # Compute fingerprints
        fingerprints: list[list[float] | None] = []
        for rec in group:
            fp = _compute_distance_fingerprint(rec, cutoff=cutoff)
            fingerprints.append(fp)

        # Track which records are merged into another
        merged_into: list[int] = [-1] * len(group)

        for i in range(len(group)):
            if merged_into[i] >= 0:
                continue
            if fingerprints[i] is None:
                continue

            fi = fingerprints[i]
            nfi = len(fi)
            # Mean minimum distance (smallest non-zero distance)
            min_dist_i = fi[0] if fi else 1.0

            for j in range(i + 1, len(group)):
                if merged_into[j] >= 0:
                    continue
                if fingerprints[j] is None:
                    continue

                fj = fingerprints[j]

                # Scale factors to account for volume differences
                vol_i = group[i].volume_per_fu
                vol_j = group[j].volume_per_fu
                scale_a = ((vol_j + vol_i) / (2.0 * vol_i)) ** (1.0 / 3.0)
                scale_b = ((vol_j + vol_i) / (2.0 * vol_j)) ** (1.0 / 3.0)

                # Compare fingerprints up to min(nfi*form_j, nfj*form_i) entries
                # (matches cryan's cross-formula-unit comparison)
                nfi_rec = group[i].n_formula_units
                nfj_rec = group[j].n_formula_units
                n_compare = min(nfi * nfj_rec, len(fj) * nfi_rec)

                min_dist_j = fj[0] if fj else 1.0
                mean_min = (min_dist_i * scale_a + min_dist_j * scale_b) / 2.0

                # Compute max absolute difference
                max_diff = 0.0
                too_different = False
                for k in range(n_compare):
                    # Map k to indices in fi and fj with formula-unit scaling
                    ii = k // nfj_rec
                    jj = k // nfi_rec
                    if ii >= nfi or jj >= len(fj):
                        break
                    if fi[ii] < 1e-10:
                        continue
                    diff = abs(fi[ii] * scale_a - fj[jj] * scale_b)
                    if diff > threshold * mean_min:
                        too_different = True
                        break
                    if diff > max_diff:
                        max_diff = diff

                if not too_different:
                    merged_into[j] = i

        # Accumulate copies
        for i in range(len(group)):
            if merged_into[i] >= 0:
                target = merged_into[i]
                # Walk to root
                while merged_into[target] >= 0:
                    target = merged_into[target]
                group[target].copies += group[i].copies

        # Keep only unmerged records
        for i in range(len(group)):
            if merged_into[i] < 0:
                result.append(group[i])

    return result


# ---------------------------------------------------------------------------
# Ranking and output
# ---------------------------------------------------------------------------


def rank_structures(
    records: list[StructureRecord],
    delta_e: float | None = None,
    formula_filter: str | None = None,
    top_n: int | None = None,
    absolute: bool = False,
) -> list[dict]:
    """Rank structures by enthalpy per formula unit.

    Returns a list of output dicts with keys needed for formatting.
    """
    if formula_filter:
        records = [r for r in records if r.reduced_formula == formula_filter]

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
                f"{'nfu':>7s}",
                f"{'formula':>18s}",
                f"{'space group':>11s}",
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
            f"{'nfu':>7s}",
            f"{'formula':>18s}",
            f"{'space group':>11s}",
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
    label = rec["label"] if long_labels else rec["label"][:20]

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
            f"{rec['nfu']:>7d}",
            f"{rec['formula']:>18s}",
            f"{rec['symm'].strip('()'):>11s}",
            f"{rec['copies']:>5d}",
        ]
    )
    if summary_mode:
        parts.append(f"{rec.get('group_copies', rec['copies']):>6d}")

    return " ".join(parts)

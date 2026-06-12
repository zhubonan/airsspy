"""Reusable helpers for AIRSS search workflows."""

from __future__ import annotations

import math
import random
import re
from collections import Counter
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np
from ase import Atoms
from pymatgen.core import Composition
from pymatgen.io.ase import AseAtomsAdaptor

from .ranking import (
    StructureRecord,
    _compute_distance_fingerprint,
)
from .restools import RESFile

DEFAULT_FORMULA_REMOVE_DIRECTIVES = (
    "NATOM",
    "SPECIES",
    "FORMULA",
    "VARVOL",
    "TARGVOL",
)


@dataclass
class FormulaSamplingOptions:
    """Options for sampling buildcell ``#FORMULA`` directives."""

    formulas: Sequence[str] = ()
    elements: Sequence[str] = ()
    max_coeff: int = 6
    target_atom_volumes: dict[str, float] = field(default_factory=dict)
    oxidation_states: dict[str, Sequence[int]] = field(default_factory=dict)
    require_charge_neutral: bool = True
    remove_directives: Sequence[str] = DEFAULT_FORMULA_REMOVE_DIRECTIVES


@dataclass
class FormulaSamplingContext:
    """Precomputed state for formula sampling."""

    formulas: list[str]
    remove_directives: tuple[str, ...] = DEFAULT_FORMULA_REMOVE_DIRECTIVES
    varvol_by_formula: dict[str, float] = field(default_factory=dict)

    def sample(
        self,
        seed_text: str,
        rng: random.Random | None = None,
    ) -> tuple[str, str, float | None]:
        """Return ``(new_seed_text, formula, varvol)`` for one sampled formula."""
        if not self.formulas:
            raise ValueError("No formulas are available for sampling")
        chooser = rng.choice if rng is not None else random.choice
        formula = chooser(self.formulas)
        varvol = self.varvol_by_formula.get(formula)
        return (
            inject_formula_directive(
                seed_text,
                formula,
                varvol=varvol,
                remove_directives=self.remove_directives,
            ),
            formula,
            varvol,
        )


@dataclass
class RssPruneOptions:
    """Options for pruning relaxed RSS candidates before final output."""

    enabled: bool = False
    pool_size: int = 100
    keep_fraction: float = 0.10
    dedup_tol: float = 0.10
    fingerprint_cutoff: float = 5.0
    min_stable_pool_size: int = 50
    stable_window: int = 20
    mean_abs_tol: float = 1e-3
    median_abs_tol: float = 1e-3
    zweight: bool = False


@dataclass
class RssCandidate:
    """A relaxed structure considered for pruning."""

    label: str
    energy: float
    atoms: Atoms | None = None
    res_path: Path | None = None
    record: StructureRecord | None = None
    source: Any = None

    @property
    def natoms(self) -> int:
        """Number of atoms in this candidate."""
        if self.record is not None:
            return self.record.natoms
        if self.atoms is not None:
            return len(self.atoms)
        return 0

    @property
    def energy_per_atom(self) -> float:
        """Energy or enthalpy per atom."""
        return self.energy / self.natoms if self.natoms > 0 else self.energy

    @property
    def reduced_formula(self) -> str:
        """Reduced formula for the candidate."""
        return _candidate_record(self).reduced_formula


def canonicalize_formula(formula: str) -> str:
    """Return pymatgen's canonical reduced formula string."""
    return Composition(formula).reduced_formula


def inject_formula_directive(
    seed_text: str,
    formula: str,
    varvol: float | None = None,
    remove_directives: Sequence[str] = DEFAULT_FORMULA_REMOVE_DIRECTIVES,
) -> str:
    """Inject ``#FORMULA`` and optional ``#VARVOL`` directives into seed text."""
    remove_patterns = [
        re.compile(rf"^\s*#{re.escape(key)}\s*=", re.IGNORECASE)
        for key in remove_directives
    ]
    filtered = [
        line
        for line in seed_text.splitlines()
        if not any(pattern.search(line) for pattern in remove_patterns)
    ]

    lines = [f"#FORMULA={canonicalize_formula(formula)}"]
    if varvol is not None:
        lines.append(f"#VARVOL={_format_seed_number(varvol)}")
    lines.extend(filtered)
    return "\n".join(lines)


def build_formula_sampling_context(
    options: FormulaSamplingOptions,
    seed_text: str | None = None,
) -> FormulaSamplingContext:
    """Build a formula sampling context from user options."""
    formulas = _resolve_formula_pool(options, seed_text=seed_text)
    target_volumes = {
        str(key): float(value) for key, value in options.target_atom_volumes.items()
    }
    varvol_by_formula = (
        _build_varvol_map(formulas, target_volumes) if target_volumes else {}
    )
    return FormulaSamplingContext(
        formulas=formulas,
        remove_directives=tuple(options.remove_directives),
        varvol_by_formula=varvol_by_formula,
    )


def make_seed_text_transform(
    context: FormulaSamplingContext,
    rng: random.Random | None = None,
) -> Callable[[str], str]:
    """Return a callable that rewrites seed text for one sampled formula."""

    def transform(seed_text: str) -> str:
        sampled_seed, _, _ = context.sample(seed_text, rng=rng)
        return sampled_seed

    return transform


def parse_key_float(text: str) -> tuple[str, float]:
    """Parse ``KEY=value`` into a string key and float value."""
    key, value = _split_key_value(text)
    return key, float(value)


def parse_key_ints(text: str) -> tuple[str, list[int]]:
    """Parse ``KEY=a,b`` into a string key and list of integers."""
    key, value = _split_key_value(text)
    return key, [int(item.strip()) for item in value.split(",") if item.strip()]


def validate_prune_options(options: RssPruneOptions) -> RssPruneOptions:
    """Validate pruning options and return *options* for convenient chaining."""
    if options.pool_size < 1:
        raise ValueError("pool_size must be >= 1")
    if (
        not math.isfinite(options.keep_fraction)
        or not 0.0 < options.keep_fraction <= 1.0
    ):
        raise ValueError("keep_fraction must be in (0, 1]")
    if not math.isfinite(options.dedup_tol) or options.dedup_tol < 0.0:
        raise ValueError("dedup_tol must be finite and >= 0")
    if (
        not math.isfinite(options.fingerprint_cutoff)
        or options.fingerprint_cutoff <= 0.0
    ):
        raise ValueError("fingerprint_cutoff must be finite and > 0")
    if options.min_stable_pool_size < 1:
        raise ValueError("min_stable_pool_size must be >= 1")
    if options.stable_window < 0:
        raise ValueError("stable_window must be >= 0")
    if not math.isfinite(options.mean_abs_tol) or options.mean_abs_tol < 0.0:
        raise ValueError("mean_abs_tol must be finite and >= 0")
    if not math.isfinite(options.median_abs_tol) or options.median_abs_tol < 0.0:
        raise ValueError("median_abs_tol must be finite and >= 0")
    return options


def pool_statistics(candidates: Sequence[RssCandidate]) -> dict[str, float]:
    """Return mean and median energy-per-atom statistics for a candidate pool."""
    values = [candidate.energy_per_atom for candidate in candidates]
    if not values:
        raise ValueError("Cannot compute statistics for an empty candidate pool")
    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
    }


def should_flush_prune_pool(
    candidates: Sequence[RssCandidate],
    stats_history: Sequence[dict[str, float]],
    options: RssPruneOptions,
) -> tuple[str, bool]:
    """Return ``(reason, should_flush)`` for a pruning pool."""
    if len(candidates) >= max(options.pool_size, 1):
        return "pool_size", True
    if options.stable_window <= 0:
        return "none", False
    if len(candidates) < options.min_stable_pool_size:
        return "none", False
    if len(stats_history) <= options.stable_window:
        return "none", False

    current = stats_history[-1]
    previous = stats_history[-options.stable_window - 1]
    mean_stable = abs(current["mean"] - previous["mean"]) <= options.mean_abs_tol
    median_stable = (
        abs(current["median"] - previous["median"]) <= options.median_abs_tol
    )
    return "stable_statistics", mean_stable and median_stable


def select_pruned_candidates(
    candidates: Sequence[RssCandidate],
    options: RssPruneOptions,
    remaining: int | None = None,
) -> tuple[list[RssCandidate], list[RssCandidate]]:
    """Select low-energy unique candidates and return ``(kept, rejected)``."""
    validate_prune_options(options)
    if not candidates:
        return [], []

    limit = len(candidates) if remaining is None else max(remaining, 0)
    if limit <= 0:
        return [], list(candidates)

    kept: list[RssCandidate] = []
    rejected: list[RssCandidate] = []
    for group in _group_candidates(candidates).values():
        quota = min(
            limit - len(kept),
            len(group),
            max(1, math.ceil(options.keep_fraction * len(group))),
        )
        if quota <= 0:
            rejected.extend(group)
            continue
        selected, not_selected = _select_group(group, options, quota)
        kept.extend(selected)
        rejected.extend(not_selected)
        if len(kept) >= limit:
            seen_ids = {id(candidate) for candidate in kept + rejected}
            rejected.extend(cand for cand in candidates if id(cand) not in seen_ids)
            break

    kept.sort(key=lambda cand: cand.energy_per_atom)
    rejected.sort(key=lambda cand: cand.energy_per_atom)
    return kept, rejected


def prune_relaxed_pool(
    candidates: Sequence[RssCandidate],
    options: RssPruneOptions,
    remaining: int | None = None,
) -> list[RssCandidate]:
    """Return the candidates kept after post-relax pruning."""
    kept, _ = select_pruned_candidates(candidates, options, remaining=remaining)
    return kept


def candidate_from_res(path: str | Path) -> RssCandidate:
    """Create a pruning candidate from an AIRSS ``.res`` file."""
    res_path = Path(path)
    res = RESFile.from_file(str(res_path), include_structure=True)
    if res.label is None or res.enthalpy is None:
        raise ValueError(f"Cannot create RSS candidate from {res_path}")

    species_counts = Counter(str(site.specie.symbol) for site in res.structure)
    record = StructureRecord(
        label=res.label,
        pressure=float(res.pressure or 0.0),
        volume=float(res.volume or 0.0),
        enthalpy=float(res.enthalpy),
        spin=float(res.spin or 0.0),
        spin_abs=float(res.spin_abs or 0.0),
        natoms=int(res.natoms or len(res.structure)),
        symm=str(res.symm or ""),
        species_counts=dict(species_counts),
        source=str(res_path),
        _raw_lines=res.lines,
    )
    return RssCandidate(
        label=res.label,
        energy=float(res.enthalpy),
        res_path=res_path,
        record=record,
        source=res,
    )


def _split_key_value(text: str) -> tuple[str, str]:
    if "=" not in text:
        raise ValueError(f"Expected KEY=value, got {text!r}")
    key, value = text.split("=", 1)
    key = key.strip()
    value = value.strip()
    if not key or not value:
        raise ValueError(f"Expected KEY=value, got {text!r}")
    return key, value


def _format_seed_number(value: float) -> str:
    value = float(value)
    nearest = round(value)
    if math.isclose(
        value, nearest, rel_tol=math.sqrt(np.finfo(float).eps), abs_tol=0.0
    ):
        return str(int(nearest))
    return str(value)


def _resolve_formula_pool(
    options: FormulaSamplingOptions,
    seed_text: str | None = None,
) -> list[str]:
    if options.max_coeff < 1:
        raise ValueError("max_coeff must be >= 1")

    if options.formulas:
        formulas = _canonical_formula_list(options.formulas)
    else:
        if not options.elements:
            raise ValueError("Either formulas or elements must be provided")
        formulas = _enumerate_reduced_formulas(options.elements, options.max_coeff)

    if options.elements:
        allowed = set(options.elements)
        formulas = [
            formula
            for formula in formulas
            if set(Composition(formula).as_dict()).issubset(allowed)
        ]

    if seed_text is not None:
        constraints = _extract_seed_constraints(seed_text)
        formulas = [
            formula
            for formula in formulas
            if _formula_fits_seed_constraints(formula, constraints)
        ]

    if options.require_charge_neutral and options.oxidation_states:
        formulas = [
            formula
            for formula in formulas
            if _has_neutral_oxidation_state(formula, options.oxidation_states)
        ]

    if not formulas:
        raise ValueError("No valid formulas remain after filtering")
    return sorted(formulas)


def _canonical_formula_list(formulas: Sequence[str]) -> list[str]:
    return sorted(
        {
            canonicalize_formula(formula.strip())
            for formula in formulas
            if formula.strip()
        }
    )


def _enumerate_reduced_formulas(elements: Sequence[str], max_coeff: int) -> list[str]:
    clean_elements = [str(element) for element in elements]
    formulas: set[str] = set()
    for coeffs in product(range(max_coeff + 1), repeat=len(clean_elements)):
        if all(coeff == 0 for coeff in coeffs):
            continue
        parts = {
            element: coeff
            for element, coeff in zip(clean_elements, coeffs)
            if coeff > 0
        }
        formulas.add(Composition(parts).reduced_formula)
    return sorted(formulas)


def _extract_seed_constraints(seed_text: str) -> dict[str, tuple[int, int] | None]:
    return {
        "natom": _parse_seed_int_range(_extract_seed_directive(seed_text, "NATOM")),
        "nform": _parse_seed_int_range(_extract_seed_directive(seed_text, "NFORM"))
        or (1, 10**9),
    }


def _extract_seed_directive(seed_text: str, key: str) -> str | None:
    pattern = re.compile(rf"^\s*#{re.escape(key)}\s*=\s*(.*)", re.IGNORECASE)
    for line in seed_text.splitlines():
        match = pattern.search(line)
        if match:
            return match.group(1).strip()
    return None


def _parse_seed_int_range(text: str | None) -> tuple[int, int] | None:
    if text is None:
        return None
    cleaned = text.strip().replace("{", "").replace("}", "")
    if not cleaned:
        return None
    if "-" in cleaned:
        match = re.match(r"^\s*(-?\d+)\s*-\s*(-?\d+)\s*$", cleaned)
        if match:
            first = int(match.group(1))
            second = int(match.group(2))
            return min(first, second), max(first, second)
    if "," in cleaned:
        values = [int(item.strip()) for item in cleaned.split(",") if item.strip()]
        if values:
            return min(values), max(values)
        return None
    if re.match(r"^-?\d+$", cleaned):
        value = int(cleaned)
        return value, value
    return None


def _formula_fits_seed_constraints(
    formula: str,
    constraints: dict[str, tuple[int, int] | None],
) -> bool:
    formula_atoms = int(sum(Composition(formula).as_dict().values()))
    if formula_atoms <= 0:
        return False
    nform_min, nform_max = constraints["nform"] or (1, 10**9)
    natom_range = constraints["natom"]
    if natom_range is None:
        return nform_min <= nform_max
    natom_min, natom_max = natom_range
    lower = max(nform_min, math.ceil(natom_min / formula_atoms))
    upper = min(nform_max, math.floor(natom_max / formula_atoms))
    return lower <= upper


def _has_neutral_oxidation_state(
    formula: str,
    oxidation_states: dict[str, Sequence[int]],
) -> bool:
    composition = Composition(formula).as_dict()
    elements = sorted(composition)
    ox_lists: list[Sequence[int]] = []
    counts: list[int] = []
    for element in elements:
        states = oxidation_states.get(element)
        if not states:
            return False
        ox_lists.append(states)
        counts.append(int(composition[element]))

    def search(index: int, charge: int) -> bool:
        if index == len(elements):
            return charge == 0
        return any(
            search(index + 1, charge + counts[index] * ox) for ox in ox_lists[index]
        )

    return search(0, 0)


def _build_varvol_map(
    formulas: Sequence[str],
    target_atom_volumes: dict[str, float],
) -> dict[str, float]:
    result: dict[str, float] = {}
    for formula in formulas:
        composition = Composition(formula).as_dict()
        missing = sorted(
            element for element in composition if element not in target_atom_volumes
        )
        if missing:
            missing_text = ", ".join(missing)
            raise ValueError(
                f"Missing target atom volumes for {formula}: {missing_text}"
            )
        total_volume = sum(
            amount * target_atom_volumes[element]
            for element, amount in composition.items()
        )
        result[formula] = total_volume / sum(composition.values()) * len(composition)
    return result


def _candidate_record(candidate: RssCandidate) -> StructureRecord:
    if candidate.record is not None:
        return candidate.record
    if candidate.atoms is None:
        raise ValueError(f"Candidate {candidate.label!r} has no structure data")

    species_counts = Counter(candidate.atoms.get_chemical_symbols())
    structure = AseAtomsAdaptor.get_structure(candidate.atoms)
    candidate.record = StructureRecord(
        label=candidate.label,
        pressure=float(candidate.atoms.info.get("pressure", 0.0)),
        volume=float(candidate.atoms.get_volume()),
        enthalpy=float(candidate.energy),
        natoms=len(candidate.atoms),
        symm=str(candidate.atoms.info.get("symm", "")),
        species_counts=dict(species_counts),
        source=str(candidate.res_path or ""),
        _atoms=candidate.atoms,
    )
    candidate.record.volume = float(structure.volume)
    return candidate.record


def _group_candidates(
    candidates: Sequence[RssCandidate],
) -> dict[str, list[RssCandidate]]:
    groups: dict[str, list[RssCandidate]] = {}
    for candidate in candidates:
        groups.setdefault(candidate.reduced_formula, []).append(candidate)
    return groups


def _select_group(
    candidates: Sequence[RssCandidate],
    options: RssPruneOptions,
    quota: int,
) -> tuple[list[RssCandidate], list[RssCandidate]]:
    ranked = sorted(candidates, key=lambda cand: cand.energy_per_atom)
    kept: list[RssCandidate] = []
    rejected: list[RssCandidate] = []
    kept_fps: list[np.ndarray] = []

    for candidate in ranked:
        fp = _compute_distance_fingerprint(
            _candidate_record(candidate),
            cutoff=options.fingerprint_cutoff,
            zweight=options.zweight,
        )
        if fp is None:
            kept.append(candidate)
        elif _is_unique_fingerprint(kept_fps, fp, options.dedup_tol):
            kept.append(candidate)
            kept_fps.append(fp)
        else:
            rejected.append(candidate)

        if len(kept) >= quota:
            break

    seen_ids = {id(candidate) for candidate in kept + rejected}
    rejected.extend(candidate for candidate in ranked if id(candidate) not in seen_ids)
    return kept, rejected


def _is_unique_fingerprint(
    references: Sequence[np.ndarray],
    fp: np.ndarray,
    tolerance: float,
) -> bool:
    for ref in references:
        if len(ref) == 0 or len(fp) == 0:
            continue
        n_compare = min(len(ref), len(fp))
        scale = min(ref[0], fp[0]) if min(ref[0], fp[0]) > 0 else 1.0
        distance = float(np.max(np.abs(ref[:n_compare] - fp[:n_compare])) / scale)
        if distance < tolerance:
            return False
    return True

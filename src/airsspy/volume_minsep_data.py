"""Dataset helpers for volume and minsep predictors."""

from __future__ import annotations

import json
import math
import warnings
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from functools import cache
from itertools import combinations_with_replacement
from pathlib import Path
from statistics import median
from typing import Any

from monty.serialization import loadfn
from pymatgen.core import Composition, Structure
from tqdm.auto import tqdm


def _mad(values: list[float], center: float) -> float:
    """Return the median absolute deviation for a non-empty list."""
    return float(median([abs(value - center) for value in values]))


def canonical_pair(symbol_a: str, symbol_b: str) -> tuple[str, str]:
    """Return a deterministic canonical pair tuple."""
    return tuple(sorted((symbol_a, symbol_b)))


def canonical_pair_key(symbol_a: str, symbol_b: str) -> str:
    """Return a canonical pair key formatted as ``A-B``."""
    left, right = canonical_pair(symbol_a, symbol_b)
    return f"{left}-{right}"


@cache
def _parsed_composition(formula: str) -> Composition:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return Composition(formula).reduced_composition


def composition_dict_from_formula(formula: str) -> dict[str, float]:
    """Parse a reduced formula into a composition dictionary."""
    parsed = _parsed_composition(formula)
    return {symbol: float(amount) for symbol, amount in parsed.as_dict().items()}


def formula_num_atoms(formula: str) -> float:
    """Return the number of atoms in one reduced formula unit."""
    return float(sum(composition_dict_from_formula(formula).values()))


def formula_species(formula: str) -> tuple[str, ...]:
    """Return species present in a reduced formula."""
    return tuple(sorted(composition_dict_from_formula(formula)))


def enumerate_formula_pairs(formula: str) -> list[tuple[str, str]]:
    """Enumerate all canonical species pairs implied by a formula."""
    return list(combinations_with_replacement(formula_species(formula), 2))


def load_dataset_records(path: str | Path) -> list[dict]:
    """Load volume/minsep training records from disk."""
    with Path(path).open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, list):
        raise TypeError(
            f"Expected a list of records in {path}, found {type(data).__name__}"
        )
    return data


def volume_per_atom_from_record(record: dict) -> float:
    """Normalize a structure's cell volume by the number of atoms in the cell."""
    composition = {
        symbol: float(amount)
        for symbol, amount in record["composition"].items()
        if not symbol.startswith("@")
    }
    total_atoms = sum(composition.values())
    if total_atoms <= 0:
        raise ValueError(
            f"Invalid atom count for record {record.get('material_id')}: {total_atoms}"
        )
    return float(record["volume"]) / total_atoms


def canonical_minsep_map(record: dict) -> dict[str, float]:
    """Collapse mirrored minsep entries like ``A-B`` and ``B-A`` into one key."""
    values: dict[str, float] = {}
    for key, value in record["minsep"].items():
        left, right = key.split("-", maxsplit=1)
        canonical_key = canonical_pair_key(left, right)
        value = float(value)
        current = values.get(canonical_key)
        if current is None or value < current:
            values[canonical_key] = value
    return values


@dataclass(frozen=True)
class VolumeAggregate:
    """Aggregated volume target for one reduced formula."""

    reduced_formula: str
    num_atoms_in_formula: float
    n_structures: int
    volume_per_atom: float
    volume_per_atom_min: float
    volume_per_atom_max: float
    volume_per_atom_mad: float


@dataclass(frozen=True)
class PairAggregate:
    """Aggregated minsep target for one formula and canonical pair."""

    reduced_formula: str
    pair_key: str
    num_atoms_in_formula: float
    n_observations: int
    minsep: float
    minsep_min: float
    minsep_max: float
    minsep_mad: float


def aggregate_training_targets(
    records: Iterable[dict],
) -> tuple[list[VolumeAggregate], list[PairAggregate], dict]:
    """Aggregate structure-level records into formula-level model targets."""
    volume_groups: dict[str, list[float]] = {}
    pair_groups: dict[tuple[str, str], list[float]] = {}
    pair_key_set: set[str] = set()
    element_set: set[str] = set()
    n_records = 0

    for record in records:
        n_records += 1
        formula = str(record["reduced_formula"])
        composition = composition_dict_from_formula(formula)
        element_set.update(composition)

        volume_groups.setdefault(formula, []).append(
            volume_per_atom_from_record(record)
        )

        for pair_key, value in canonical_minsep_map(record).items():
            pair_key_set.add(pair_key)
            pair_groups.setdefault((formula, pair_key), []).append(value)

    volume_rows: list[VolumeAggregate] = []
    for formula, values in sorted(volume_groups.items()):
        center = float(median(values))
        volume_rows.append(
            VolumeAggregate(
                reduced_formula=formula,
                num_atoms_in_formula=formula_num_atoms(formula),
                n_structures=len(values),
                volume_per_atom=center,
                volume_per_atom_min=float(min(values)),
                volume_per_atom_max=float(max(values)),
                volume_per_atom_mad=_mad(values, center),
            )
        )

    pair_rows: list[PairAggregate] = []
    for (formula, pair_key), values in sorted(pair_groups.items()):
        center = float(median(values))
        pair_rows.append(
            PairAggregate(
                reduced_formula=formula,
                pair_key=pair_key,
                num_atoms_in_formula=formula_num_atoms(formula),
                n_observations=len(values),
                minsep=center,
                minsep_min=float(min(values)),
                minsep_max=float(max(values)),
                minsep_mad=_mad(values, center),
            )
        )

    summary = {
        "n_structure_records": n_records,
        "n_unique_formulas": len(volume_rows),
        "n_pair_rows": len(pair_rows),
        "n_unique_elements": len(element_set),
        "n_unique_pairs": len(pair_key_set),
        "max_species_per_formula": max(
            (len(formula_species(row.reduced_formula)) for row in volume_rows),
            default=0,
        ),
        "max_pairs_per_formula": max(
            (len(enumerate_formula_pairs(row.reduced_formula)) for row in volume_rows),
            default=0,
        ),
        "element_list": sorted(element_set),
        "pair_key_list": sorted(pair_key_set),
        "median_volume_per_atom": float(
            median([row.volume_per_atom for row in volume_rows])
        )
        if volume_rows
        else math.nan,
        "median_minsep": float(median([row.minsep for row in pair_rows]))
        if pair_rows
        else math.nan,
    }

    return volume_rows, pair_rows, summary


def volume_rows_to_dicts(rows: Iterable[VolumeAggregate]) -> list[dict]:
    """Serialize volume aggregate rows to dictionaries."""
    return [asdict(row) for row in rows]


def pair_rows_to_dicts(rows: Iterable[PairAggregate]) -> list[dict]:
    """Serialize pair aggregate rows to dictionaries."""
    return [asdict(row) for row in rows]


def get_minsep_dict(structure: Structure, cutoff: float = 5.0) -> dict[str, float]:
    """Return the minimum observed distance for each ordered species pair."""
    neighbor_list = structure.get_neighbor_list(cutoff)
    symbols = [site.specie.name for site in structure.sites]
    minsep_dict: dict[str, float] = {}
    for left_index, right_index, _image, distance in zip(*neighbor_list):
        symbol_left = symbols[left_index]
        symbol_right = symbols[right_index]
        left_key = f"{symbol_left}-{symbol_right}"
        right_key = f"{symbol_right}-{symbol_left}"
        if left_key not in minsep_dict or float(distance) < minsep_dict[left_key]:
            minsep_dict[left_key] = float(distance)
            minsep_dict[right_key] = float(distance)
    return minsep_dict


def build_dataset_rows(
    docs: Iterable[dict[str, Any]],
    *,
    show_progress: bool = True,
) -> list[dict[str, Any]]:
    """Build raw minsep/volume rows from Materials Project-like documents."""
    rows: list[dict[str, Any]] = []
    items = docs if isinstance(docs, list) else list(docs)
    iterator = tqdm(
        items,
        total=len(items),
        disable=not show_progress,
        desc="Build minsep/volume dataset",
        dynamic_ncols=True,
    )
    for doc in iterator:
        structure = doc["structure"]
        if not isinstance(structure, Structure):
            structure = Structure.from_dict(structure)
        rows.append(
            {
                "minsep": get_minsep_dict(structure),
                "material_id": str(doc["material_id"]),
                "volume": float(structure.volume),
                "composition": structure.composition.as_dict(),
                "reduced_formula": structure.composition.reduced_formula,
                "chemical_formula": structure.composition.formula,
            }
        )
    return rows


def build_minsep_volume_dataset(
    *,
    mp_docs_path: str | Path,
    output_path: str | Path,
    show_progress: bool = True,
) -> dict[str, Any]:
    """Build a raw minsep/volume dataset from a Materials Project document dump."""
    docs = loadfn(str(mp_docs_path))
    rows = build_dataset_rows(docs, show_progress=show_progress)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    return {
        "input_docs": len(rows),
        "output_path": str(output_path.resolve()),
    }


def _load_json(path: str | Path):
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def curate_minsep_volume_dataset(
    *,
    minsep_dataset_path: str | Path,
    mp_docs_path: str | Path,
    output_path: str | Path,
) -> dict:
    """Select the lowest-energy-above-hull row for each reduced formula."""
    rows = _load_json(minsep_dataset_path)
    docs = _load_json(mp_docs_path)
    energy_above_hull = {
        doc["material_id"]: float(doc["energy_above_hull"]) for doc in docs
    }

    by_formula: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        row = dict(row)
        row["energy_above_hull"] = energy_above_hull[row["material_id"]]
        by_formula[row["reduced_formula"]].append(row)

    curated_rows = []
    duplicate_formula_count = 0
    tie_count = 0
    spread_values = []
    for _formula, formula_rows in sorted(by_formula.items()):
        if len(formula_rows) > 1:
            duplicate_formula_count += 1
            sorted_rows = sorted(
                formula_rows,
                key=lambda row: (row["energy_above_hull"], row["material_id"]),
            )
            if (
                abs(
                    sorted_rows[0]["energy_above_hull"]
                    - sorted_rows[1]["energy_above_hull"]
                )
                < 1e-12
            ):
                tie_count += 1
            spread_values.append(
                sorted_rows[-1]["energy_above_hull"]
                - sorted_rows[0]["energy_above_hull"]
            )
            chosen = sorted_rows[0]
        else:
            chosen = formula_rows[0]
        curated_rows.append(chosen)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(curated_rows, handle, indent=2)

    summary = {
        "input_rows": len(rows),
        "curated_rows": len(curated_rows),
        "unique_formulas": len(by_formula),
        "duplicate_formulas": duplicate_formula_count,
        "exact_energy_ties": tie_count,
        "average_duplicate_ehull_spread": (sum(spread_values) / len(spread_values))
        if spread_values
        else 0.0,
        "output_path": str(output_path),
    }
    with output_path.with_suffix(".summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    return summary

"""Volume/minsep estimates and buildcell seed helpers."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from statistics import median

import numpy as np
from ase import Atoms
from pymatgen.core import Composition

from .volume_minsep_data import (
    canonical_minsep_map,
    canonical_pair_key,
    composition_dict_from_formula,
    enumerate_formula_pairs,
    load_dataset_records,
    volume_per_atom_from_record,
)

DEFAULT_VOLUME_SCALE = 1.1
DEFAULT_MINSEP_SCALE_LOW = 0.9
DEFAULT_MINSEP_SCALE_HIGH = 1.1
DEFAULT_MINSEP_DECIMALS = 2
DEFAULT_BASELINE_MINSEP_RANGE = (0.5, 1.0)


@dataclass(frozen=True)
class EstimateSource:
    """Provenance for an estimated value."""

    method: str
    n_observations: int


@dataclass(frozen=True)
class VolumeMinsepEstimate:
    """Formula-level volume and minsep estimate."""

    input_formula: str
    reduced_formula: str
    composition: dict[str, float]
    atoms_per_formula_unit: float
    volume_per_atom: float
    total_volume: float
    canonical_minsep: dict[str, float]
    volume_source: EstimateSource
    minsep_sources: dict[str, EstimateSource]
    formula_match_count: int
    material_id: str | None = None
    energy_above_hull: float | None = None
    predictor_type: str | None = None
    model_artifact: str | None = None
    model_metrics: dict | None = None


def normalize_formula(formula: str) -> str:
    """Return pymatgen's reduced formula for *formula*."""
    return Composition(formula).reduced_formula


def required_pair_keys(formula: str) -> list[str]:
    """Return canonical pair keys required by a formula."""
    return [
        canonical_pair_key(left_symbol, right_symbol)
        for left_symbol, right_symbol in enumerate_formula_pairs(formula)
    ]


@cache
def _cached_dataset_records(path: str) -> tuple[dict, ...]:
    return tuple(load_dataset_records(path))


def load_volume_minsep_dataset_index(dataset_path: str | Path) -> dict[str, dict]:
    """Load a one-row-per-formula curated dataset index."""
    records = _cached_dataset_records(str(Path(dataset_path).resolve()))
    index: dict[str, dict] = {}
    duplicates: list[str] = []
    for record in records:
        formula = str(record["reduced_formula"])
        if formula in index:
            duplicates.append(formula)
            continue
        index[formula] = record
    if duplicates:
        sample = ", ".join(sorted(set(duplicates))[:10])
        raise ValueError(
            f"Expected one row per reduced formula in {dataset_path}, "
            f"found duplicates for: {sample}"
        )
    return index


def volume_minsep_estimate_from_exact_record(
    *,
    formula: str,
    record: dict,
) -> VolumeMinsepEstimate:
    """Build an estimate from one curated dataset record."""
    reduced_formula = normalize_formula(formula)
    record_formula = normalize_formula(str(record["reduced_formula"]))
    if record_formula != reduced_formula:
        raise ValueError(
            f"Dataset row formula {record_formula!r} does not match requested "
            f"formula {reduced_formula!r}"
        )

    composition = composition_dict_from_formula(reduced_formula)
    atoms_per_formula_unit = float(sum(composition.values()))
    canonical_minsep = canonical_minsep_map(record)
    required_pairs = required_pair_keys(reduced_formula)
    missing_pairs = [pair_key for pair_key in required_pairs if pair_key not in canonical_minsep]
    if missing_pairs:
        joined_pairs = ", ".join(sorted(missing_pairs))
        raise ValueError(
            f"Exact dataset row for {reduced_formula} is missing minsep "
            f"observations for required pairs: {joined_pairs}"
        )

    volume_per_atom = volume_per_atom_from_record(record)
    return VolumeMinsepEstimate(
        input_formula=formula,
        reduced_formula=reduced_formula,
        composition=composition,
        atoms_per_formula_unit=atoms_per_formula_unit,
        volume_per_atom=volume_per_atom,
        total_volume=volume_per_atom * atoms_per_formula_unit,
        canonical_minsep={pair_key: canonical_minsep[pair_key] for pair_key in required_pairs},
        volume_source=EstimateSource(method="exact_curated_row", n_observations=1),
        minsep_sources={
            pair_key: EstimateSource(method="exact_curated_row", n_observations=1)
            for pair_key in required_pairs
        },
        formula_match_count=1,
        material_id=str(record.get("material_id"))
        if record.get("material_id") is not None
        else None,
        energy_above_hull=float(record["energy_above_hull"])
        if record.get("energy_above_hull") is not None
        else None,
    )


def lookup_exact_volume_minsep_estimate(
    *,
    formula: str,
    dataset_path: str | Path,
) -> VolumeMinsepEstimate:
    """Look up an exact curated dataset estimate for *formula*."""
    reduced_formula = normalize_formula(formula)
    index = load_volume_minsep_dataset_index(dataset_path)
    try:
        record = index[reduced_formula]
    except KeyError as exc:
        raise ValueError(
            f"Could not find reduced formula {reduced_formula!r} in {dataset_path}"
        ) from exc
    return volume_minsep_estimate_from_exact_record(formula=formula, record=record)


def volume_minsep_estimate_from_prediction(
    *,
    formula: str,
    prediction: dict,
) -> VolumeMinsepEstimate:
    """Convert a baseline prediction payload into an estimate."""
    reduced_formula = normalize_formula(formula)
    composition = composition_dict_from_formula(reduced_formula)
    atoms_per_formula_unit = float(sum(composition.values()))
    required_pairs = required_pair_keys(reduced_formula)
    canonical_minsep = {
        pair_key: float(prediction["canonical_minsep"][pair_key])
        for pair_key in required_pairs
    }
    volume_per_atom = float(prediction["volume_per_atom"])
    predictor_type = str(prediction["predictor_type"])
    return VolumeMinsepEstimate(
        input_formula=formula,
        reduced_formula=reduced_formula,
        composition=composition,
        atoms_per_formula_unit=atoms_per_formula_unit,
        volume_per_atom=volume_per_atom,
        total_volume=volume_per_atom * atoms_per_formula_unit,
        canonical_minsep=canonical_minsep,
        volume_source=EstimateSource(method=f"baseline_{predictor_type}", n_observations=1),
        minsep_sources={
            pair_key: EstimateSource(
                method=f"baseline_{predictor_type}", n_observations=1
            )
            for pair_key in required_pairs
        },
        formula_match_count=1,
        predictor_type=predictor_type,
        model_artifact=str(prediction.get("model_artifact"))
        if prediction.get("model_artifact")
        else None,
        model_metrics=dict(prediction.get("model_metrics") or {}),
    )


def predict_baseline_volume_minsep_estimate(
    *,
    formula: str,
    bundle_path: str | Path,
) -> VolumeMinsepEstimate:
    """Predict an estimate using a baseline bundle."""
    from .volume_minsep_model import BaselineFormulaPredictor

    prediction = BaselineFormulaPredictor(bundle_path=bundle_path).predict(formula)
    return volume_minsep_estimate_from_prediction(formula=formula, prediction=prediction)


def reference_volume_minsep_estimate(
    *,
    formula: str,
    references: Sequence[Atoms],
) -> VolumeMinsepEstimate:
    """Build an estimate from reference ASE structures."""
    if not references:
        raise ValueError(f"No reference structures available for {formula}")
    reduced_formula = normalize_formula(formula)
    composition = composition_dict_from_formula(reduced_formula)
    atoms_per_formula_unit = float(sum(composition.values()))
    pair_values: dict[str, list[float]] = {
        pair_key: [] for pair_key in required_pair_keys(reduced_formula)
    }
    volume_per_atom_values: list[float] = []
    for atoms in references:
        ref_formula = normalize_formula(atoms.get_chemical_formula())
        if ref_formula != reduced_formula:
            raise ValueError(
                f"Reference formula mismatch: expected {reduced_formula}, got {ref_formula}"
            )
        volume_per_atom_values.append(float(atoms.get_volume()) / max(len(atoms), 1))
        minima = _reference_pair_minima(atoms)
        for pair_key in pair_values:
            if pair_key in minima:
                pair_values[pair_key].append(float(minima[pair_key]))

    missing_pairs = [pair_key for pair_key, values in pair_values.items() if not values]
    if missing_pairs:
        joined_pairs = ", ".join(sorted(missing_pairs))
        raise ValueError(
            f"Reference structures for {reduced_formula} are missing pair minima "
            f"for: {joined_pairs}"
        )
    canonical_minsep = {
        pair_key: float(median(values)) for pair_key, values in pair_values.items()
    }
    volume_per_atom = float(median(volume_per_atom_values))
    return VolumeMinsepEstimate(
        input_formula=formula,
        reduced_formula=reduced_formula,
        composition=composition,
        atoms_per_formula_unit=atoms_per_formula_unit,
        volume_per_atom=volume_per_atom,
        total_volume=volume_per_atom * atoms_per_formula_unit,
        canonical_minsep=canonical_minsep,
        volume_source=EstimateSource(
            method="reference_structure",
            n_observations=len(volume_per_atom_values),
        ),
        minsep_sources={
            pair_key: EstimateSource(
                method="reference_structure", n_observations=len(values)
            )
            for pair_key, values in pair_values.items()
        },
        formula_match_count=len(references),
    )


def resolve_nform(
    *,
    atoms_per_formula_unit: float,
    max_atoms: int,
    max_nform: int = 8,
    nform_override: str | None = None,
) -> str | dict[str, list[int]]:
    """Resolve the buildcell ``#NFORM`` value from atom budget constraints."""
    if nform_override:
        return nform_override
    if max_atoms <= 0:
        raise ValueError("max_atoms must be positive")
    if max_nform <= 0:
        raise ValueError("max_nform must be positive")

    max_formula_units = int(max_atoms // atoms_per_formula_unit)
    if max_formula_units < 1:
        raise ValueError(
            f"Formula unit contains {atoms_per_formula_unit:g} atoms, which exceeds "
            f"max_atoms={max_atoms}"
        )

    allowed_nforms = [2, 3, 4, 6, 8]
    candidates = [nform for nform in allowed_nforms if nform <= max_nform]
    candidates = [nform for nform in candidates if nform <= max_formula_units]
    if not candidates:
        candidates = [1]
    return {"random": candidates}


def apply_minsep_headroom(
    canonical_minsep: dict[str, float],
    *,
    low_scale: float = DEFAULT_MINSEP_SCALE_LOW,
    high_scale: float = DEFAULT_MINSEP_SCALE_HIGH,
    digits: int = DEFAULT_MINSEP_DECIMALS,
) -> dict[str, tuple[float, float]]:
    """Convert canonical minsep values into buildcell ranges."""
    expanded: dict[str, tuple[float, float]] = {}
    for pair_key, value in canonical_minsep.items():
        low = round(float(value) * low_scale, digits)
        high = round(float(value) * high_scale, digits)
        expanded[pair_key] = (low, high)
    return expanded


def build_seed_text_from_estimate(
    seed_text: str,
    *,
    estimate: VolumeMinsepEstimate,
    volume_scale: float = DEFAULT_VOLUME_SCALE,
    minsep_scale_low: float = DEFAULT_MINSEP_SCALE_LOW,
    minsep_scale_high: float = DEFAULT_MINSEP_SCALE_HIGH,
    nform: str | dict[str, list[int]] = "1-2",
) -> str:
    """Inject formula, volume, minsep, and nform directives into seed text."""
    from .search import inject_buildcell_estimate_directives

    varvol = estimate.total_volume * float(volume_scale)
    minsep = apply_minsep_headroom(
        estimate.canonical_minsep,
        low_scale=minsep_scale_low,
        high_scale=minsep_scale_high,
    )
    return inject_buildcell_estimate_directives(
        seed_text,
        formula=estimate.reduced_formula,
        varvol=varvol,
        minsep=minsep,
        nform=nform,
    )


def _reference_pair_minima(
    atoms: Atoms,
    *,
    repeat: tuple[int, int, int] = (2, 2, 2),
) -> dict[str, float]:
    repeated = atoms.repeat(repeat)
    if len(repeated) < 2:
        return {}
    distances = np.asarray(repeated.get_all_distances(mic=False), dtype=float)
    symbols = repeated.get_chemical_symbols()
    minima: dict[str, float] = {}
    for left_index in range(len(repeated)):
        for right_index in range(left_index + 1, len(repeated)):
            distance = float(distances[left_index, right_index])
            if distance <= 1e-8:
                continue
            pair_key = canonical_pair_key(symbols[left_index], symbols[right_index])
            current = minima.get(pair_key)
            if current is None or distance < current:
                minima[pair_key] = distance
    return minima


def validate_estimate_options(
    *,
    volume_scale: float,
    minsep_scale_low: float,
    minsep_scale_high: float,
    max_atoms: int,
    max_nform: int,
) -> None:
    """Validate common CLI estimate options."""
    if not math.isfinite(volume_scale) or volume_scale <= 0:
        raise ValueError("volume_scale must be finite and > 0")
    if not math.isfinite(minsep_scale_low) or minsep_scale_low <= 0:
        raise ValueError("minsep_scale_low must be finite and > 0")
    if not math.isfinite(minsep_scale_high) or minsep_scale_high <= 0:
        raise ValueError("minsep_scale_high must be finite and > 0")
    if minsep_scale_low > minsep_scale_high:
        raise ValueError("minsep_scale_low must be <= minsep_scale_high")
    if max_atoms <= 0:
        raise ValueError("max_atoms must be positive")
    if max_nform <= 0:
        raise ValueError("max_nform must be positive")

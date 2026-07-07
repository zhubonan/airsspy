"""Lightweight baseline predictor for volume and minsep estimates."""

from __future__ import annotations

import hashlib
import json
import math
import warnings
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from pymatgen.core import Composition, Element

from .volume_minsep_data import (
    PairAggregate,
    VolumeAggregate,
    aggregate_training_targets,
    composition_dict_from_formula,
    enumerate_formula_pairs,
    load_dataset_records,
    pair_rows_to_dicts,
    volume_rows_to_dicts,
)

PROPERTY_NAMES = (
    "Z",
    "X",
    "atomic_mass",
    "atomic_radius",
    "atomic_radius_calculated",
    "metallic_radius",
    "average_ionic_radius",
    "mendeleev_no",
    "row",
    "group",
    "min_oxidation_state",
    "max_oxidation_state",
)


def _coerce_float(value: object) -> float:
    if value is None:
        return 0.0
    if isinstance(value, bool):
        return float(value)
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return 0.0
    if not math.isfinite(numeric):
        return 0.0
    return numeric


def _element_property_vector(symbol: str) -> np.ndarray:
    element = Element(symbol)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        values = [_coerce_float(getattr(element, name, None)) for name in PROPERTY_NAMES]
        values.extend(
            (
                float(element.is_metal),
                float(element.is_transition_metal),
                float(element.is_post_transition_metal),
                float(element.is_metalloid),
                float(element.is_alkali),
                float(element.is_alkaline),
                float(element.is_halogen),
                float(element.is_chalcogen),
                float(element.is_lanthanoid),
                float(element.is_actinoid),
            )
        )
    return np.asarray(values, dtype=float)


@dataclass(frozen=True)
class FeatureBuilder:
    """Build reusable formula and pair feature vectors."""

    element_list: tuple[str, ...]

    def __post_init__(self) -> None:
        index = {symbol: position for position, symbol in enumerate(self.element_list)}
        object.__setattr__(self, "element_index", index)
        object.__setattr__(
            self,
            "property_cache",
            {symbol: _element_property_vector(symbol) for symbol in self.element_list},
        )
        object.__setattr__(self, "_formula_cache", {})

    @classmethod
    def from_formulas(cls, formulas: Iterable[str]) -> FeatureBuilder:
        elements = set()
        for formula in formulas:
            elements.update(composition_dict_from_formula(formula))
        sorted_elements = tuple(sorted(elements, key=lambda symbol: Element(symbol).Z))
        return cls(element_list=sorted_elements)

    @property
    def formula_feature_names(self) -> list[str]:
        names = [
            "n_species",
            "n_atoms_reduced",
            "log1p_n_atoms_reduced",
            "composition_entropy",
            "max_fraction",
            "min_fraction_nonzero",
        ]
        names.extend([f"fraction::{symbol}" for symbol in self.element_list])
        stat_names = ("mean", "std", "min", "max", "range")
        base_names = list(PROPERTY_NAMES) + [
            "is_metal",
            "is_transition_metal",
            "is_post_transition_metal",
            "is_metalloid",
            "is_alkali",
            "is_alkaline",
            "is_halogen",
            "is_chalcogen",
            "is_lanthanoid",
            "is_actinoid",
        ]
        for stat_name in stat_names:
            names.extend(
                [f"{stat_name}::{feature_name}" for feature_name in base_names]
            )
        return names

    @property
    def pair_feature_names(self) -> list[str]:
        names = list(self.formula_feature_names)
        names.extend([f"pair_left::{symbol}" for symbol in self.element_list])
        names.extend([f"pair_right::{symbol}" for symbol in self.element_list])
        names.extend(
            [
                "pair_same_element",
                "pair_fraction_left",
                "pair_fraction_right",
                "pair_fraction_sum",
                "pair_fraction_abs_diff",
            ]
        )
        base_names = list(PROPERTY_NAMES) + [
            "is_metal",
            "is_transition_metal",
            "is_post_transition_metal",
            "is_metalloid",
            "is_alkali",
            "is_alkaline",
            "is_halogen",
            "is_chalcogen",
            "is_lanthanoid",
            "is_actinoid",
        ]
        for op_name in ("pair_mean", "pair_abs_diff", "pair_min", "pair_max"):
            names.extend([f"{op_name}::{feature_name}" for feature_name in base_names])
        return names

    def _composition_items(self, formula: str) -> tuple[dict[str, float], float]:
        composition = composition_dict_from_formula(formula)
        total_atoms = float(sum(composition.values()))
        return composition, total_atoms

    def formula_features(self, formula: str) -> np.ndarray:
        cached = self._formula_cache.get(formula)
        if cached is not None:
            return cached

        composition, total_atoms = self._composition_items(formula)
        fractions = np.zeros(len(self.element_list), dtype=float)
        present_symbols: list[str] = []
        present_fractions: list[float] = []
        for symbol, amount in composition.items():
            fraction = float(amount) / total_atoms
            fractions[self.element_index[symbol]] = fraction
            present_symbols.append(symbol)
            present_fractions.append(fraction)

        weight_vector = np.asarray(present_fractions, dtype=float)
        property_matrix = np.vstack(
            [self.property_cache[symbol] for symbol in present_symbols]
        )
        weighted_mean = weight_vector @ property_matrix
        centered = property_matrix - weighted_mean
        weighted_var = weight_vector @ np.square(centered)
        weighted_std = np.sqrt(np.maximum(weighted_var, 0.0))
        weighted_min = property_matrix.min(axis=0)
        weighted_max = property_matrix.max(axis=0)
        weighted_range = weighted_max - weighted_min

        entropy = -sum(
            fraction * math.log(fraction)
            for fraction in present_fractions
            if fraction > 0
        )
        result = np.concatenate(
            [
                np.asarray(
                    [
                        float(len(composition)),
                        total_atoms,
                        math.log1p(total_atoms),
                        entropy,
                        float(max(present_fractions)),
                        float(min(present_fractions)),
                    ],
                    dtype=float,
                ),
                fractions,
                weighted_mean,
                weighted_std,
                weighted_min,
                weighted_max,
                weighted_range,
            ]
        )
        self._formula_cache[formula] = result
        return result

    def pair_features(self, formula: str, pair_key: str) -> np.ndarray:
        formula_features = self.formula_features(formula)
        left_symbol, right_symbol = pair_key.split("-", maxsplit=1)

        left_indicator = np.zeros(len(self.element_list), dtype=float)
        right_indicator = np.zeros(len(self.element_list), dtype=float)
        left_indicator[self.element_index[left_symbol]] = 1.0
        right_indicator[self.element_index[right_symbol]] = 1.0

        composition, total_atoms = self._composition_items(formula)
        left_fraction = float(composition[left_symbol]) / total_atoms
        right_fraction = float(composition[right_symbol]) / total_atoms

        left_properties = self.property_cache[left_symbol]
        right_properties = self.property_cache[right_symbol]
        pair_stats = np.concatenate(
            [
                0.5 * (left_properties + right_properties),
                np.abs(left_properties - right_properties),
                np.minimum(left_properties, right_properties),
                np.maximum(left_properties, right_properties),
            ]
        )
        pair_scalars = np.asarray(
            [
                float(left_symbol == right_symbol),
                left_fraction,
                right_fraction,
                left_fraction + right_fraction,
                abs(left_fraction - right_fraction),
            ],
            dtype=float,
        )
        return np.concatenate(
            [
                formula_features,
                left_indicator,
                right_indicator,
                pair_scalars,
                pair_stats,
            ]
        )


def split_name_for_formula(
    formula: str,
    *,
    seed: int,
    train_ratio: float = 0.8,
    val_ratio: float = 0.1,
) -> str:
    """Assign a deterministic split label to a reduced formula."""
    encoded = f"{seed}:{formula}".encode()
    digest = hashlib.sha256(encoded).hexdigest()
    fraction = int(digest[:12], 16) / float(16**12)
    if fraction < train_ratio:
        return "train"
    if fraction < train_ratio + val_ratio:
        return "val"
    return "test"


def _stack_features(items: Sequence, feature_fn: Callable[[object], np.ndarray]):
    if not items:
        raise ValueError("Cannot stack features for an empty item list")
    return np.vstack([feature_fn(item) for item in items]).astype(float)


def _safe_log(values: np.ndarray) -> np.ndarray:
    if np.any(values <= 0):
        raise ValueError("Targets must be strictly positive before applying log transform")
    return np.log(values)


def _metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    error = y_pred - y_true
    abs_error = np.abs(error)
    return {
        "mae": float(np.mean(abs_error)),
        "rmse": float(np.sqrt(np.mean(np.square(error)))),
        "median_absolute_error": float(np.median(abs_error)),
    }


@dataclass
class BaselineRegressor:
    """Ridge regressor with explicit feature standardization."""

    alpha: float
    transform: str
    feature_mean: np.ndarray
    feature_scale: np.ndarray
    target_mean: float
    coefficients: np.ndarray
    intercept: float
    feature_names: list[str]

    def predict_raw(self, features: np.ndarray) -> np.ndarray:
        scaled = (features - self.feature_mean) / self.feature_scale
        return scaled @ self.coefficients + self.intercept

    def predict(self, features: np.ndarray) -> np.ndarray:
        raw = self.predict_raw(features)
        if self.transform == "log":
            return np.exp(raw)
        return raw

    def to_payload(self) -> dict:
        return {
            "alpha": self.alpha,
            "transform": self.transform,
            "feature_names": self.feature_names,
            "feature_mean": self.feature_mean.tolist(),
            "feature_scale": self.feature_scale.tolist(),
            "target_mean": self.target_mean,
            "coefficients": self.coefficients.tolist(),
            "intercept": self.intercept,
        }

    @classmethod
    def from_payload(cls, payload: dict) -> BaselineRegressor:
        return cls(
            alpha=float(payload["alpha"]),
            transform=str(payload["transform"]),
            feature_mean=np.asarray(payload["feature_mean"], dtype=float),
            feature_scale=np.asarray(payload["feature_scale"], dtype=float),
            target_mean=float(payload["target_mean"]),
            coefficients=np.asarray(payload["coefficients"], dtype=float),
            intercept=float(payload["intercept"]),
            feature_names=list(payload["feature_names"]),
        )


def fit_ridge_regressor(
    features: np.ndarray,
    targets: np.ndarray,
    *,
    alpha: float,
    feature_names: list[str],
    sample_weight: np.ndarray | None = None,
    transform: str = "log",
) -> BaselineRegressor:
    """Fit a ridge regressor in closed form."""
    if transform == "log":
        transformed_targets = _safe_log(targets)
    elif transform == "identity":
        transformed_targets = targets.astype(float, copy=True)
    else:
        raise ValueError(f"Unsupported transform: {transform}")

    if sample_weight is None:
        sample_weight = np.ones(features.shape[0], dtype=float)
    sample_weight = np.asarray(sample_weight, dtype=float)
    if np.any(sample_weight <= 0):
        raise ValueError("Sample weights must be positive")

    normalized_weight = sample_weight / sample_weight.sum()
    feature_mean = normalized_weight @ features
    centered_features = features - feature_mean
    feature_var = normalized_weight @ np.square(centered_features)
    feature_scale = np.sqrt(np.maximum(feature_var, 1e-12))
    scaled_features = centered_features / feature_scale

    target_mean = float(normalized_weight @ transformed_targets)
    centered_targets = transformed_targets - target_mean

    sqrt_weight = np.sqrt(sample_weight)[:, None]
    lhs = (scaled_features * sqrt_weight).T @ (scaled_features * sqrt_weight)
    rhs = (scaled_features * sqrt_weight).T @ (centered_targets * sqrt_weight[:, 0])
    lhs += alpha * np.eye(lhs.shape[0], dtype=float)
    coefficients = np.linalg.solve(lhs, rhs)
    intercept = target_mean

    return BaselineRegressor(
        alpha=alpha,
        transform=transform,
        feature_mean=feature_mean,
        feature_scale=feature_scale,
        target_mean=target_mean,
        coefficients=coefficients,
        intercept=intercept,
        feature_names=feature_names,
    )


def evaluate_regressor(
    model: BaselineRegressor,
    features: np.ndarray,
    targets: np.ndarray,
) -> dict[str, float]:
    """Return simple regression metrics."""
    return _metrics(targets, model.predict(features))


def _subsample(items: Sequence, limit: int | None) -> list:
    if limit is None or limit >= len(items):
        return list(items)
    return list(items[:limit])


def _rebalance_small_splits(split_items: dict[str, list]) -> None:
    if not split_items["train"]:
        raise ValueError("Training split is empty; increase the dataset size")
    for split in ("val", "test"):
        if not split_items[split] and len(split_items["train"]) > 1:
            split_items[split].append(split_items["train"].pop())


def train_baseline_regressor(
    items: Sequence[VolumeAggregate] | Sequence[PairAggregate],
    *,
    feature_fn: Callable[[object], np.ndarray],
    target_fn: Callable[[object], float],
    weight_fn: Callable[[object], float],
    feature_names: list[str],
    seed: int,
    alpha_grid: Sequence[float],
    max_rows: int | None = None,
) -> tuple[BaselineRegressor, dict]:
    """Train a deterministic ridge model with validation alpha selection."""
    selected_items = _subsample(items, max_rows)
    split_items = {"train": [], "val": [], "test": []}
    for item in selected_items:
        split_items[split_name_for_formula(item.reduced_formula, seed=seed)].append(
            item
        )
    _rebalance_small_splits(split_items)

    feature_cache = {
        split: _stack_features(group, feature_fn) if group else None
        for split, group in split_items.items()
    }
    target_cache = {
        split: np.asarray([target_fn(item) for item in group], dtype=float)
        for split, group in split_items.items()
    }
    weight_cache = {
        split: np.asarray([weight_fn(item) for item in group], dtype=float)
        for split, group in split_items.items()
    }

    best_alpha = None
    best_model = None
    best_val_rmse = math.inf
    for alpha in alpha_grid:
        model = fit_ridge_regressor(
            feature_cache["train"],
            target_cache["train"],
            alpha=float(alpha),
            feature_names=feature_names,
            sample_weight=weight_cache["train"],
        )
        val_metrics = (
            evaluate_regressor(model, feature_cache["val"], target_cache["val"])
            if split_items["val"]
            else evaluate_regressor(model, feature_cache["train"], target_cache["train"])
        )
        if val_metrics["rmse"] < best_val_rmse:
            best_alpha = float(alpha)
            best_model = model
            best_val_rmse = val_metrics["rmse"]

    assert best_model is not None
    metrics = {
        "selected_alpha": best_alpha,
        "n_rows": len(selected_items),
        "split_counts": {split: len(group) for split, group in split_items.items()},
        "train_metrics": evaluate_regressor(
            best_model, feature_cache["train"], target_cache["train"]
        ),
        "val_metrics": evaluate_regressor(
            best_model, feature_cache["val"], target_cache["val"]
        )
        if split_items["val"]
        else None,
        "test_metrics": evaluate_regressor(
            best_model, feature_cache["test"], target_cache["test"]
        )
        if split_items["test"]
        else None,
    }
    return best_model, metrics


def save_baseline_bundle(
    output_dir: str | Path,
    *,
    builder: FeatureBuilder,
    volume_model: BaselineRegressor,
    minsep_model: BaselineRegressor,
    metrics: dict,
    dataset_summary: dict,
) -> None:
    """Write a JSON baseline bundle."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    payload = {
        "builder": {"element_list": list(builder.element_list)},
        "volume_model": volume_model.to_payload(),
        "minsep_model": minsep_model.to_payload(),
        "metrics": metrics,
        "dataset_summary": dataset_summary,
    }
    with (output_path / "baseline_bundle.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def load_baseline_bundle(
    path: str | Path,
) -> tuple[FeatureBuilder, BaselineRegressor, BaselineRegressor, dict]:
    """Load a JSON baseline bundle."""
    with Path(path).open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    builder = FeatureBuilder(tuple(payload["builder"]["element_list"]))
    volume_model = BaselineRegressor.from_payload(payload["volume_model"])
    minsep_model = BaselineRegressor.from_payload(payload["minsep_model"])
    metadata = {
        "metrics": payload["metrics"],
        "dataset_summary": payload["dataset_summary"],
    }
    return builder, volume_model, minsep_model, metadata


def _compact_summary(summary: dict) -> dict:
    compact = dict(summary)
    compact.pop("element_list", None)
    compact.pop("pair_key_list", None)
    return compact


def train_baseline_models(
    *,
    dataset_path: str | Path,
    output_dir: str | Path,
    seed: int = 17,
    pair_min_count: int = 1,
    volume_alpha_grid: Sequence[float] = (1e-6, 1e-4, 1e-2, 1.0, 1e2),
    minsep_alpha_grid: Sequence[float] = (1e-6, 1e-4, 1e-2, 1.0, 1e2),
    max_volume_rows: int | None = None,
    max_pair_rows: int | None = None,
) -> dict:
    """Train baseline volume and minsep regressors from a dataset JSON file."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    records = load_dataset_records(dataset_path)
    volume_rows, pair_rows, dataset_summary = aggregate_training_targets(records)
    pair_rows = [row for row in pair_rows if row.n_observations >= pair_min_count]

    builder = FeatureBuilder.from_formulas([row.reduced_formula for row in volume_rows])
    volume_model, volume_metrics = train_baseline_regressor(
        volume_rows,
        feature_fn=lambda row: builder.formula_features(row.reduced_formula),
        target_fn=lambda row: row.volume_per_atom,
        weight_fn=lambda row: float(row.n_structures),
        feature_names=builder.formula_feature_names,
        seed=seed,
        alpha_grid=volume_alpha_grid,
        max_rows=max_volume_rows,
    )
    minsep_model, minsep_metrics = train_baseline_regressor(
        pair_rows,
        feature_fn=lambda row: builder.pair_features(row.reduced_formula, row.pair_key),
        target_fn=lambda row: row.minsep,
        weight_fn=lambda row: float(row.n_observations),
        feature_names=builder.pair_feature_names,
        seed=seed,
        alpha_grid=minsep_alpha_grid,
        max_rows=max_pair_rows,
    )

    metrics = {"baseline": {"volume": volume_metrics, "minsep": minsep_metrics}}
    dataset_summary = _compact_summary(dataset_summary)
    dataset_summary["pair_rows_after_min_count_filter"] = len(pair_rows)
    dataset_summary["pair_min_count"] = pair_min_count

    save_baseline_bundle(
        output_path / "baseline",
        builder=builder,
        volume_model=volume_model,
        minsep_model=minsep_model,
        metrics=metrics,
        dataset_summary=dataset_summary,
    )

    with (output_path / "aggregated_volume_targets.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(volume_rows_to_dicts(volume_rows), handle, indent=2)
    with (output_path / "aggregated_minsep_targets.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(pair_rows_to_dicts(pair_rows), handle, indent=2)
    with (output_path / "training_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {"dataset_summary": dataset_summary, "metrics": metrics},
            handle,
            indent=2,
            sort_keys=True,
        )

    return {
        "dataset_summary": dataset_summary,
        "metrics": metrics,
        "output_dir": str(output_path),
    }


def _validate_supported_formula(formula: str, builder: FeatureBuilder) -> dict[str, float]:
    composition = composition_dict_from_formula(formula)
    missing = sorted(set(composition) - set(builder.element_index))
    if missing:
        supported = ", ".join(builder.element_list)
        missing_text = ", ".join(missing)
        raise ValueError(
            f"Formula {formula} contains unsupported elements: {missing_text}. "
            f"Supported elements: {supported}"
        )
    return composition


def _prediction_payload(
    *,
    formula: str,
    composition: dict[str, float],
    volume_per_atom: float,
    canonical_minsep: dict[str, float],
    model_metrics: dict,
    model_artifact: str,
) -> dict:
    atoms_per_formula_unit = float(sum(composition.values()))
    return {
        "reduced_formula": formula,
        "composition": composition,
        "volume_per_atom": volume_per_atom,
        "volume_per_formula_unit": volume_per_atom * atoms_per_formula_unit,
        "canonical_minsep": canonical_minsep,
        "predictor_type": "baseline",
        "model_artifact": model_artifact,
        "model_metrics": model_metrics,
    }


@dataclass
class BaselineFormulaPredictor:
    """Predict formula volume and minsep values from a baseline bundle."""

    bundle_path: str | Path

    def __post_init__(self) -> None:
        self.bundle_path = Path(self.bundle_path)
        (
            self.builder,
            self.volume_model,
            self.minsep_model,
            self.metadata,
        ) = load_baseline_bundle(self.bundle_path)

    def predict(self, formula: str) -> dict:
        reduced_formula = Composition(formula).reduced_formula
        composition = _validate_supported_formula(reduced_formula, self.builder)
        volume_features = self.builder.formula_features(reduced_formula)[None, :]
        volume_per_atom = float(self.volume_model.predict(volume_features)[0])

        canonical_minsep = {}
        for left_symbol, right_symbol in enumerate_formula_pairs(reduced_formula):
            pair_key = f"{left_symbol}-{right_symbol}"
            pair_features = self.builder.pair_features(reduced_formula, pair_key)[None, :]
            canonical_minsep[pair_key] = float(self.minsep_model.predict(pair_features)[0])

        return _prediction_payload(
            formula=reduced_formula,
            composition=composition,
            volume_per_atom=volume_per_atom,
            canonical_minsep=canonical_minsep,
            model_metrics=self.metadata["metrics"],
            model_artifact=str(self.bundle_path),
        )

    def predict_many(self, formulas: Sequence[str]) -> list[dict]:
        """Predict many formulas in order."""
        return [self.predict(formula) for formula in formulas]

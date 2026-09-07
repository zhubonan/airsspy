from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
from ase import Atoms
from pymatgen.core import Lattice, Structure

from airsspy.volume_minsep import (
    build_seed_text_from_estimate,
    lookup_exact_volume_minsep_estimate,
    reference_volume_minsep_estimate,
    resolve_nform,
)
from airsspy.volume_minsep_data import (
    aggregate_training_targets,
    build_dataset_rows,
    curate_minsep_volume_dataset,
    get_minsep_dict,
)
from airsspy.volume_minsep_model import (
    BaselineFormulaPredictor,
    train_baseline_models,
)


def _record(
    *,
    material_id: str,
    reduced_formula: str,
    composition: dict[str, float],
    volume: float,
    minsep: dict[str, float],
    energy_above_hull: float | None = None,
) -> dict:
    record = {
        "material_id": material_id,
        "reduced_formula": reduced_formula,
        "chemical_formula": reduced_formula,
        "composition": dict(composition),
        "volume": volume,
        "minsep": dict(minsep),
    }
    if energy_above_hull is not None:
        record["energy_above_hull"] = energy_above_hull
    return record


def _toy_records() -> list[dict]:
    return [
        _record(
            material_id="toy-1",
            reduced_formula="SiO2",
            composition={"Si": 1.0, "O": 2.0},
            volume=45.0,
            minsep={"Si-O": 1.60, "O-Si": 1.60, "O-O": 2.55, "Si-Si": 3.05},
        ),
        _record(
            material_id="toy-2",
            reduced_formula="SiO2",
            composition={"Si": 2.0, "O": 4.0},
            volume=91.2,
            minsep={"Si-O": 1.58, "O-Si": 1.58, "O-O": 2.50, "Si-Si": 3.10},
        ),
        _record(
            material_id="toy-3",
            reduced_formula="MgO",
            composition={"Mg": 1.0, "O": 1.0},
            volume=22.4,
            minsep={"Mg-O": 1.95, "O-Mg": 1.95, "Mg-Mg": 3.02, "O-O": 3.02},
        ),
        _record(
            material_id="toy-4",
            reduced_formula="NaCl",
            composition={"Na": 1.0, "Cl": 1.0},
            volume=34.0,
            minsep={
                "Na-Cl": 2.81,
                "Cl-Na": 2.81,
                "Na-Na": 3.95,
                "Cl-Cl": 3.95,
            },
        ),
    ]


def _write_toy_dataset(path: Path) -> None:
    path.write_text(json.dumps(_toy_records(), indent=2), encoding="utf-8")


def test_dataset_builder_extracts_minsep_and_curation_selects_lowest_ehull(
    tmp_path: Path,
) -> None:
    structure = Structure(
        lattice=Lattice.cubic(5.64),
        species=["Na", "Cl"],
        coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    rows = build_dataset_rows(
        [{"material_id": "mp-test", "structure": structure}],
        show_progress=False,
    )

    assert rows[0]["material_id"] == "mp-test"
    assert rows[0]["reduced_formula"] == "NaCl"
    assert rows[0]["volume"] == pytest.approx(structure.volume)
    assert get_minsep_dict(structure)["Cl-Na"] > 0

    raw_path = tmp_path / "raw.json"
    mp_docs_path = tmp_path / "mp_docs.json"
    curated_path = tmp_path / "curated.json"
    raw_path.write_text(
        json.dumps(
            [
                _record(
                    material_id="mp-high",
                    reduced_formula="SiO2",
                    composition={"Si": 1.0, "O": 2.0},
                    volume=60.0,
                    minsep={"Si-O": 1.6, "O-O": 2.4, "Si-Si": 3.0},
                ),
                _record(
                    material_id="mp-low",
                    reduced_formula="SiO2",
                    composition={"Si": 1.0, "O": 2.0},
                    volume=45.0,
                    minsep={"Si-O": 1.5, "O-O": 2.3, "Si-Si": 2.9},
                ),
            ]
        ),
        encoding="utf-8",
    )
    mp_docs_path.write_text(
        json.dumps(
            [
                {"material_id": "mp-high", "energy_above_hull": 0.05},
                {"material_id": "mp-low", "energy_above_hull": 0.0},
            ]
        ),
        encoding="utf-8",
    )

    summary = curate_minsep_volume_dataset(
        minsep_dataset_path=raw_path,
        mp_docs_path=mp_docs_path,
        output_path=curated_path,
    )

    curated = json.loads(curated_path.read_text(encoding="utf-8"))
    assert summary["curated_rows"] == 1
    assert curated[0]["material_id"] == "mp-low"
    assert curated[0]["energy_above_hull"] == 0.0
    assert curated_path.with_suffix(".summary.json").exists()


def test_aggregation_training_and_baseline_prediction_do_not_import_torch(
    tmp_path: Path,
) -> None:
    had_torch = "torch" in sys.modules
    dataset_path = tmp_path / "toy_dataset.json"
    output_dir = tmp_path / "artifacts"
    _write_toy_dataset(dataset_path)

    records = json.loads(dataset_path.read_text(encoding="utf-8"))
    volume_rows, pair_rows, summary = aggregate_training_targets(records)

    si_row = next(row for row in volume_rows if row.reduced_formula == "SiO2")
    assert si_row.volume_per_atom == pytest.approx(15.1)
    assert {row.pair_key for row in pair_rows if row.reduced_formula == "SiO2"} == {
        "O-O",
        "O-Si",
        "Si-Si",
    }
    assert summary["n_unique_formulas"] == 3

    result = train_baseline_models(
        dataset_path=dataset_path,
        output_dir=output_dir,
        seed=7,
        pair_min_count=1,
        volume_alpha_grid=[1e-4, 1e-2, 1.0],
        minsep_alpha_grid=[1e-4, 1e-2, 1.0],
    )

    bundle_path = output_dir / "baseline" / "baseline_bundle.json"
    assert bundle_path.exists()
    assert (output_dir / "aggregated_volume_targets.json").exists()
    assert (output_dir / "aggregated_minsep_targets.json").exists()
    assert result["dataset_summary"]["pair_rows_after_min_count_filter"] > 0

    prediction = BaselineFormulaPredictor(bundle_path=bundle_path).predict("SiO2")
    assert prediction["composition"] == {"Si": 1.0, "O": 2.0}
    assert set(prediction["canonical_minsep"]) == {"O-O", "O-Si", "Si-Si"}
    assert prediction["volume_per_atom"] > 0
    if not had_torch:
        assert "torch" not in sys.modules


def test_exact_and_reference_estimates_build_seed_text(tmp_path: Path) -> None:
    dataset_path = tmp_path / "curated.json"
    dataset_path.write_text(
        json.dumps(
            [
                _record(
                    material_id="mp-sio2",
                    reduced_formula="SiO2",
                    composition={"Si": 1.0, "O": 2.0},
                    volume=45.0,
                    minsep={"Si-O": 1.6, "O-O": 2.5, "Si-Si": 3.0},
                    energy_above_hull=0.0,
                )
            ]
        ),
        encoding="utf-8",
    )

    estimate = lookup_exact_volume_minsep_estimate(
        formula="O2Si",
        dataset_path=dataset_path,
    )
    assert estimate.reduced_formula == "SiO2"
    assert estimate.volume_per_atom == pytest.approx(15.0)
    assert estimate.total_volume == pytest.approx(45.0)
    assert resolve_nform(
        atoms_per_formula_unit=estimate.atoms_per_formula_unit,
        max_atoms=12,
        max_nform=3,
    ) == {"random": [2, 3]}
    assert resolve_nform(
        atoms_per_formula_unit=estimate.atoms_per_formula_unit,
        max_atoms=5,
        max_nform=3,
    ) == {"random": [1]}

    out = build_seed_text_from_estimate(
        "#FORMULA=Si\n#VARVOL=999\n#MINSEP=9\n#NFORM=1\n#SLACK=0.25\n",
        estimate=estimate,
        volume_scale=1.1,
        minsep_scale_low=0.9,
        minsep_scale_high=1.1,
        nform={"random": [2, 3]},
    )

    assert "#FORMULA=SiO2" in out
    assert "#VARVOL=49.5" in out
    assert "#MINSEP=0.5-1 O-O=2.25-2.75 O-Si=1.44-1.76 Si-Si=2.7-3.3" in out
    assert "#NFORM={2,3}" in out
    assert "#SLACK=0.25" in out
    assert "#VARVOL=999" not in out

    atoms = Atoms(
        symbols=["Na", "Cl"],
        positions=[[0.0, 0.0, 0.0], [2.82, 2.82, 2.82]],
        cell=[5.64, 5.64, 5.64],
        pbc=True,
    )
    reference = reference_volume_minsep_estimate(formula="NaCl", references=[atoms])
    assert reference.volume_source.method == "reference_structure"
    assert reference.canonical_minsep["Cl-Na"] > 0

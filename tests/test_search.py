import random

import pytest
from ase import Atoms

from airsspy.search import (
    FormulaSamplingOptions,
    RssCandidate,
    RssPruneOptions,
    build_formula_sampling_context,
    inject_formula_directive,
    make_seed_text_transform,
    pool_statistics,
    prune_relaxed_pool,
    select_pruned_candidates,
    should_flush_prune_pool,
    validate_prune_options,
)


def test_inject_formula_directive_removes_conflicts():
    seed = "\n".join(
        [
            "#SPECIES=Si,O",
            "#NATOM=4-16",
            "#FORMULA=SiO",
            "#VARVOL=999",
            "#SLACK=0.25",
        ]
    )

    out = inject_formula_directive(seed, "SiO2", varvol=20.0)

    assert "#FORMULA=SiO2" in out
    assert "#VARVOL=20" in out
    assert "#SPECIES=" not in out
    assert "#NATOM=" not in out
    assert "#FORMULA=SiO\n" not in out
    assert "#SLACK=0.25" in out


def test_formula_context_canonicalizes_and_samples_deterministically():
    context = build_formula_sampling_context(
        FormulaSamplingOptions(
            formulas=["O2Si", "SiO2", "Li2TiO3"],
            target_atom_volumes={"Si": 10.0, "O": 8.0, "Li": 12.0, "Ti": 11.0},
        )
    )

    assert context.formulas == ["Li2TiO3", "SiO2"]
    assert context.varvol_by_formula["SiO2"] == pytest.approx(52.0 / 3.0)
    assert context.varvol_by_formula["Li2TiO3"] == pytest.approx(59.0 / 6.0 * 3.0)

    transform = make_seed_text_transform(context, rng=random.Random(1))
    out = transform("#SPECIES=Li,Ti,O\n#NFORM=1\n")

    assert "#FORMULA=" in out
    assert "#SPECIES=" not in out


def test_target_volume_scales_formula_average_by_species_types():
    sio2_context = build_formula_sampling_context(
        FormulaSamplingOptions(
            formulas=["SiO2"],
            target_atom_volumes={"Si": 10.0, "O": 6.0},
        )
    )
    srtio3_context = build_formula_sampling_context(
        FormulaSamplingOptions(
            formulas=["SrTiO3"],
            target_atom_volumes={"Sr": 20.0, "Ti": 15.0, "O": 10.0},
        )
    )

    assert sio2_context.varvol_by_formula["SiO2"] == pytest.approx(22.0 / 3.0 * 2.0)
    assert srtio3_context.varvol_by_formula["SrTiO3"] == pytest.approx(65.0 / 5.0 * 3.0)


def test_formula_context_filters_by_seed_constraints_and_oxidation_states():
    context = build_formula_sampling_context(
        FormulaSamplingOptions(
            elements=["Li", "Ti", "O"],
            max_coeff=3,
            oxidation_states={"Li": [1], "Ti": [4], "O": [-2]},
        ),
        seed_text="#NATOM=4-24\n#NFORM=1-2\n",
    )

    assert "Li2TiO3" in context.formulas
    for formula in context.formulas:
        assert "O" in formula


def test_formula_context_errors_on_missing_target_volume():
    with pytest.raises(ValueError, match="Missing target atom volumes"):
        build_formula_sampling_context(
            FormulaSamplingOptions(
                formulas=["SiO2"],
                target_atom_volumes={"Si": 10.0},
            )
        )


def _candidate(label, distance, energy):
    atoms = Atoms(
        "H2",
        positions=[[0.0, 0.0, 0.0], [distance, 0.0, 0.0]],
        cell=[10.0, 10.0, 10.0],
        pbc=True,
    )
    return RssCandidate(label=label, atoms=atoms, energy=energy)


def test_prune_selects_low_energy_unique_candidates():
    options = RssPruneOptions(
        enabled=True,
        pool_size=4,
        keep_fraction=0.5,
        dedup_tol=1e-6,
        fingerprint_cutoff=5.0,
    )
    candidates = [
        _candidate("low-a", 1.0, -4.0),
        _candidate("low-duplicate", 1.0, -3.0),
        _candidate("next-unique", 2.0, -2.0),
        _candidate("high-unique", 3.0, -1.0),
    ]

    kept, rejected = select_pruned_candidates(candidates, options)

    assert [candidate.label for candidate in kept] == ["low-a", "next-unique"]
    assert "low-duplicate" in {candidate.label for candidate in rejected}


def test_prune_relaxed_pool_groups_by_formula():
    options = RssPruneOptions(enabled=True, keep_fraction=1.0)
    h2 = _candidate("h2", 1.0, -2.0)
    he = RssCandidate(
        label="he",
        atoms=Atoms("He", positions=[[0.0, 0.0, 0.0]], cell=[8.0, 8.0, 8.0], pbc=True),
        energy=-1.0,
    )

    kept = prune_relaxed_pool([h2, he], options)

    assert {candidate.label for candidate in kept} == {"h2", "he"}


def test_prune_flush_triggers_and_validation():
    options = RssPruneOptions(
        enabled=True,
        pool_size=3,
        min_stable_pool_size=2,
        stable_window=1,
        mean_abs_tol=1.0,
        median_abs_tol=1.0,
    )
    candidates = [
        _candidate("a", 1.0, -2.0),
        _candidate("b", 2.0, -2.1),
        _candidate("c", 3.0, -2.2),
    ]
    stats = [pool_statistics(candidates[: i + 1]) for i in range(len(candidates))]

    assert should_flush_prune_pool(candidates, stats, options) == ("pool_size", True)

    with pytest.raises(ValueError, match="pool_size"):
        validate_prune_options(RssPruneOptions(pool_size=0))

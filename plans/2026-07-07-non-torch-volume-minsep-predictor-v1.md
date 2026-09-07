# Non-Torch Volume and Minsep Predictor Integration

## Objective

Integrate the non-torch volume and minsep prediction functionality from the sibling
`airss-arena` project directly into `airsspy`, without adding `airss-arena` as a
runtime dependency and without enabling the torch/deep predictor path in this
first slice. The outcome should let users derive buildcell `#VARVOL`, pairwise
`#MINSEP`, and automatic `#NFORM` settings from a small/simple baseline
volume/minsep predictor, an exact user-supplied curated dataset row, or reference
structures, then apply those settings through existing `ap run search` formula
sampling.

This plan assumes the first implementation keeps Python 3.9 compatibility and
does not add torch as a required dependency. It may use dependencies already in
`airsspy`, especially `numpy`, `pymatgen`, and `ase` from `pyproject.toml:32`.

## Implementation Plan

- [ ] 1. Add package-local runtime resources for non-torch volume/minsep estimates.
  Copy only small runtime artifacts into a new `src/airsspy/data/` tree. The
  trained baseline bundle in arena is small enough to consider packaging, but
  the curated minsep/volume dataset should not be bundled by default because
  `../airss-arena/src/airss_arena/data/minsep_vol_dataset_curated.json` is a Git
  LFS object with recorded size `62112633` bytes. Treat curated datasets as
  user-supplied or generated training inputs instead. Arena packages resources
  under `src/airss_arena/data/` and resolves them through
  `src/airss_arena/_resources.py:9`; add corresponding package-data entries to
  `pyproject.toml:78` only for the small bundle and any lightweight metadata.
  Keep deep `.pt` artifacts out of scope for this slice.

- [ ] 2. Introduce an `airsspy` resource helper module.
  Add a small module such as `src/airsspy/resources.py` with helpers for
  locating `data/`, the default baseline bundle, and optional user-provided
  curated dataset paths. This should mirror the useful shape of arena's resource helpers in
  `../airss-arena/src/airss_arena/_resources.py:9`, but use `airsspy` package
  paths and avoid arena imports. Tests should verify the packaged paths exist
  when installed from the local source tree.

- [ ] 3. Port formula/pair data utilities into `airsspy`.
  Add a focused module such as `src/airsspy/volume_minsep_data.py` containing canonical
  pair handling, formula composition parsing, required pair enumeration, dataset
  JSON loading, volume-per-atom extraction, and canonical minsep-map
  normalization. The source behavior is in
  `../airss-arena/src/airss_arena/formula_model/data.py:13`, including canonical pair
  keys, composition dictionaries, formula pair enumeration, and canonicalizing
  ordered minsep entries to the shortest unordered pair value.

- [ ] 4. Port dataset generation and curation as offline tooling.
  Add developer/user-facing utilities for creating the large curated training
  dataset rather than copying that dataset into the package. The raw dataset
  builder is `../airss-arena/src/airss_arena/formula_model/dataset_builder.py:15`
  and converts Materials Project document dumps into rows containing `minsep`,
  `material_id`, `volume`, `composition`, `reduced_formula`, and
  `chemical_formula`; its entry point is
  `../airss-arena/src/airss_arena/formula_model/dataset_builder.py:52`. The
  curation script is
  `../airss-arena/src/airss_arena/formula_model/curate_dataset.py:16` and keeps
  one row per reduced formula by selecting the lowest `energy_above_hull` entry
  from the MP docs. Expose this in `airsspy` with neutral volume/minsep naming,
  not arena's broad "prior" wording.

- [ ] 5. Port the lightweight baseline predictor only.
  Add a non-torch predictor module, for example
  `src/airsspy/volume_minsep_model.py`, containing the deterministic feature builder and
  JSON-loaded ridge baseline regressor. Use the arena baseline implementation as
  the behavioral reference:
  `../airss-arena/src/airss_arena/formula_model/features.py:68` and
  `../airss-arena/src/airss_arena/formula_model/baseline.py:62`. Do not port
  arena's deep predictor or torch model classes in this slice.

- [ ] 6. Add a lightweight baseline training module.
  Add a baseline-only training module or CLI command that consumes a raw or
  curated minsep/volume dataset and writes a small `baseline_bundle.json`.
  Reuse the non-deep flow from
  `../airss-arena/src/airss_arena/formula_model/train.py:54`: load dataset
  records, aggregate formula and pair targets via
  `../airss-arena/src/airss_arena/formula_model/data.py:120`, create a
  `FeatureBuilder`, train separate ridge regressors for volume-per-atom and
  minsep, save aggregated target summaries, and write training metrics. Do not
  migrate arena's `train.py` wholesale because it imports
  `../airss-arena/src/airss_arena/formula_model/deep.py` at module import time
  even when `train_deep` is false. The `airsspy` trainer should import only
  numpy/pymatgen-based modules and should not import torch.

- [ ] 7. Add the central volume/minsep estimate API.
  Add a module such as `src/airsspy/volume_minsep.py` with `EstimateSource`,
  `VolumeMinsepEstimate`, `normalize_formula`, `flatten_formula`,
  `lookup_exact_volume_minsep_estimate`, `volume_minsep_estimate_from_prediction`,
  `predict_baseline_volume_minsep_estimate`, `reference_volume_minsep_estimate`,
  `resolve_nform`, `apply_minsep_headroom`, and `build_seed_atoms_from_estimate`.
  The API should be modeled on arena's estimate-producing flow in
  `../airss-arena/src/airss_arena/generation_backend.py:91`,
  `../airss-arena/src/airss_arena/generation_backend.py:206`,
  `../airss-arena/src/airss_arena/generation_backend.py:253`,
  `../airss-arena/src/airss_arena/generation_backend.py:297`, and
  `../airss-arena/src/airss_arena/generation_backend.py:338`. Map estimates into
  existing `SeedAtoms` tags rather than inventing a new serializer, because
  `BuildcellParam` already supports `nform`, nested `minsep`, `varvol`, and
  `seed` in `src/airsspy/seed.py:493`.

- [ ] 8. Preserve existing manual formula sampling behavior.
  Extend, rather than replace, `FormulaSamplingOptions` and
  `FormulaSamplingContext` in `src/airsspy/search.py:35`. Existing users can
  already inject `#FORMULA` and optional manual `#VARVOL` through
  `inject_formula_directive()` at `src/airsspy/search.py:131`, and tests cover
  this behavior in `tests/test_search.py:21`. The new estimator path should be
  optional and should not change output when no volume/minsep estimator is requested.

- [ ] 9. Extend seed text transformation to inject minsep and nform.
  Add a structured seed transform path that can remove conflicting directives
  and prepend `#FORMULA`, `#VARVOL`, `#MINSEP`, and `#NFORM` from resolved
  estimates. Existing `DEFAULT_FORMULA_REMOVE_DIRECTIVES` in
  `src/airsspy/search.py:26` removes `NATOM`, `SPECIES`, `FORMULA`, `VARVOL`,
  and `TARGVOL`; extend or parameterize this for `MINSEP` and `NFORM` only when
  auto-estimates are active. Keep manual `#SLACK`, `#OVERLAP`, and other user seed
  directives unless explicitly superseded. Verify the nested minsep string
  matches current serializer behavior covered by `tests/test_seed.py:52`.

- [ ] 10. Add CLI controls for non-torch volume/minsep estimator selection.
  Extend `ap run search` in `src/airsspy/cli/cmd_run.py:927` with options for
  selecting an estimate source and tuning headroom. Suggested flags:
  `--volume-minsep-source {none,dataset,baseline,reference}`,
  `--volume-minsep-dataset`, `--volume-minsep-bundle`, `--volume-scale`,
  `--minsep-scale-low`,
  `--minsep-scale-high`, `--max-atoms`, `--max-nform`, and an optional
  `--reference-structure` repeatable path for reference-derived estimates. The
  first implementation should require `--formula` or resolved formula sampling
  for estimate lookup; formula enumeration through `--elements` can be supported by
  resolving estimates for each generated formula after filtering.

- [ ] 11. Add CLI controls for dataset build, curation, and baseline training.
  Add a small `ap` command group or subcommands with neutral naming, for example
  volume/minsep dataset build, dataset curate, and baseline train. Arena's CLI
  currently exposes these as `tools build-prior-dataset`,
  `tools curate-prior-dataset`, and `tools train-prior-model` in
  `../airss-arena/src/airss_arena/cli.py:627`,
  `../airss-arena/src/airss_arena/cli.py:635`, and
  `../airss-arena/src/airss_arena/cli.py:649`. The `airsspy` commands should
  avoid "prior" terminology, should clearly mark MP-docs input as user supplied,
  and should support baseline-only training without any deep-model flags.

- [ ] 12. Improve `--diagnose` output for estimator-driven seeds.
  Extend the existing diagnostic branch in `src/airsspy/cli/cmd_run.py:1204` to
  print estimate provenance, resolved volume per atom, total volume, scaled
  `#VARVOL`, canonical minsep values, minsep ranges, and resolved `#NFORM`.
  This creates a safe preview mode for buildcell seed changes before users run a
  long search. Keep output human-readable and avoid requiring buildcell.

- [ ] 13. Add unit tests for dataset build, curation, and training.
  Add small synthetic tests that exercise `get_minsep_dict`-equivalent behavior,
  raw dataset row construction from pymatgen structures, curation by
  `energy_above_hull`, and baseline training on a toy dataset. Arena's smoke
  coverage is in `../airss-arena/tests/test_formula_model.py:77` and
  `../airss-arena/tests/test_formula_model.py:108`. The tests should confirm the
  trainer writes `baseline_bundle.json`, aggregated target JSON, and a training
  summary, all without importing torch.

- [ ] 14. Add unit tests for data and baseline prediction.
  Add tests covering canonical pair normalization, exact dataset lookup, missing
  pairs, baseline bundle loading, unsupported element errors, and successful
  volume/minsep prediction for formulas present in the bundled model. Arena's comparable
  tests live in `../airss-arena/tests/test_formula_model.py:92` and
  `../airss-arena/tests/test_formula_model.py:108`. These tests should not
  import `airss_arena` and should not require torch.

- [ ] 15. Add unit tests for estimate-to-seed mapping.
  Add tests that construct a small `VolumeMinsepEstimate` object and verify the
  generated `SeedAtoms` or seed text includes expected `#VARVOL`, pairwise
  `#MINSEP`, randomized `#NFORM`, and stable atom tag names. Existing seed
  serialization tests in `tests/test_seed.py:24` and `tests/test_seed.py:166`
  provide the local style to extend.

- [ ] 16. Add CLI tests for estimator-driven formula sampling.
  Extend CLI coverage around `ap run search --diagnose` so tests can validate
  the resolved buildcell text without invoking buildcell or external DFT codes.
  The command already supports a diagnose-only path at
  `src/airsspy/cli/cmd_run.py:1204`, making it the safest integration test
  point. Cover dataset/baseline source selection, conflicting directive removal,
  and preservation of unrelated seed settings.

- [ ] 17. Document the new non-torch volume/minsep estimator workflow.
  Update CLI/reference documentation, likely `docs/reference/cli.md`, and add a
  short how-to under `docs/how-to/` describing automatic volume and minsep
  configuration. Explain the distinction between manual `--target-volume` from
  `src/airsspy/cli/cmd_run.py:1019` and auto estimates, and state that deep/torch
  prediction is intentionally out of scope for this first migration. Include the
  offline workflow for generating a raw minsep/volume dataset from MP docs,
  curating one row per formula, and training the small baseline bundle.

- [ ] 18. Run focused verification.
  Run the new tests plus existing related coverage: `tests/test_search.py`,
  `tests/test_seed.py`, and CLI tests around `run search`. Run `ruff check` on
  touched modules. If package-data changes are made, verify the resource helpers
  can load the packaged JSON from the source checkout. Do not run buildcell
  end-to-end unless the executable is available; buildcell-dependent tests
  should remain guarded as in `tests/test_build.py:45`.

- [ ] 19. Request read-only review before finalizing implementation.
  After implementation, spawn the required read-only reviewer subagent from
  `AGENTS.md:214`. Ask it to inspect the branch diff for correctness,
  regressions, missing tests, race conditions, security issues, and unrelated
  edits. Fix high-confidence issues before final response.

## Verification Criteria

- The package can resolve a bundled baseline bundle path and user-supplied
  curated dataset paths without importing `airss_arena`.
- The package does not include the large curated minsep/volume dataset by
  default; users can generate or provide it explicitly.
- Offline dataset tools can build raw minsep/volume rows from MP docs, curate one
  row per formula by lowest `energy_above_hull`, and emit a summary JSON.
- Offline baseline training can consume the dataset and write a small
  `baseline_bundle.json` plus aggregated target and training summary files.
- Exact dataset estimates return volume per atom, total volume, required canonical
  minsep pairs, and source metadata for a requested formula.
- Baseline estimates run without torch and return the same output shape as exact
  dataset estimates.
- Reference-derived estimates compute median volume per atom and pair minima from
  ASE `Atoms` inputs.
- Auto `NFORM` respects atom budget and max-nform constraints, matching arena's
  intended behavior from `../airss-arena/src/airss_arena/generation_backend.py:338`.
- `ap run search --diagnose --volume-minsep-source ...` prints buildcell text containing
  formula, varvol, minsep, and nform directives while preserving unrelated seed
  directives.
- Existing manual `--target-volume` behavior remains unchanged when
  `--volume-minsep-source none` or no estimator options are provided.
- No torch import is required for dataset, baseline, reference, CLI diagnose, or
  seed-transform tests.
- `ruff check` passes for touched Python files.

## Potential Risks and Mitigations

1. **Buildcell syntax regressions for nested minsep or random nform.**
   Mitigation: route all generated seed objects through existing `SeedAtoms` and
   `BuildcellParam` serialization where possible, and add explicit tests for
   the final directive text.

2. **Changing legacy formula sampling behavior.**
   Mitigation: make estimate resolution opt-in, preserve existing tests, and add
   regression tests showing manual `--target-volume` output is unchanged.

3. **Unsupported formulas in the bundled baseline model.**
   Mitigation: raise clear `ValueError`/`ClickException` messages naming missing
   elements and supported elements, following arena's validation behavior in
   `../airss-arena/src/airss_arena/formula_model/inference.py:29`.

4. **Package-data omissions in wheels or source distributions.**
   Mitigation: add resource existence tests for the small runtime bundle and
   update both wheel and sdist packaging rules in `pyproject.toml`. Keep the
   large curated dataset outside the wheel by default.

5. **Dependency creep from accidentally porting torch/deep code.**
   Mitigation: keep deep predictor code out of the runtime module, avoid imports
   from arena's `deep.py` and `inference.py` torch path, and add a test or import
   check that baseline estimate loading does not import torch.

6. **Ambiguity between `VARVOL` semantics and total target volume.**
   Mitigation: document the chosen mapping and retain diagnostics showing
   volume per atom, total formula-unit volume, and final scaled `VARVOL`.
   Compare generated values against arena's `build_seed_atoms()` behavior in
   `../airss-arena/src/airss_arena/generation_backend.py:403`.

7. **Reference-derived minsep misses periodic-image contacts.**
   Mitigation: port arena's current repeat-based behavior first for parity, then
   consider a later improvement using ASE/pymatgen neighbor-list MIC logic if
   tests reveal gaps.

8. **Training script accidentally imports torch through arena's combined train module.**
   Mitigation: implement a baseline-only `airsspy` trainer that imports only the
   ported data, feature, and ridge modules. Add an import smoke test confirming
   torch is not imported for baseline training.

## Alternative Approaches

1. **Dataset-only first.**
   This is the smallest slice and avoids model inference entirely, but it would
   only work for formulas present in the curated dataset and would not satisfy
   the "simple model" direction as well as baseline inference.

2. **Baseline model through an optional extra.**
   This keeps core smaller, but the baseline path uses only existing core
   dependencies plus JSON assets, so an extra may add unnecessary friction.

3. **Keep estimate generation in `airsspy.search` only.**
   This reduces new modules but would make `search.py` too broad. A separate
   estimator module gives a cleaner Python API and keeps CLI, formula filtering,
   and estimate inference testable independently.

4. **Import `airss_arena` directly.**
   This is explicitly rejected for this implementation. It would pull in arena's
   project assumptions, Python 3.11 requirement, and torch dependency path,
   conflicting with `airsspy` packaging and the requested direct migration.

5. **Bundle the full curated dataset.**
   This would make exact lookup available out of the box, but it would add about
   59 MB of package data and couple release artifacts to MP-derived data. The
   preferred plan is to bundle only a small trained baseline bundle and provide
   reproducible dataset generation/curation commands.

## Out Of Scope For This Plan

- Deep/torch volume/minsep predictor inference.
- `ap run search --code ml` integrated buildcell-plus-ML relaxation.
- TorchSim autobatching, CuEq, constant-volume, and hydrostatic-strain controls.
- Arena project directory orchestration, reference pulling, matching, and
  analysis pipelines.

# Symmetrix MACE Backend

## Objective

Support Symmetrix as an explicit ML backend for ASE-based geometry optimisation
and single-point calculations in airsspy, while preserving the existing
torch-sim default and ASE fallback behavior. The intended user-facing syntax is
`symmetrix:mace:<model>` so it is visibly parallel to the existing plain
`mace:<model>` torch-sim path and the explicit `ase:mace:<model>` fallback
documented in `docs/reference/cli.md:91` and routed in
`src/airsspy/cli/cmd_run.py:740`. Only MACE models should be accepted by this
backend, and MACE model identifiers should use the same names already accepted
by the torch-sim MACE path, such as `medium`, `medium-mpa-0`, and local model
paths handled by `src/airsspy/jf/ml_runners.py:692`.

Assumptions: the request's "purchasing backend" refers to the existing
torch-sim backend, "torture scene" refers to torch-sim model naming, and
"metrics" refers to Symmetrix. The prior Symmetrix exploration found that the
ASE calculator is `Symmetrix(model_file, dtype="float64", use_kokkos=True,
**kwargs)` in `/home/bonan/appdir/symmetrix/symmetrix/source/symmetrix/symmetrix_calc.py:39`,
that it exposes energy, forces, and stress in
`/home/bonan/appdir/symmetrix/symmetrix/source/symmetrix/symmetrix_calc.py:36`,
and that it can convert torch MACE checkpoints to Symmetrix JSON through
`species`, `head`, and `num_spline_points` kwargs in
`/home/bonan/appdir/symmetrix/symmetrix/source/symmetrix/symmetrix_calc.py:54`.

## Implementation Plan

- [ ] 1. Define the backend syntax and parser contract.
  Treat `symmetrix:mace:<model>` as the only supported Symmetrix form, where
  the `mace` segment is mandatory and `<model>` follows the same
  foundation-model or file-path vocabulary as the existing torch-sim MACE
  loader in `src/airsspy/jf/ml_runners.py:684`. The parser should reject
  `symmetrix:sevennet:*`, `symmetrix:mattersim:*`, bare
  `symmetrix:<model>`, and unknown sub-backends with clear CLI errors so users
  do not assume Symmetrix is a general torch-sim replacement.
- [ ] 2. Add backend classification helpers in `cmd_run.py`.
  Place them near `_is_torchsim_model` at `src/airsspy/cli/cmd_run.py:740`.
  The helpers should distinguish torch-sim specs, explicit ASE fallback specs,
  and explicit Symmetrix specs before `_ensure_torchsim_available` runs. The
  intent is that `symmetrix:mace:medium` never requires `torch_sim`, while
  plain `mace:medium` keeps requiring torch-sim as it does at
  `src/airsspy/cli/cmd_run.py:750`.
- [ ] 3. Add a Symmetrix normalization path in `ml_runners.py`.
  This path should construct an ASE calculator spec for `symmetrix:Symmetrix`
  and calculator kwargs for MACE-only model resolution, dtype, Kokkos mode,
  species extraction, head selection, and spline resolution. Keep the current
  generic `_resolve_calculator` behavior at `src/airsspy/jf/ml_runners.py:59`
  for explicit `ase:` fallback only.
- [ ] 4. Share MACE model-name resolution with torch-sim.
  For named foundation models, resolve through the same MACE foundation-model
  selectors currently used when torch-sim loads `mace:<model>` in
  `src/airsspy/jf/ml_runners.py:706`, while preserving local checkpoint paths
  recognized at `src/airsspy/jf/ml_runners.py:697`. Avoid introducing
  Symmetrix-specific aliases unless they are only internal canonicalization to
  the same names users already pass to torch-sim.
- [ ] 5. Make Symmetrix JSON extraction cacheable and deterministic.
  Because Symmetrix JSON files are element-set-specific according to
  `/home/bonan/appdir/symmetrix/symmetrix/README.md:34`, derive the required
  species from each structure or batch before constructing the calculator.
  Cache converted JSON models under a stable airsspy cache directory keyed by
  model identity, species, head, dtype-relevant settings, and spline count, and
  reuse the JSON across structures with the same element set. Local checkpoint
  paths should be keyed by resolved path and file metadata to avoid stale cache
  reuse.
- [ ] 6. Reuse the existing ASE ML relax and single-point runners.
  `AirssMlRelaxRunner` already attaches an ASE calculator, runs FIRE or BFGS,
  applies optional pressure through `ExpCellFilter`, writes extxyz, and records
  relaxation metadata at `src/airsspy/jf/ml_runners.py:250`.
  `AirssMlSinglePointRunner` already writes energy, forces, stress, and
  pressure metadata at `src/airsspy/jf/ml_runners.py:182`. Use these existing
  runners after calculator normalization rather than adding a separate batch
  execution path.
- [ ] 7. Thread Symmetrix handling through all local run commands.
  The commands currently branch to a torch-sim batch path when
  `_is_torchsim_model` is true in `src/airsspy/cli/cmd_run.py:1931`,
  `src/airsspy/cli/cmd_run.py:2281`, and `src/airsspy/cli/cmd_run.py:1584`.
  Symmetrix specs should stay on the ASE runner path while still accepting the
  same `--optimizer`, `--fmax`, `--max-iterations`, and `--pressure` controls.
  `--device` and `--batch-size` should either be ignored with explicit debug
  logging or rejected for Symmetrix, with a preference for clear rejection if
  no GPU-backed Symmetrix device selection is available.
- [ ] 8. Preserve task-document and RES output behavior.
  `compose_ml_task_doc` already reads extxyz, calculates pressure and enthalpy,
  emits REM metadata, and writes force annotations into `.res` files at
  `src/airsspy/jf/ml_runners.py:359`. Extend REM metadata so outputs identify
  the normalized calculator as Symmetrix MACE while retaining the user's
  original model spec, allowing downstream ranking and conversion to behave
  unchanged.
- [ ] 9. Add focused unit tests for parser and normalization behavior.
  Extend `tests/test_cli.py:1207` style coverage to verify that plain
  `mace:medium` remains torch-sim, `ase:mace:medium` remains generic ASE
  fallback, `symmetrix:mace:medium` routes to the ASE runner path, non-MACE
  Symmetrix specs fail clearly, and Symmetrix specs do not call
  `_ensure_torchsim_available`.
- [ ] 10. Add runner tests with fake modules.
  Extend `tests/test_jf_ml_runners.py:35` and the torch-sim fake-module style
  around `tests/test_jf_ml_runners.py:386` to inject a fake
  `symmetrix.Symmetrix` calculator, verify model path and kwargs propagation,
  verify species-aware cache key construction if caching is implemented in
  airsspy, and verify single-point plus relax outputs keep energy, forces,
  stress, pressure, and relaxation metadata.
- [ ] 11. Add a guarded optional integration test.
  The test should be skipped unless `symmetrix`, `mace-torch`, and a cached or
  local MACE model are available, then perform a tiny periodic structure
  single-point and a short relaxation. This mirrors the repository's convention
  that optional executable or ML backends should be guarded, as described in
  `AGENTS.md:220`.
- [ ] 12. Update user-facing documentation and dependency notes.
  Add CLI examples beside the existing ML examples in `docs/reference/cli.md:91`,
  `docs/reference/cli.md:108`, and `docs/reference/cli.md:120`; explain that
  `symmetrix:mace:<model>` is MACE-only, uses the same MACE model names as
  `mace:<model>`, and uses ASE optimizers rather than torch-sim batching.
  Update installation docs near `docs/getting-started/installation.md:49` to
  note Symmetrix is an optional extra or external install if it cannot be
  declared on PyPI cleanly.

## Verification Criteria

- `ap run relax --cell "*.cell" --code ml --calculator symmetrix:mace:medium --optimizer FIRE` uses the ASE ML relax runner and does not try to import `torch_sim`.
- `ap run sp --cell "*.res" --code ml --calculator symmetrix:mace:medium-mpa-0` produces extxyz and `.res` outputs with energy, forces when available, stress when available, pressure metadata, and REM lines identifying Symmetrix MACE.
- `ap run crud --workdir . --code ml --calculator symmetrix:mace:medium --nostop` consumes hopper jobs through the non-batch ASE runner path and preserves existing success/failure artifact movement.
- `symmetrix:mace:<model>` accepts the same MACE model identifiers as the torch-sim MACE path, including named foundation models and local checkpoint files, while `symmetrix:sevennet:*` and `symmetrix:<model>` fail with actionable error messages.
- Existing torch-sim behavior remains unchanged: plain `mace:medium` still routes to torch-sim and still errors clearly when torch-sim is missing.
- Existing explicit ASE behavior remains unchanged: `ase:mace:medium` still normalizes to the MACE ASE calculator as currently tested in `tests/test_cli.py:1212`.
- Unit tests cover backend parsing, calculator normalization, MACE-only validation, cache key behavior, single-point metadata, relaxation metadata, and the three CLI command entry points.
- Documentation shows all three ML backend forms clearly: torch-sim default `mace:<model>`, generic ASE fallback `ase:mace:<model>`, and Symmetrix MACE `symmetrix:mace:<model>`.

## Potential Risks and Mitigations

1. **Symmetrix model conversion is species-specific and may be expensive.**
   Mitigation: cache converted JSON files by model identity and species set, and construct one calculator per species set rather than repeatedly converting per structure.

2. **Foundation-model selector behavior could drift from torch-sim MACE naming.**
   Mitigation: centralize MACE model-name resolution in one helper used by both torch-sim MACE loading and Symmetrix conversion, and add tests for the model names already present in CLI tests, especially `medium` and `medium-mpa-0`.

3. **Symmetrix may not support every MACE architecture accepted by MACE itself.**
   Mitigation: surface conversion failures as clear backend errors that include the model id and note that Symmetrix currently supports only compatible MACE models; keep the generic ASE fallback available for unsupported MACE variants.

4. **Device and batch options could mislead users.**
   Mitigation: document that Symmetrix uses the ASE runner path, validate or explicitly ignore `--device` and `--batch-size`, and keep torch-sim batching only for plain `mace:<model>` specs.

5. **Optional dependency availability may vary across Python versions and systems.**
   Mitigation: keep Symmetrix imports lazy, guard integration tests, and document the external install path if packaging cannot provide a reliable optional extra.

6. **Cache invalidation mistakes could use a stale JSON model after a checkpoint changes.**
   Mitigation: include resolved checkpoint path, file size, modification time, species, head, and spline count in cache keys for local files; include canonical model id and extraction parameters for named foundation models.

## Alternative Approaches

1. **Use `ase:symmetrix:Symmetrix@...` as the only interface.**
   This requires fewer parser changes but fails the MACE-only requirement, exposes too much generic constructor detail, and does not provide model-name consistency with the torch-sim MACE path.

2. **Make `mace:<model>` automatically choose Symmetrix when torch-sim is unavailable.**
   This is convenient but changes the meaning of the existing default backend and could surprise users who expect torch-sim batching or GPU behavior. It also weakens error messages around optional dependencies.

3. **Add a separate top-level code such as `--code symmetrix`.**
   This makes backend selection explicit but duplicates the ML runner surface and would require parallel `run relax`, `run sp`, and CRUD logic even though Symmetrix naturally fits the existing ASE calculator path.

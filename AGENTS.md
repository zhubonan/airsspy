# AGENTS.md

This file provides guidance to coding agents when working with code in this
repository.

Venv for development is at .venv, activate with `source .venv/bin/activate` before any command. Always use `uv pip install` instead of `pip install` for installing packages.

## Worktree virtualenv bootstrap

Worktree checkouts may start with an empty or partial `.venv`. When that
happens, mirror the dependency versions from the main repository instead of
letting uv resolve a fresh environment. This is especially important for torch:
the worktree must use the same torch build as the main repo.

Recommended workflow:

```bash
# From the worktree checkout
uv pip freeze --python /home/bonan/appdir/airsspy/.venv/bin/python > /tmp/airsspy-main-freeze.txt
sed 's#^-e file:///home/bonan/appdir/airsspy$#-e file://'"$(pwd)"'#' \
  /tmp/airsspy-main-freeze.txt > /tmp/airsspy-worktree-freeze.txt
uv venv --python /home/bonan/appdir/airsspy/.venv/bin/python --clear .venv
uv pip install --python .venv/bin/python --no-deps -r /tmp/airsspy-worktree-freeze.txt
```

If running inside a sandboxed worktree, run the final `uv pip install` outside
the sandbox when necessary so uv can use its normal cache. Verify torch before
running ML-related tests:

```bash
/home/bonan/appdir/airsspy/.venv/bin/python -c "import torch; print(torch.__version__)"
.venv/bin/python -c "import torch; print(torch.__version__)"
```

## Project Overview

airsspy is a Python library providing an ASE-based interface for Ab initio Random Structure Searching (AIRSS). It wraps the external `buildcell` executable (not bundled) and lets users construct search seeds programmatically, generate random structures, parse `.res` output files, run distributed searches via jobflow, and analyse results. Licensed under GPLv2.

## Build, Install, and Test Commands

```bash
# Install in editable mode with dev dependencies (includes test tools)
uv pip install -e ".[dev]"

# Install with docs dependencies
uv pip install -e ".[docs]"

# Run the full test suite
pytest tests/ -q

# Run a single test file
pytest tests/test_seed.py

# Run a specific test
pytest tests/test_seed.py::test_some_function

# Lint and format (ruff, configured for line-length=88, target py39)
ruff check src/
ruff format src/

# Type check
mypy src/ --ignore-missing-imports

# Run all pre-commit hooks
pre-commit run --all-files
```

## Architecture

### Core modules (`src/airsspy/`)

- **`seed.py`** — `SeedAtoms` (extends `ase.Atoms`) holds a `BuildcellParam` instance (`gentags`) for cell-level buildcell parameters and per-atom `SeedAtomTag` entries stored in the `atom_gentags` array. Tag descriptors (`BoolTag`, `GenericTag`, `RangeTag`, `NestedRangeTag`) map Python attribute assignments to buildcell keyword syntax. `SeedAtoms.write_seed()` serializes to a `.cell` file that `buildcell` reads.

- **`build.py`** — `Buildcell` wraps the external `buildcell` binary via `subprocess.Popen`. `Buildcell.generate()` pipes a seed's `.cell` content to `buildcell` and parses the output back into an ASE `Atoms` object.

- **`restools.py`** — `RESFile` class and helper functions for parsing CASTEP `.res` files (TITL/CELL/SFAC blocks), extracting energies/volumes/spacegroups, and converting to ASE or pymatgen structures. Supports `from_packed()` for concatenated RES files, `from_file()`, `from_string()`, and round-trip via `to_res_lines()`.

- **`common.py`** — Defines `BuildcellError` and `RelaxError` exceptions.
- **`log.py`** — Shared CLI logging setup.
- **`utils.py`** — General utilities (string parsing, file pattern matching, k-point calculation).
- **`search.py`** — Reusable local-search helpers: formula sampling/injection for buildcell seeds, charge-neutral formula filtering, target-volume handling, and RSS post-relax pruning based on energy-per-atom and fingerprints.

### Code-specific tools (`src/airsspy/`)

- **`casteptools.py`** — CASTEP output parsing and file management: `parse_dot_castep()` (geometry convergence), `parse_param()`, `RASH_prepare_seed()`, `extract_REM()`, `extract_result()`, `write_converge()`, `push_cell()`, `castep_finish_ok()`, `gulp_relax_finish_ok()`. Exception classes: `CastepRunError`, `CastepSkip`, `CastepManualTimedout`.

- **`abacustools.py`** — ABACUS utilities: parses ABACUS running logs and STRU files, converts CASTEP `.cell` content to ABACUS `STRU`, detects output logs, extracts REM metadata, and composes AIRSS-compatible result documents / `.res` output.

- **`scftools.py`** — `SCFInfo` class for parsing SCF convergence data from `.castep` files. Uses `ScfData` namedtuples per geometry step. Provides `get_summary()` statistics and `plot_scf()`/`plot_conv()` matplotlib plotting methods.

- **`gulptools.py`** — GULP utilities: `geom_opt_progress()` parses geometry optimisation, `check_gulp()` monitors for divergence/overflow, `guarded_gulp()` runs GULP with automatic termination.

- **`fullrelax.py`** — `FullRelax` class implementing self-consistent CASTEP relaxation with restart capability, alternating cell constraints, and state persistence. Includes `parse_geom_text_output()` for `.geom` files and `geom_to_cell()` for converting last configuration to cell blocks.

- **`tools/modcell.py`** — `replace_block()` and `modify_cell()` for programmatically modifying CASTEP cell files.

### Scheduler (`src/airsspy/`)

- **`scheduler.py`** — Unified scheduler interface: `Scheduler` base class, `Slurm`, `SGE`, and `Dummy` implementations. `Scheduler.get_scheduler()` factory auto-detects the current environment.

### Analysis sub-package (`src/airsspy/analysis/`)

- **`collect.py`** — DataFrame collection utilities: `collect_res_in_df()`, `combine_res_cryan()` (cryan CLI wrapper), `read_ca()` (ca command parser), `read_stream()`, `get_minsep_range()`, `get_entry()` (ComputedEntry factory), `export_dataframe_as_res()`, `get_pressure_gpa()`.

- **`hull.py`** — Phase diagram plotting: `PlotlyPDPlotter` extends pymatgen's `PDPlotter` with Plotly backend for ternary/binary convex hulls.

- **`query.py`** — Bridge from jobflow store to analysis: `collect_results_df()` collects results from a `SearchStore` into an analysis-ready DataFrame.

### Jobflow integration (`src/airsspy/jf/`)

- **`documents.py`** — Pydantic output models: `AirssJobDoc` (one per job, contains N results) and `AirssResultDoc` (per-structure data). `RelaxOutcome` enum for relaxation status. Both makers produce the same `AirssJobDoc` type.

- **`runners.py`** — Pure computation classes: `AirssCastepRelaxRunner`, `AirssCastepSinglePointRunner`, `AirssGulpRelaxRunner`, `AirssPp3RelaxRunner`, `AirssAbacusRelaxRunner`, `AirssAbacusSinglePointRunner`, `AirssScriptRelaxRunner`, `run_buildcell()`, `clean_files()`, and `compose_task_doc()`. Usable standalone or within Makers.

- **`ml_runners.py`** — Optional machine-learning potential runners. `AirssMlRelaxRunner` and `AirssMlSinglePointRunner` support ASE-compatible calculators and optional `torch_sim` batch backends; `compose_ml_task_doc()` writes `.res` output from extxyz results with energy/forces/stress metadata.

- **`jobs.py`** — Jobflow Makers: `AirssSearchMaker` (build+relax N structures per job), `AirssRelaxMaker` (relax N provided structures), `AirssValidateMaker`. All produce `AirssJobDoc` output.

- **`store.py`** — `SearchStore` query layer wrapping maggma `MongoStore`. Methods: `retrieve_project()`, `retrieve_project_df()`, `list_projects()`, `list_seeds()`, `show_struct_counts()`, `throughput_summary()`.

### CLI (`src/airsspy/cli/`)

Entry point: `ap` command (registered in `pyproject.toml`).

- **`main.py`** — Top-level `airss` Click group registered as the `ap` entry point. Global options include `--db-host`, `--db-port`, `--db-name`, `-v/--verbose`, and `-q/--quiet`.
- **`cmd_deploy.py`** — `deploy search` and `deploy relax` commands (both with `--dryrun`).
- **`cmd_db.py`** — `db list-projects`, `db list-seeds`, `db summary`, `db throughput`, `db retrieve-project`.
- **`cmd_check.py`** — `check airss`, `check scheduler`, `check database`.
- **`cmd_tools.py`** — `tools modcell`.
- **`cmd_run.py`** — Local non-jobflow AIRSS runner, similar in role to `airss.pl`. Provides `run search`, `run relax`, and `run sp`; search supports CASTEP/GULP/PP3/ABACUS, while relax supports CASTEP/GULP/PP3/ABACUS/ML and single-point supports CASTEP/ABACUS/ML. Also supports build-only mode, packing results, MPI launcher options, walltime-buffer stopping, formula sampling (`--formula`, `--elements`, `--target-volume`, oxidation states), and RSS pruning (`--prune*` options).
- **`cmd_rank.py`** — `rank` command for ranking structures by enthalpy per formula unit. Reads `.res` (packed) or extxyz from stdin/file args. Core options: `-t/--top`, `-de/--delta-e`, `-f/--formula` (exact formula, comma-separated elements like `Si,O`, or glob like `Si*`), `-nr/--absolute`, `-l/--long-labels`, `-s/--summary`, `-u/--unite` (fingerprint merging), `--input-format`, `--fingerprint-cutoff` (default 10.0 Å), `--unite-no-zweight` (Z-weighting is on by default). Pipeline controls: `--filter-name` (glob on label), `-p/--pressure` (external pressure in GPa, applies PV correction), `--unite-ethresh` (energy threshold in eV/atom for pre-merge filtering, default 0.1), `--unite-output` (directory for merged groups), `--unite-output-format` (`res` or `extxyz`).
- **`cmd_pack.py`** — `pack` and `unpack` commands for text-based file concatenation and splitting. `pack` concatenates files (glob args or `--from-dir`); `unpack` splits a packed file into a directory. RES unpack splits on TITL→END blocks; extxyz unpack splits per structure. `--format` flag controls output format (`res` or `extxyz`).
- **`cmd_convert.py`** — `convert` command for RES ↔ extxyz conversion. Bulk: `ap convert packed.res output.xyz` or `ap convert input.xyz output_dir/`. Single: `ap convert -l Si-002 packed.res Si-002.res`.

### Ranking module (`src/airsspy/ranking.py`)

Standalone ranking logic with no pandas dependency. Key components:

- **`StructureRecord`** dataclass — lightweight record with TITL metadata + species counts. Properties: `reduced_formula` (Hill system with cryan ordering: C first, H second if C present, O always last, rest by atomic number), `n_formula_units` (GCD of species counts), `enthalpy_per_fu`, `volume_per_fu`. Private `_merged_peers` field tracks structures merged into this record by `eliminate_similar()` (used for `--unite-output`).
- **`_parse_res_fast()`** — parses only TITL + species counts, avoids pymatgen overhead.
- **`rank_structures()`** — groups by formula, computes relative energies, filters by delta_e.
- **`summary_structures()`** — most stable per composition (cryan `-s` equivalent).
- **`apply_external_pressure(records, pressure_gpa)`** — adds PV correction to enthalpy (`H = E + PV`), enabling ranking at non-zero external pressure.
- **`filter_by_name(records, pattern)`** — glob-based filter on structure labels (TITL label field).
- **`filter_by_formula(records, formula)`** — formula filter supporting three modes: exact reduced formula match, comma-separated element set (e.g. `Si,O` matches any composition containing only Si and O), or glob pattern on reduced formula string.
- **`prefilter_records(records, ethresh)`** — removes structures whose energy per atom is more than `ethresh` eV/atom above the best in their composition group. Applied before merge to avoid expensive fingerprint computation on clearly unstable structures.
- **`_compute_distance_fingerprint()`** — computes sorted distance fingerprint using pymatgen's `get_all_neighbors(cutoff)` (all periodic images within cutoff). Returns `np.ndarray | None`. Optional Z-weighting scales each distance as `d * zmax^2 / (Z_i * Z_j)` to distinguish atom-type pairs. Z-weighting is now enabled by default.
- **`eliminate_similar()`** — fingerprint-based merging (cryan `-u` equivalent). Accepts a `cutoff` parameter (default 4.0 A) controlling the neighbour search radius. Uses numpy-vectorised fingerprint comparison. Records which structures were merged via `_merged_peers` on surviving records.
- **`format_header()`/`format_rank_line()`** — cryan-compatible tabular output.

#### Ranking pipeline order

The `ap rank` command processes structures through a fixed pipeline. The order matters because each stage reduces the work for subsequent stages:

1. **Read** — parse input (packed `.res` or extxyz)
2. **Apply pressure** (`-p`) — add PV correction if external pressure specified
3. **Filter** — `--filter-name` (glob on label) and `-f/--formula` (composition filter)
4. **Pre-rank** (`--unite-ethresh`) — `prefilter_records()` removes high-energy outliers per composition
5. **Merge** (`-u`) — `eliminate_similar()` merges fingerprint-duplicate structures
6. **Post-rank** (`-de`, `-t`) — `rank_structures()` computes relative energies and applies final delta_e/top-N filtering
7. **Display** — tabular output, summary, or file export

### Conversion module (`src/airsspy/convert.py`)

Lossless RES ↔ extxyz conversion with force support:

- **`res_to_extxyz()`** — packed .res → single extxyz. Forces from atom columns 8-10 → `SinglePointCalculator`.
- **`extxyz_to_res()`** — extxyz → unpacked individual .res files (one per structure, named `<label>.res`).
- **`extract_structure()`** — extract single structure by label from packed .res or extxyz. Fast TITL-only scan first, then loads only the match. Output format from extension.
- **`_atoms_to_res_lines()`** — builds RES lines from ASE Atoms, appends force columns `fx fy fz` if present.
- **`_parse_res_forces()`** — reads force columns 8-10 from .res atom lines.

Force column placement in .res: `Symbol index x y z occ [spin] [fx fy fz]`. Safe because cryan reads only columns 1-5, cabal reads up to column 7 (spin).

### Search helpers (`src/airsspy/search.py`)

Local AIRSS search support outside jobflow:

- **Formula sampling** — `FormulaSamplingOptions`, `FormulaSamplingContext`, `build_formula_sampling_context()`, and `make_seed_text_transform()` inject sampled `#FORMULA`/`#VARVOL` directives into seed text. Formula pools can come from explicit formulas or enumerated element coefficient combinations and can be filtered by seed `#NATOM`/`#NFORM` constraints and oxidation-state charge neutrality.
- **RSS pruning** — `RssPruneOptions`, `RssCandidate`, `candidate_from_res()`, `pool_statistics()`, `should_flush_prune_pool()`, `select_pruned_candidates()`, and `prune_relaxed_pool()` retain low-energy, fingerprint-unique candidates from relaxed pools before final output.

### ABACUS and ML support

- ABACUS integration lives in both `abacustools.py` and `jf/runners.py`. ABACUS runners expect an ABACUS input file suffix from `cmd_deploy.SUFFIX_MAP`, generate/consume `STRU` and ABACUS output directories, then compose AIRSS-style `.res` files.
- ML support is optional under the `[ml]` extra (`torch-sim`). The default ML path uses dynamic ASE calculator specs such as `module.path:ClassName@model`; when `torch_sim` is installed, `ml_runners.py` has batch relax/static helpers for model specs like `mace:medium`.

### Reference implementations

The AIRSS reference Fortran code lives at `~/appdir/airss-git/`:

- **`internal/cryan/src/cryan.f90`** — Structure ranking/analysis tool. Key routines: `read_res()` (line 293), `rank()` (line 1667), `compositions()` (line 852), `summary()`, `eliminate()` (line 2200). Reads atom lines as `symbol, index, x, y, z` only (ignores columns 6+). Formula ordering uses `elements_alpha` with O→huge, C→0, H→0.1 if C present.
- **`internal/cabal/src/cabal.f90`** — Structure manipulation tool. `read_res()` (line 557) and `write_res()` (line 664) show full RES read/write. Atom format: `Symbol index x y z occ [spin]`, with write format `a4,i4,3f17.13,f4.1[,f7.2]`. TITL format: `label P V H spin spin_abs [dos] nat (symm) n - copies`.

### Key design patterns

- `SeedAtoms` extends `ase.Atoms`; its `build` property returns the `BuildcellParam` object for setting cell-level tags like `varvol`, `symmops`, `numat`. Per-atom tags (e.g., `posamp`, `tagname`) are accessed via `seed.atom_tags[i]`.
- The tag descriptor system uses Python descriptors (`__set__`/`__get__`) that validate and serialize values into the buildcell cell-file format.
- `Buildcell` requires the AIRSS `buildcell` executable on `$PATH`. Tests that call `buildcell` will fail if it is not installed.
- **Multi-structure jobs**: Both `AirssSearchMaker` and `AirssRelaxMaker` produce a single `AirssJobDoc` containing N `AirssResultDoc` entries, reducing JobStore document count for high-throughput searches.
- **Runners are jobflow-free**: `AirssCastepRelaxRunner` and `run_buildcell()` can be used standalone without jobflow orchestration.
- Local CLI search uses the shared runner layer rather than jobflow. It writes `.cell`, relaxation outputs, and `.res` files in the selected work directory and honors a `stop` file plus scheduler walltime checks.

### Dependency tiers

```toml
# Core (always installed)
ase>=3.17, castepinput>=0.1, click>=8.0, jobflow>=0.1, maggma>=0.50, numpy>=1.20,
pandas>=1.3.0, plotly, pymatgen>=2022.0.0, spglib>=1.16, tabulate, tqdm

# Optional
[dev]  pytest, pytest-cov, ruff, mypy, pre-commit, twine
[ml]   torch-sim
[docs] sphinx, pydata-sphinx-theme, myst-nb, sphinx-autodoc2, sphinx-design, etc.
```

Python >=3.9.

## Testing

Tests live in `tests/` with `conftest.py` providing fixtures (`al_atoms`, `tmpfile`). pytest is configured in `pyproject.toml` with `--strict-markers --strict-config`; the `e2e` marker is reserved for tests that run real CASTEP/ABACUS executables. Tests that depend on external programs such as `buildcell`, CASTEP, ABACUS, GULP, PP3, or optional ML backends should be skipped/guarded when the executable or package is unavailable. The suite covers core seed/build/RES utilities, CASTEP/GULP/ABACUS tools, ranking/conversion, search helpers, scheduler, jobflow documents/jobs/runners/store, ML runner helpers, and CLI commands (using `click.testing.CliRunner`).

## Mandatory final review subagent

After implementing changes and before final response:

1. Spawn a read-only reviewer subagent.
2. Ask it to review the current branch diff against the base branch.
3. The reviewer must focus on:
   - correctness
   - regressions
   - missing tests
   - race conditions
   - security issues
   - unintended unrelated edits
4. Wait for the reviewer result.
5. Fix any high-confidence issue.
6. Summarize the review result in the final answer.

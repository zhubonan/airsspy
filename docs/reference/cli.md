# Command Line Interface

airsspy installs the `ap` console script. The Click command group is named
`airss` internally, but the packaged entry point is `ap`.

```bash
ap --help
ap --version
```

Global options control logging and the MongoDB store used by the jobflow
commands:

```bash
ap -v --db-host localhost --db-port 27017 --db-name airss db summary
```

Use `-v` or `-vv` for more logging, and `-q`, `-qq`, or `-qqq` for less.

## Command Overview

| Command | Purpose |
| --- | --- |
| `ap check` | Check AIRSS executables, scheduler detection, or database access. |
| `ap run` | Run local non-jobflow searches, relaxations, CRUD workers, and single-points. |
| `ap rank` | Rank `.res` or extxyz structures by enthalpy, merge duplicates, or run convex-hull analysis. |
| `ap convert` | Convert packed `.res` files to extxyz, extxyz to `.res`, or extract one labelled structure. |
| `ap pack` / `ap unpack` | Concatenate or split `.res`, `.xyz`, and `.extxyz` structure files. |
| `ap deploy` | Build jobflow search or relaxation flows and store results in MongoDB. |
| `ap db` | Query and retrieve jobflow results from MongoDB. |
| `ap tools modcell` | Replace the structure in a CASTEP `.cell` template with an ASE-readable structure. |
| `ap tools volume-minsep-*` | Build/curate volume-minsep datasets and train the lightweight baseline predictor. |

## Environment Checks

```bash
ap check airss
ap check scheduler --allow-dummy
ap --db-host localhost --db-name airss check database
```

`check airss` looks for `buildcell`, `cabal`, and `castep_relax` on `PATH`.
`check scheduler` auto-detects Slurm, SGE, or the dummy fallback. `check
database` opens the configured jobflow store and lists visible projects.

## Local Runs

`ap run` is the local AIRSS interface. It writes files in `--workdir`, honors a
`stop` file for searches, checks scheduler walltime with `--walltime-buffer`,
and can pack successful `.res` outputs into `packed.res`.

### Random Searches

```bash
ap run search --seed Si --code castep --nmax 100 --pack
ap run search --seed Si --code gulp --exe ggulp --cluster
ap run search --seed Si --code vasp --potcar-dir /path/to/potpaw --potcar-map Si=Si
ap run search --seed seed --build-only --nmax 20
```

Required files:

| Code | Required input files |
| --- | --- |
| `castep` | `<seed>.cell`, `<seed>.param` |
| `gulp` | `<seed>.cell`, `<seed>.lib` |
| `pp3` | `<seed>.cell`, `<seed>.pp` |
| `abacus` | `<seed>.cell`, `<seed>.INPUT` |
| `vasp` | `<seed>.cell`, `<seed>.INCAR`; optional `<seed>.KPOINTS` |

Useful options include `--pressure`, `--max-iterations`, `--build-timeout`,
`--mpinp`, `--keep`, and `--diagnose`.

Formula sampling can inject sampled `#FORMULA` and `#VARVOL` directives into
the seed sent to buildcell:

```bash
ap run search --seed seed --formula SiO2,Si2O3 --target-volume Si=20 --target-volume O=12
ap run search --seed seed --elements Li,P,O --max-coeff 4 --oxidation-state Li=1,P=5,O=-2
ap run search --seed seed --elements Li,P,O --diagnose 3
```

Volume/minsep estimates can additionally inject `#MINSEP` and automatic
`#NFORM` directives. The large curated dataset is user-supplied or generated
offline; airsspy does not bundle it by default.

```bash
ap run search --seed seed --formula SiO2 \
  --volume-minsep-source dataset \
  --volume-minsep-dataset generation/minsep_vol_dataset_curated.json \
  --diagnose 1

ap run search --seed seed --formula SiO2 \
  --volume-minsep-source baseline \
  --volume-minsep-bundle formula_model_artifacts/baseline/baseline_bundle.json
```

Use `--volume-scale`, `--minsep-scale-low`, `--minsep-scale-high`,
`--max-atoms`, and `--max-nform` to tune generated directives.

Post-relax RSS pruning is available with `--prune`. The main controls are
`--prune-pool-size`, `--prune-keep-fraction`, `--prune-dedup-tol`,
`--prune-fingerprint-cutoff`, and `--prune-zweight`.

### Volume/Minsep Dataset Tools

The volume/minsep tools are offline utilities for reproducing a curated dataset
and training the small non-torch baseline bundle:

```bash
ap tools volume-minsep-build-dataset --mp-docs mp-stable-docs.json --output minsep_vol_dataset.json
ap tools volume-minsep-curate-dataset --dataset minsep_vol_dataset.json --mp-docs mp-stable-docs.json --output minsep_vol_dataset_curated.json
ap tools volume-minsep-train-baseline --dataset minsep_vol_dataset_curated.json --output-dir formula_model_artifacts
```

The curation step selects one row per reduced formula by lowest
`energy_above_hull`.

### Relax Existing Structures

```bash
ap run relax --cell "*.cell" --seed Si --code castep --pack
ap run relax --cell "*.res" --seed Si --code vasp --potcar-dir /path/to/potpaw
ap run relax --cell "*.res" --code ml --calculator mace:medium --device cuda --batch-size 16
ap run relax --cell "*.cell" --code ml --calculator ase:mace:medium --optimizer BFGS
```

`run relax` accepts single-structure `.cell` and `.res` inputs. Packed `.res`
files matched by the glob are skipped. For `.res` inputs with non-ML codes,
airsspy looks for a root template such as `Si.cell` so it can preserve the
non-structural cell settings while replacing the lattice and positions.

Supported relaxation codes are `castep`, `gulp`, `pp3`, `abacus`, `vasp`, and
`ml`. The `--singlepoint` flag reuses this command for single-point calculations
where supported: `castep`, `abacus`, `vasp`, and `ml`.

### Single-Point Calculations

```bash
ap run sp --cell "*.cell" --seed Si --code castep
ap run sp --cell "*.res" --code ml --calculator mace:medium --batch-size 32
ap run sp --cell "*.res" --seed Si --code vasp --potcar-dir /path/to/potpaw
```

`run sp` supports `castep`, `abacus`, `vasp`, and `ml`. RES input is currently
supported for VASP and ML single-points.

### CRUD Queue Worker

```bash
ap run crud --workdir . --code castep
ap run crud --workdir . --code vasp --singlepoint --potcar-dir /path/to/potpaw
ap run crud --workdir . --code ml --calculator mace:medium --batch-size 8 --nostop
```

The CRUD worker consumes queued `hopper/*-*.res` files, converts each claimed
RES structure into backend inputs, and moves outputs into `good_castep/` or
`bad_castep/`. Use `STOP_CRUD` to stop a long-running worker. `--cycle` can
requeue CASTEP-like electronic-minimisation failures.

## Ranking And Hull Analysis

`ap rank` reads packed `.res` from files or standard input. It also reads extxyz
files when the extension or `--input-format extxyz` says so.

```bash
cat *.res | ap rank
ap rank packed.res -t 20 -de 0.05
ap rank structures.xyz --input-format extxyz --label-field structure_id
ap rank packed.res -f SiO2 --filter-name "SiO2-*"
```

Important filters and display modes:

| Option | Behaviour |
| --- | --- |
| `-t, --top` | Keep only the top N structures per ranking output. |
| `-de, --delta-e` | Filter by relative energy per atom in eV. |
| `-f, --formula` | Match an exact formula, a comma-separated element set, or a glob. |
| `-fu, --formula-unit` | Filter by exact number of formula units. |
| `-sn, --speciesnumber` | Filter by exact number of species. |
| `-in, --ionsnumber` | Filter by ion count; negative values mean up to `abs(N)`. |
| `-nr, --absolute` | Show absolute enthalpy instead of relative enthalpy. |
| `-l, --long-labels` | Do not shorten labels in output. |
| `-s, --summary` | Show the most stable structure per composition. |

Duplicate merging uses distance fingerprints:

```bash
ap rank packed.res -u 0.1 --unite-ethresh 0.05
ap rank packed.res -u 0.1 --unite-output united --unite-output-format res
ap rank packed.res -u 0.1 --distance 8.0 --unite-zweight
```

Pathological low-energy pruning is opt-in:

```bash
ap rank packed.res --prune-pathological --pathology-sigma-factor 3
```

Convex-hull, or Maxwell construction, mode uses pymatgen:

```bash
ap rank packed.res -m -el Si,O
ap rank packed.res -m -el Si,O --plot hull.html
ap rank packed.res -m -el Si,O -s
```

`-f/--formula` cannot be combined with `-m/--maxwell`; use `-el/--element-list`
to choose the phase diagram elements.

## Convert, Pack, And Unpack

```bash
ap convert packed.res structures.xyz
ap convert structures.xyz res_dir/
ap convert packed.res res_dir/
ap convert -l Si-002 packed.res Si-002.res
ap convert -l Si-002 packed.res Si-002.xyz
```

The conversion path preserves RES lines when possible. Forces in RES atom
columns 8-10 are transferred to extxyz calculators, and extxyz forces are
written back to RES force columns.

```bash
ap pack *.res packed.res
ap pack --from-dir res_files packed.res
ap pack --from-dir xyz_files packed.xyz
ap unpack packed.res split_res/
ap unpack packed.xyz split_xyz/
```

`pack` concatenates files. `unpack` splits RES on `TITL` ... `END` blocks and
splits extxyz by atom-count structure blocks. Use `ap convert` when changing
between RES and extxyz formats.

## Jobflow Deployment

`ap deploy` creates jobflow flows and runs them locally with a MongoDB-backed
`JobStore`.

```bash
ap --db-name airss deploy search --seed Si --project test --num 100 --code castep
ap deploy search --seed Si --project test --num 20 --n-structures 5 --dryrun
ap deploy search --seed Si --project test --num 10 --code vasp --dryrun
```

`deploy search` uses the same seed companion suffixes as local search, except
ML is not a deploy backend. Supported deploy codes are `castep`, `gulp`, `pp3`,
`abacus`, and `vasp`.

```bash
ap deploy relax --project test --cell "*.cell" --param Si.param --code castep
ap deploy relax --project test --cell "*.xyz" --param Si.INCAR --code vasp --dryrun
```

`deploy relax` can read `.cell` files directly or ASE-readable structures. With
`--base-cell`, it overlays each input structure onto a CASTEP cell template.

## Database Queries

The `db` commands read from the MongoDB database selected by the global
`--db-host`, `--db-port`, and `--db-name` options:

```bash
ap db list-projects
ap db list-seeds --project test
ap db summary --project test
ap db throughput --past-days 7
ap db retrieve-project --project test --output res_out
```

`retrieve-project` writes stored RES contents as `<struct_name>.res`.

## Utility Tools

`tools modcell` prints a CASTEP `.cell` file made from a template plus a new
structure readable by ASE:

```bash
ap tools modcell template.cell relaxed.xyz > new.cell
```

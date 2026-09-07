# Auto Volume And Minsep

`ap run search` can derive buildcell `#VARVOL`, `#MINSEP`, and `#NFORM`
directives from a volume/minsep estimate while preserving unrelated seed
settings.

The first implementation is intentionally lightweight. It supports exact lookup
from a user-supplied curated JSON dataset, a baseline JSON bundle trained from
that dataset, or reference structures. It does not require torch, and the large
curated dataset is not bundled with airsspy.

## Preview Generated Seed Text

Use `--diagnose` before a long search:

```bash
ap run search --seed seed --build-only --formula SiO2 \
  --volume-minsep-source dataset \
  --volume-minsep-dataset generation/minsep_vol_dataset_curated.json \
  --diagnose 3
```

The diagnostic output shows the selected formula, estimate provenance, volume
per atom, total volume, generated minsep ranges, resolved `#NFORM`, and the
final buildcell input.

## Use A Baseline Bundle

After training a baseline bundle, point search at the generated
`baseline_bundle.json`:

```bash
ap run search --seed seed --formula SiO2 \
  --volume-minsep-source baseline \
  --volume-minsep-bundle artifacts/baseline/baseline_bundle.json
```

`--volume-scale`, `--minsep-scale-low`, `--minsep-scale-high`, `--max-atoms`,
and `--max-nform` tune the generated buildcell directives.

## Generate Training Data

If you have a Materials Project document dump with structures and
`energy_above_hull`, build and curate the dataset offline:

```bash
ap tools volume-minsep-build-dataset \
  --mp-docs generation/mp-stable-docs.json \
  --output generation/minsep_vol_dataset.json

ap tools volume-minsep-curate-dataset \
  --dataset generation/minsep_vol_dataset.json \
  --mp-docs generation/mp-stable-docs.json \
  --output generation/minsep_vol_dataset_curated.json
```

The curation step keeps the lowest `energy_above_hull` row for each reduced
formula.

## Train The Baseline Predictor

Train the non-torch ridge baseline from the raw or curated dataset:

```bash
ap tools volume-minsep-train-baseline \
  --dataset generation/minsep_vol_dataset_curated.json \
  --output-dir formula_model_artifacts
```

The command writes `baseline/baseline_bundle.json`, aggregated target JSON
files, and a training summary under the output directory.

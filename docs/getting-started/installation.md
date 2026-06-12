# Installation

This guide covers how to install airsspy and its dependencies.

## Requirements

- Python 3.9 or later
- pip or uv package manager
- AIRSS buildcell executable (for structure generation)

## Installing airsspy

### Using pip (recommended for users)

```bash
pip install airsspy
```

### Using uv (faster, modern alternative)

```bash
uv pip install airsspy
```

### Development Installation

For development or to get the latest version from source:

```bash
git clone https://github.com/zhubonan/airsspy.git
cd airsspy
pip install -e ".[dev,docs]"
```

Or with uv:

```bash
git clone https://github.com/zhubonan/airsspy.git
cd airsspy
uv pip install -e ".[dev,docs]"
```

## Optional Dependencies

airsspy has several optional dependency groups:

- **`dev`**: Development and test tools (pytest, pytest-cov, ruff, mypy,
  pre-commit, twine)
- **`ml`**: Optional machine-learning potential support through torch-sim-atomistic
- **`docs`**: Documentation building tools (Sphinx, themes, extensions)

Install with optional dependencies:

```bash
pip install "airsspy[ml]"
pip install "airsspy[docs]"
```

## Installing AIRSS

To generate random structures, you need the AIRSS buildcell executable.

### Download and Build AIRSS

1. Download AIRSS from the [official website](http://www.mtg.msm.cam.ac.uk/Codes/AIRSS)

2. Extract and build:

```bash
tar -xzf airss-*.tgz
cd airss-*
make
make install
```

3. Add AIRSS binaries to your PATH:

```bash
export PATH=<airss root folder>/bin:$PATH
```

Add this line to your `~/.bashrc` or `~/.zshrc` to make it permanent.

### Verify Installation

Check that buildcell is accessible:

```bash
which buildcell
```

## Dependencies

airsspy depends on:

- **ase** - Atomic Simulation Environment
- **castepinput** (≥0.1) - CASTEP input file handling
- **click** - Command-line interface
- **jobflow** and **maggma** - Workflow and MongoDB-backed result storage
- **numpy** (≥1.20) - Numerical computing
- **pandas** - DataFrame analysis helpers
- **plotly** - Plotting support
- **pymatgen** (≥2022.0.0) - Interface for reading SHELX files
- **spglib** - Symmetry detection
- **tabulate** and **tqdm** - CLI tables and progress display

These are automatically installed when you install airsspy.

## Verifying Installation

Test your installation:

```python
import airsspy
print(airsspy.__version__)

# Test creating a seed
from airsspy import SeedAtoms
seed = SeedAtoms('C4', cell=[5, 5, 5], pbc=True)
print("airsspy is working!")
```

## Troubleshooting

### buildcell not found

If you get an error about buildcell not being found:

1. Check that AIRSS is installed: `which buildcell`
2. Ensure AIRSS bin directory is in your PATH
3. Try specifying the full path to buildcell if needed

### Import errors

If you encounter import errors:

1. Verify you're using Python 3.9 or later: `python --version`
2. Try installing in a fresh virtual environment
3. Ensure optional dependencies are installed when needed, for example
   `pip install "airsspy[ml]"` for ML workflows

## Next Steps

Now that airsspy is installed, continue to the [Quickstart Tutorial](quickstart.md) to perform your first structure search.

# API Reference

Complete API documentation for airsspy modules and classes.

## Core Modules

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} airsspy.seed
:link: ../apidocs/airsspy/airsspy.seed
:link-type: doc

SeedAtoms, BuildcellParam, and tag classes for seed generation.
:::

:::{grid-item-card} airsspy.build
:link: ../apidocs/airsspy/airsspy.build
:link-type: doc

Buildcell class for interfacing with AIRSS buildcell executable.
:::

:::{grid-item-card} airsspy.restools
:link: ../apidocs/airsspy/airsspy.restools
:link-type: doc

RESFile and utilities for handling CASTEP .res files.
:::

:::{grid-item-card} airsspy.utils
:link: ../apidocs/airsspy/airsspy.utils
:link-type: doc

Utility functions for file parsing and data processing.
:::

:::{grid-item-card} Command Line
:link: cli
:link-type: doc

Installed `ap` command groups and workflow examples.
:::
::::

## Quick Navigation

### Main Classes

- {py:class}`~airsspy.SeedAtoms` - Core class for creating structure seeds
- {py:class}`~airsspy.Buildcell` - Wrapper for buildcell executable
- {py:class}`~airsspy.RESFile` - RES file handler with structure and metadata

### Tag Classes

- {py:class}`~airsspy.seed.BuildcellParam` - Global buildcell parameters
- {py:class}`~airsspy.seed.SeedAtomTag` - Per-atom generation tags
- {py:class}`~airsspy.seed.SeedAtom` - Individual atom with tags

### RES File Functions

- {py:func}`~airsspy.save_airss_res` - Save structure in RES format
- {py:func}`~airsspy.read_res_atoms` - Read RES file to Atoms
- {py:func}`~airsspy.read_res_pmg` - Read RES file to pymatgen Structure
- {py:func}`~airsspy.extract_res` - Extract metadata from RES file
- {py:func}`~airsspy.parse_titl` - Parse TITL line

### Utility Functions

- {py:func}`~airsspy.get_spacegroup_atoms` - Detect space group
- {py:func}`~airsspy.get_minsep` - Calculate minimum separations
- {py:func}`~airsspy.format_minsep` - Format minsep dictionary
- {py:func}`~airsspy.calc_kpt_tuple_recip` - Calculate k-point mesh

### Command Line Modules

- {py:mod}`airsspy.cli.main` - Top-level `ap` command and global options
- {py:mod}`airsspy.cli.cmd_run` - Local search, relaxation, CRUD, and single-point commands
- {py:mod}`airsspy.cli.cmd_rank` - Structure ranking and convex-hull command
- {py:mod}`airsspy.cli.cmd_convert` - RES/extxyz conversion command
- {py:mod}`airsspy.cli.cmd_pack` - Pack and unpack commands
- {py:mod}`airsspy.cli.cmd_deploy` - jobflow deployment commands
- {py:mod}`airsspy.cli.cmd_db` - jobflow database query commands

## Module Index

The generated module index is available at [API Modules](../apidocs/index).

## Usage Examples

### Creating Seeds

```python
from airsspy import SeedAtoms

# Create a seed
seed = SeedAtoms('Si8', cell=[5, 5, 5], pbc=True)
seed[0].num = 8
seed.gentags.minsep = 2.0
seed.gentags.varvol = 20

# Generate structure
atoms = seed.build_random_atoms()
```

### Using Buildcell

```python
from airsspy import Buildcell

# Create Buildcell instance
builder = Buildcell(seed)

# Generate with custom timeout
atoms = builder.generate(timeout=30)
```

### Working with RES Files

```python
from airsspy import RESFile, save_airss_res

# Read RES file
res = RESFile.from_file('structure.res')
print(f"Formula: {res.formula}")
print(f"Enthalpy: {res.enthalpy:.4f} eV")

# Save structure
info_dict = {
    'uid': 'test-1',
    'P': 0.0,
    'V': atoms.get_volume(),
    'H': atoms.get_potential_energy(),
    'nat': len(atoms),
    'sym': 'P1'
}
save_airss_res(atoms, info_dict, 'output.res')
```

## Type Hints

airsspy uses type hints throughout for better IDE support and code clarity. Import types from modules:

```python
from airsspy.seed import SeedAtoms, BuildcellParam
from airsspy.build import Buildcell
from airsspy.restools import RESFile, TitlInfo
```

## See Also

- [Getting Started](../getting-started/index.md) - Learn the basics
- [How-To Guides](../how-to/index.md) - Task-oriented guides
- [Explanations](../explanation/index.md) - Understand concepts

```{toctree}
:hidden:

cli
../apidocs/index
```

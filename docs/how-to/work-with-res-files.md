---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# How to Work with RES Files

Learn how to read, write, and analyze AIRSS .res output files.

## What are RES Files?

RES files are the standard output format from AIRSS searches. They contain:
- Structure information (atoms, cell, symmetry)
- Computed properties (energy, volume, pressure)
- Metadata (unique ID, space group)

## Setup: Create Sample RES Files

First, let's create some sample RES files for demonstration:

```{code-cell} ipython3
:tags: [hide-output]

from airsspy import SeedAtoms, save_airss_res
from ase.calculators.lj import LennardJones
from ase.filters import UnitCellFilter
from ase.optimize import BFGS
import tempfile
import os

# Create temporary directory for examples
temp_dir = tempfile.gettempdir()
example_dir = os.path.join(temp_dir, 'res_examples')
os.makedirs(example_dir, exist_ok=True)

# Generate and save 3 example structures
seed = SeedAtoms('Al', cell=[2, 2, 2], pbc=True)
seed[0].num = 4
seed.gentags.minsep = 1.5

calc = LennardJones()

for i in range(3):
    atoms = seed.build_random_atoms()
    if atoms is not None:
        atoms.set_calculator(calc)
        opt = BFGS(UnitCellFilter(atoms), logfile=None)
        opt.run(fmax=0.05)

        info_dict = {
            'uid': f'Al-example-{i+1}',
            'P': 0.0,
            'V': atoms.get_volume(),
            'H': atoms.get_potential_energy(),
            'nat': len(atoms),
            'sym': 'P1'
        }

        res_path = os.path.join(example_dir, f'Al-example-{i+1}.res')
        save_airss_res(atoms, info_dict, res_path, force_write=True)

print(f"Created 3 example RES files in {example_dir}")
```

## Reading RES Files

### Using RESFile Class

The most convenient way to work with RES files:

```{code-cell} ipython3
from airsspy import RESFile

# Load a RES file
res_path = os.path.join(example_dir, 'Al-example-1.res')
res = RESFile.from_file(res_path)

# Access properties
print(f"Formula: {res.formula}")
print(f"Enthalpy: {res.enthalpy:.4f} eV")
print(f"Volume: {res.volume:.2f} Å³")
print(f"Pressure: {res.pressure:.2f} GPa")
print(f"Space group: {res.symm}")
print(f"Number of atoms: {res.natoms}")
```

### Loading Multiple Files

Process a directory of RES files:

```{code-cell} ipython3
from pathlib import Path
from airsspy import RESFile

# Load all example RES files
res_files = Path(example_dir).glob('*.res')
structures = []

for fpath in res_files:
    res = RESFile.from_file(str(fpath))
    structures.append(res)

# Sort by enthalpy
structures.sort(key=lambda r: r.enthalpy)

print("Lowest energy structures:")
for i, res in enumerate(structures, 1):
    print(f"  {i}. {res.label}: {res.enthalpy:.4f} eV, V={res.volume:.2f} Å³")
```

### Fast Loading (Metadata Only)

Load only the TITL line without parsing structure:

```{code-cell} ipython3
res = RESFile.from_file(res_path, only_titl=True)

# Access metadata (faster)
print(f"Label: {res.label}")
print(f"Enthalpy: {res.enthalpy:.4f} eV")
print(f"Volume: {res.volume:.2f} Å³")
# Note: structure will be None with only_titl=True
print(f"Structure loaded: {res.structure is not None}")
```

### Using extract_res Function

Extract metadata as a dictionary:

```{code-cell} ipython3
from airsspy import extract_res

info = extract_res(res_path)

print("Extracted metadata:")
print(f"  UID: {info['uid']}")
print(f"  Enthalpy: {info['H']:.4f} eV")
print(f"  Volume: {info['V']:.2f} Å³")
print(f"  Pressure: {info['P']:.2f} GPa")
print(f"  N atoms: {info['nat']}")
print(f"  Space group: {info['sym']}")
```

## Working with Structures

### Get ASE Atoms Object

Convert RES structure to ASE Atoms:

```{code-cell} ipython3
from airsspy import RESFile

res = RESFile.from_file(res_path)

# Get as ASE Atoms
atoms = res.atoms

# Now use with ASE
print(f"Chemical formula: {atoms.get_chemical_formula()}")
print(f"Number of atoms: {len(atoms)}")
print(f"Cell volume: {atoms.get_volume():.2f} Å³")
print(f"Positions shape: {atoms.get_positions().shape}")
```

### Get Pymatgen Structure

Get the pymatgen Structure object:

```{code-cell} ipython3
res = RESFile.from_file(res_path)

# Access pymatgen Structure
structure = res.structure

# Use pymatgen functionality
print(f"Composition: {structure.composition}")
print(f"Reduced formula: {structure.composition.reduced_formula}")
print(f"Density: {structure.density:.2f} g/cm³")
```

## Writing RES Files

### Using save_airss_res

Save an ASE Atoms object as a RES file:

```{code-cell} ipython3
from airsspy import save_airss_res, SeedAtoms

# Create a simple structure
seed = SeedAtoms('C', cell=[3, 3, 3], pbc=True)
seed[0].num = 4
seed.gentags.minsep = 1.5

atoms = seed.build_random_atoms()

if atoms is not None:
    # Prepare metadata
    info_dict = {
        'uid': 'carbon-demo',
        'P': 0.0,                           # Pressure (GPa)
        'V': atoms.get_volume(),            # Volume (Å³)
        'H': -20.5,                         # Enthalpy/energy (eV)
        'nat': len(atoms),                  # Number of atoms
        'sym': 'P1'                         # Space group
    }

    # Save to file
    demo_path = os.path.join(example_dir, 'carbon-demo.res')
    save_airss_res(atoms, info_dict, demo_path, force_write=True)
    print(f"✓ Saved to {demo_path}")

    # Read it back to verify
    res = RESFile.from_file(demo_path)
    print(f"✓ Verified: {res.label}, {res.natoms} atoms, H={res.enthalpy:.4f} eV")
```

### Auto-naming

Let the function generate the filename from uid:

```{code-cell} ipython3
:tags: [hide-output]

if atoms is not None:
    info_dict = {
        'uid': 'auto-named-structure',
        'P': 0.0,
        'V': atoms.get_volume(),
        'H': -15.0,
        'nat': len(atoms),
        'sym': 'P1'
    }

    # Change to example directory
    import os
    orig_dir = os.getcwd()
    os.chdir(example_dir)

    # Filename will be 'auto-named-structure.res'
    save_airss_res(atoms, info_dict, force_write=True)

    os.chdir(orig_dir)
    print(f"Saved as auto-named-structure.res")
```

## Batch Processing

### Analyze Multiple Results

Process and rank all structures in a search:

```{code-cell} ipython3
from pathlib import Path
from airsspy import RESFile

# Load all RES files from example directory
results = []
for fpath in Path(example_dir).glob('*.res'):
    try:
        res = RESFile.from_file(str(fpath))
        results.append({
            'file': fpath.name,
            'label': res.label,
            'enthalpy': res.enthalpy,
            'volume': res.volume,
            'formula': res.formula,
            'symm': res.symm
        })
    except Exception as e:
        print(f"Failed to load {fpath.name}: {e}")

# Sort by enthalpy
results.sort(key=lambda x: x['enthalpy'])

# Print summary
print(f"\n{'Rank':<6} {'Label':<25} {'Enthalpy':<12} {'Volume':<10} {'Symm':<10}")
print("-" * 75)
for i, r in enumerate(results[:5], 1):
    print(f"{i:<6} {r['label']:<25} {r['enthalpy']:<12.4f} {r['volume']:<10.2f} {r['symm']:<10}")
```

## Complete Example: Post-Processing

Here's a complete example for analyzing AIRSS results:

```{code-cell} ipython3
from pathlib import Path
from airsspy import RESFile

def analyze_search_results(directory, top_n=5):
    """Analyze and summarize AIRSS search results"""

    # Load all structures
    structures = []
    search_path = Path(directory)

    for fpath in search_path.glob('*.res'):
        try:
            res = RESFile.from_file(str(fpath))
            structures.append(res)
        except Exception as e:
            print(f"Warning: Could not load {fpath.name}: {e}")

    if not structures:
        print("No structures found!")
        return

    # Sort by enthalpy
    structures.sort(key=lambda r: r.enthalpy)

    print(f"\nAnalyzed {len(structures)} structures")
    print(f"Energy range: {structures[-1].enthalpy - structures[0].enthalpy:.4f} eV")
    print(f"\nTop {min(top_n, len(structures))} structures:\n")

    # Print summary table
    print(f"{'#':<4} {'Label':<25} {'Energy':<12} {'V/atom':<10} {'Symm':<12}")
    print("-" * 75)

    for i, res in enumerate(structures[:top_n], 1):
        v_per_atom = res.volume / res.natoms
        print(f"{i:<4} {res.label:<25} {res.enthalpy:<12.4f} {v_per_atom:<10.2f} {res.symm:<12}")

    # Symmetry distribution
    symmetries = [res.symm for res in structures]
    print(f"\nSpace group distribution:")
    unique_symms = sorted(set(symmetries))
    for sg in unique_symms[:5]:  # Show top 5
        count = symmetries.count(sg)
        print(f"  {sg}: {count} structure(s)")

# Run analysis on our example directory
analyze_search_results(example_dir, top_n=5)
```

## Troubleshooting

### File Format Errors

If RESFile fails to load:

1. Check the file is a valid RES file
2. Ensure TITL line is present
3. Try `only_titl=True` for corrupted structures

### Missing spglib

Some symmetry functions require spglib:

```bash
pip install spglib
```

## See Also

- [Quickstart Tutorial](../getting-started/quickstart.md) - Complete workflow example
- [Generating Structures](generate-structures.md) - Creating structures to save
- [API: RESFile](../reference/index.md) - Full RESFile documentation

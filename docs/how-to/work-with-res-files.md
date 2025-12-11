# How to Work with RES Files

Learn how to read, write, and analyze AIRSS .res output files.

## What are RES Files?

RES files are the standard output format from AIRSS searches. They contain:
- Structure information (atoms, cell, symmetry)
- Computed properties (energy, volume, pressure)
- Metadata (unique ID, space group)

## Reading RES Files

### Using RESFile Class

The most convenient way to work with RES files:

```python
from airsspy import RESFile

# Load a RES file
res = RESFile.from_file('structure.res')

# Access properties
print(f"Formula: {res.formula}")
print(f"Enthalpy: {res.enthalpy:.4f} eV")
print(f"Volume: {res.volume:.2f} Å³")
print(f"Pressure: {res.pressure:.2f} GPa")
print(f"Space group: {res.symm}")
```

### Loading Multiple Files

Process a directory of RES files:

```python
from pathlib import Path
from airsspy import RESFile

res_files = Path('.').glob('*.res')
structures = []

for fpath in res_files:
    res = RESFile.from_file(str(fpath))
    structures.append(res)

# Sort by enthalpy
structures.sort(key=lambda r: r.enthalpy)

print("Lowest energy structures:")
for res in structures[:5]:
    print(f"{res.label}: {res.enthalpy:.4f} eV")
```

### Fast Loading (Metadata Only)

Load only the TITL line without parsing structure:

```python
res = RESFile.from_file('structure.res', only_titl=True)

# Access metadata (faster)
print(f"Label: {res.label}")
print(f"Enthalpy: {res.enthalpy}")
# Note: structure will be None with only_titl=True
```

### Using extract_res Function

Extract metadata as a dictionary:

```python
from airsspy import extract_res

info = extract_res('structure.res')

print(info['uid'])      # Unique identifier
print(info['H'])        # Enthalpy
print(info['V'])        # Volume
print(info['P'])        # Pressure
print(info['nat'])      # Number of atoms
print(info['sym'])      # Space group
print(info['rem'])      # REM lines (list)
```

## Working with Structures

### Get ASE Atoms Object

Convert RES structure to ASE Atoms:

```python
from airsspy import RESFile

res = RESFile.from_file('structure.res')

# Get as ASE Atoms
atoms = res.atoms

# Now use with ASE
print(f"Chemical formula: {atoms.get_chemical_formula()}")
print(f"Number of atoms: {len(atoms)}")
print(f"Cell volume: {atoms.get_volume():.2f} Å³")

# Can use ASE I/O
from ase.io import write
write('structure.cif', atoms)
```

### Get Pymatgen Structure

Get the pymatgen Structure object:

```python
res = RESFile.from_file('structure.res')

# Access pymatgen Structure
structure = res.structure

# Use pymatgen functionality
print(f"Composition: {structure.composition}")
print(f"Density: {structure.density:.2f} g/cm³")
print(f"Lattice: {structure.lattice}")
```

## Writing RES Files

### Using save_airss_res

Save an ASE Atoms object as a RES file:

```python
from airsspy import save_airss_res
from ase import Atoms

# Your atoms object (from calculation, generation, etc.)
atoms = Atoms('C4', positions=[[0,0,0], [1,1,1], [2,2,2], [3,3,3]],
              cell=[5, 5, 5], pbc=True)

# Prepare metadata
info_dict = {
    'uid': 'carbon-test-1',
    'P': 0.0,                           # Pressure (GPa)
    'V': atoms.get_volume(),            # Volume (Å³)
    'H': -20.5,                         # Enthalpy/energy (eV)
    'nat': len(atoms),                  # Number of atoms
    'sym': 'P1'                         # Space group
}

# Save to file
save_airss_res(atoms, info_dict, 'carbon-test-1.res')
```

### Overwrite Protection

By default, `save_airss_res` won't overwrite existing files:

```python
# This will raise FileExistsError if file exists
save_airss_res(atoms, info_dict, 'existing.res')

# Force overwrite
save_airss_res(atoms, info_dict, 'existing.res', force_write=True)
```

### Auto-naming

Let the function generate the filename from uid:

```python
info_dict = {'uid': 'my-structure', 'P': 0.0, 'V': 64.0,
             'H': -10.0, 'nat': 8, 'sym': 'Fm-3m'}

# Filename will be 'my-structure.res'
save_airss_res(atoms, info_dict)
```

## Analyzing MINSEP

Extract minimum separation information:

```python
from airsspy import RESFile, format_minsep

res = RESFile.from_file('structure.res')

# Get minsep as dictionary
minsep_dict = res.get_minsep(string=False)
print(minsep_dict)
# Example: {'C-C': 1.52, 'C-O': 1.43, 'O-O': 2.45}

# Get minsep as formatted string
minsep_str = res.get_minsep(string=True)
print(minsep_str)
# Example: "C-C=1.52 C-O=1.43 O-O=2.45"
```

## Advanced Reading

### Using read_res_atoms

Low-level function that returns atoms and TITL info:

```python
from airsspy import read_res_atoms

with open('structure.res') as f:
    lines = f.readlines()

titl_info, atoms = read_res_atoms(lines)

print(f"Label: {titl_info.label}")
print(f"Enthalpy: {titl_info.enthalpy}")
print(f"Atoms: {len(atoms)}")
```

### Using read_res_pmg

Get pymatgen Structure with additional data:

```python
from airsspy import read_res_pmg

with open('structure.res') as f:
    lines = f.readlines()

titl_info, rem_lines, structure, spins = read_res_pmg(lines)

print(f"REM lines: {rem_lines}")
print(f"Structure: {structure.composition}")
print(f"Spins: {spins}")
```

## Batch Processing

### Analyze Multiple Results

Process and rank all structures in a search:

```python
from pathlib import Path
from airsspy import RESFile

# Load all RES files
results = []
for fpath in Path('.').glob('*.res'):
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
        print(f"Failed to load {fpath}: {e}")

# Sort by enthalpy
results.sort(key=lambda x: x['enthalpy'])

# Print summary
print(f"\n{'Rank':<6} {'Label':<20} {'Enthalpy':<12} {'Volume':<10} {'Symm':<10}")
print("-" * 70)
for i, r in enumerate(results[:10], 1):
    print(f"{i:<6} {r['label']:<20} {r['enthalpy']:<12.4f} {r['volume']:<10.2f} {r['symm']:<10}")
```

### Export to DataFrame

Convert results to pandas DataFrame for analysis:

```python
import pandas as pd
from pathlib import Path
from airsspy import RESFile

# Load structures
data = []
for fpath in Path('.').glob('*.res'):
    res = RESFile.from_file(str(fpath))
    data.append({
        'label': res.label,
        'enthalpy': res.enthalpy,
        'volume': res.volume,
        'pressure': res.pressure,
        'natoms': res.natoms,
        'symm': res.symm,
        'formula': res.formula
    })

# Create DataFrame
df = pd.DataFrame(data)
df = df.sort_values('enthalpy')

print(df.head(10))

# Save to CSV
df.to_csv('results_summary.csv', index=False)
```

## Symmetry Analysis

### Using get_spacegroup_atoms

Get symmetry information for an atoms object:

```python
from airsspy import get_spacegroup_atoms

atoms = ...  # Your atoms object

# Detect space group
sg_info = get_spacegroup_atoms(atoms, symprec=0.1, angle_tolerance=5.0)
print(f"Space group: {sg_info}")
```

## Complete Example: Post-Processing

Here's a complete example for analyzing AIRSS results:

```python
from pathlib import Path
from airsspy import RESFile, format_minsep

def analyze_search_results(pattern='*.res', top_n=10):
    """Analyze and summarize AIRSS search results"""

    # Load all structures
    structures = []
    for fpath in Path('.').glob(pattern):
        try:
            res = RESFile.from_file(str(fpath))
            structures.append(res)
        except Exception as e:
            print(f"Warning: Could not load {fpath}: {e}")

    if not structures:
        print("No structures found!")
        return

    # Sort by enthalpy
    structures.sort(key=lambda r: r.enthalpy)

    print(f"\nAnalyzed {len(structures)} structures")
    print(f"Energy range: {structures[-1].enthalpy - structures[0].enthalpy:.4f} eV")
    print(f"\nTop {top_n} structures:\n")

    # Print summary table
    print(f"{'#':<4} {'Label':<25} {'Energy':<12} {'V/atom':<10} {'Symm':<12}")
    print("-" * 75)

    for i, res in enumerate(structures[:top_n], 1):
        v_per_atom = res.volume / res.natoms
        print(f"{i:<4} {res.label:<25} {res.enthalpy:<12.4f} {v_per_atom:<10.2f} {res.symm:<12}")

    # Symmetry distribution
    symmetries = [res.symm for res in structures]
    print(f"\nSpace group distribution:")
    for sg in sorted(set(symmetries)):
        count = symmetries.count(sg)
        print(f"  {sg}: {count} structures")

# Run analysis
analyze_search_results(top_n=10)
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

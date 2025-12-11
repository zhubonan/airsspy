# Quickstart Tutorial

Learn the basics of airsspy by performing a simple random structure search.

## Prerequisites

- airsspy installed ([Installation Guide](installation.md))
- AIRSS buildcell executable in PATH
- Basic familiarity with Python and ASE

## Your First Structure Search

This tutorial demonstrates a complete AIRSS workflow: creating a seed, generating random structures, and optimizing them.

### Step 1: Import Required Modules

```python
from airsspy import SeedAtoms, Buildcell
from ase.calculators.lj import LennardJones
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter
```

### Step 2: Create a Seed Structure

A seed defines the search template. Let's search for an 8-atom aluminum structure:

```python
# Define the search seed
seed = SeedAtoms('Al', cell=[2, 2, 2], pbc=True)
seed.gentags.minsep = 1.5  # Minimum separation constraint

# Define per-atom tags - request 8 Al atoms
al = seed[0]
al.num = 8
```

**What's happening here:**
- `SeedAtoms('Al', ...)` creates a template with one Al atom
- `cell=[2, 2, 2]` sets the initial cell size (will be scaled by buildcell)
- `gentags.minsep = 1.5` sets minimum atomic separation
- `al.num = 8` tells buildcell to generate 8 aluminum atoms

### Step 3: Inspect the Seed File

You can view the seed file content:

```python
print('\n'.join(seed.get_cell_inp_lines()))
```

Output:
```
%BLOCK lattice_cart
2.0000000000  0.0000000000  0.0000000000
0.0000000000  2.0000000000  0.0000000000
0.0000000000  0.0000000000  2.0000000000
%ENDBLOCK lattice_cart
%BLOCK positions_abs
Al  0.0000000000 0.0000000000 0.0000000000 # Al0 % NUM=8
%ENDBLOCK positions_abs
#MINSEP=1.5
```

### Step 4: Generate a Random Structure

```python
random_atoms = seed.build_random_atoms()

# Check the generated structure
if random_atoms is not None:
    print(f"Generated {len(random_atoms)} atoms")
    print(f"Volume: {random_atoms.get_volume():.2f} Å³")
else:
    print("Failed to generate structure - check buildcell installation")
```

The buildcell program:
- Takes your 2×2×2 cell (volume = 8 Å³)
- Multiplies by 8 atoms per cell: target volume = 64 Å³
- Randomly varies the volume (typically ±50%)
- Places 8 atoms respecting the minsep constraint

### Step 5: Perform a Structure Search

Now let's generate and relax multiple random structures:

```python
def relax_structure(seed, calculator, n_structures=10):
    """
    Generate and relax multiple random structures
    """
    results = []

    for i in range(n_structures):
        # Generate random structure
        atoms = seed.build_random_atoms()
        if atoms is None:
            continue

        # Set up calculator
        atoms.set_calculator(calculator)

        # Optimize structure (both atomic positions and cell)
        optimizer = BFGS(UnitCellFilter(atoms), logfile=None)
        optimizer.run(fmax=0.05)

        results.append(atoms)

    return results

# Run the search with Lennard-Jones potential
lj_calc = LennardJones()
structures = relax_structure(seed, lj_calc, n_structures=20)
print(f"Successfully relaxed {len(structures)} structures")
```

### Step 6: Analyze Results

Check the energies and identify the lowest-energy structure:

```python
# Get energies
energies = [atoms.get_potential_energy() for atoms in structures]

# Find lowest energy structure
min_idx = energies.index(min(energies))
best_structure = structures[min_idx]

print(f"Lowest energy: {energies[min_idx]:.4f} eV")
print(f"Energy range: {max(energies) - min(energies):.4f} eV")
```

### Step 7: Check Symmetry (Optional)

If you have spglib installed, you can analyze the symmetry:

```python
from spglib import get_spacegroup

symmetries = [get_spacegroup(atoms, symprec=0.5) for atoms in structures]
print("Space groups found:")
for sym in set(symmetries):
    count = symmetries.count(sym)
    print(f"  {sym}: {count} structures")
```

### Step 8: Save Results

Save the best structure in AIRSS .res format:

```python
from airsspy import save_airss_res

# Prepare metadata
info_dict = {
    'uid': 'Al-search-best',
    'P': 0.0,  # Pressure
    'V': best_structure.get_volume(),
    'H': best_structure.get_potential_energy(),
    'nat': len(best_structure),
    'sym': 'P1'  # Would be determined by symmetry analysis
}

save_airss_res(best_structure, info_dict, 'Al-best.res')
print("Saved best structure to Al-best.res")
```

## Complete Example Code

Here's the full working example:

```python
from airsspy import SeedAtoms
from ase.calculators.lj import LennardJones
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter

# Create seed
seed = SeedAtoms('Al', cell=[2, 2, 2], pbc=True)
seed.gentags.minsep = 1.5
seed[0].num = 8

# Set up calculator
calc = LennardJones()

# Run search
results = []
for i in range(20):
    atoms = seed.build_random_atoms()
    if atoms is None:
        continue
    atoms.set_calculator(calc)
    opt = BFGS(UnitCellFilter(atoms), logfile=None)
    opt.run(fmax=0.05)
    results.append(atoms)

# Find best
energies = [atoms.get_potential_energy() for atoms in results]
best_idx = energies.index(min(energies))
print(f"Found {len(results)} structures")
print(f"Best energy: {energies[best_idx]:.4f} eV")
```

## Understanding buildcell Behavior

When you call `build_random_atoms()`:

1. **Cell volume**: The initial cell volume (8 Å³) is multiplied by the number of atoms per formula unit
2. **Volume variation**: buildcell randomly varies the volume (typically by ±50%)
3. **Random positions**: Atoms are placed randomly, respecting minsep constraints
4. **Random cell shape**: The cell shape is randomized (within constraints)
5. **Symmetry**: Random symmetry operations may be applied (controlled by `symmops`)

## Common Parameters

Here are some commonly used seed generation tags:

```python
seed.gentags.minsep = 2.0        # Minimum separation (Å)
seed.gentags.varvol = 20         # Volume variation (%)
seed.gentags.symmops = (2, 4)    # Number of symmetry operations
seed.gentags.compact = True      # Apply Niggli reduction
seed.gentags.nform = 1           # Number of formula units
```

## Next Steps

Now that you understand the basics, explore:

- [Creating Complex Seeds](../how-to/create-seeds.md) - Advanced seed creation techniques
- [Buildcell Parameters](../explanation/buildcell-parameters.md) - Complete parameter reference
- [Working with RES Files](../how-to/work-with-res-files.md) - Reading and analyzing AIRSS output
- [API Reference](../reference/index.md) - Detailed API documentation

## Troubleshooting

### buildcell not found

If you get an error about buildcell:
- Verify AIRSS is installed: `which buildcell`
- Check PATH includes AIRSS bin directory
- See [Installation Guide](installation.md)

### build_random_atoms() returns None

This happens when buildcell times out or fails:
- Increase the timeout: `seed.build_random_atoms(timeout=30)`
- Check your minsep constraints aren't too tight
- Verify your seed file is valid

### Structures too similar

If all structures converge to the same result:
- Increase volume variation: `seed.gentags.varvol = 50`
- Generate more structures
- Try different starting cell sizes
- Add symmetry operations: `seed.gentags.symmops = (1, 4)`

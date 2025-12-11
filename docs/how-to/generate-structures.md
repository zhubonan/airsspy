---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# How to Generate Random Structures

Learn different approaches for generating random structures with airsspy.

## Basic Generation

### Using SeedAtoms Directly

The simplest way to generate a structure:

```{code-cell} ipython3
from airsspy import SeedAtoms

seed = SeedAtoms('C4', cell=[3, 3, 3], pbc=True)
seed[0].num = 4
seed.gentags.minsep = 1.5

# Generate one random structure
atoms = seed.build_random_atoms()

if atoms is not None:
    print(f"Generated {len(atoms)} atoms")
    print(f"Volume: {atoms.get_volume():.2f} Å³")
else:
    print("Generation failed")
```

### Using Buildcell Class

For more control, use the Buildcell class directly:

```{code-cell} ipython3
from airsspy import SeedAtoms, Buildcell

seed = SeedAtoms('Si8', cell=[4, 4, 4], pbc=True)
seed[0].num = 8
seed.gentags.minsep = 2.0

# Create Buildcell instance
builder = Buildcell(seed)

# Generate structure with custom timeout
atoms = builder.generate(timeout=30)

if atoms is not None:
    print(f"Generated {len(atoms)} atoms")
    print(f"Formula: {atoms.get_chemical_formula()}")
else:
    print("Generation timed out")
```

## Handling Generation Failures

Buildcell may fail or timeout. Handle these cases:

### With fail_ok (Default)

```{code-cell} ipython3
# Returns None on failure (default behavior)
atoms = seed.build_random_atoms(fail_ok=True)

if atoms is None:
    print("Generation failed or timed out")
else:
    print(f"Success! Generated {len(atoms)} atoms")
```

### Raise Exception on Failure

```{code-cell} ipython3
from airsspy.common import BuildcellError

try:
    # This might fail intentionally with tight constraints
    tight_seed = SeedAtoms('C8', cell=[1, 1, 1], pbc=True)  # Very small cell
    tight_seed[0].num = 8
    tight_seed.gentags.minsep = 3.0  # Impossible to satisfy

    atoms = tight_seed.build_random_atoms(fail_ok=False, timeout=5)
except BuildcellError as e:
    print(f"Buildcell failed as expected: {type(e).__name__}")
    print("Constraints were too tight for the cell size")
```

## Batch Generation

### Generate Multiple Structures

Generate several structures in a loop:

```{code-cell} ipython3
seed = SeedAtoms('Al', cell=[2, 2, 2], pbc=True)
seed[0].num = 4
seed.gentags.minsep = 1.5

# Generate 3 structures (using fewer for faster docs build)
structures = []
for i in range(3):
    atoms = seed.build_random_atoms()
    if atoms is not None:
        structures.append(atoms)
        print(f"Structure {i+1}: {len(atoms)} atoms, V={atoms.get_volume():.2f} Å³")

print(f"\nTotal generated: {len(structures)}/{3}")
```

### With Success Rate Tracking

Track how often generation succeeds:

```{code-cell} ipython3
seed = SeedAtoms('C', cell=[3, 3, 3], pbc=True)
seed[0].num = 6
seed.gentags.minsep = 1.3

n_attempts = 5
successes = 0
failures = 0

for i in range(n_attempts):
    atoms = seed.build_random_atoms(timeout=10)
    if atoms is not None:
        successes += 1
    else:
        failures += 1

print(f"Success rate: {successes}/{n_attempts} ({100*successes/n_attempts:.1f}%)")
if failures > 0:
    print(f"Tip: If success rate is low, try:")
    print(f"  - Increasing cell size")
    print(f"  - Reducing minsep constraint")
    print(f"  - Increasing timeout")
```

## Advanced Options

### Custom Timeout

Control how long buildcell tries:

```{code-cell} ipython3
:tags: [hide-output]

# Short timeout for quick attempts
quick_atoms = seed.build_random_atoms(timeout=5)

# Longer timeout for difficult cases
patient_atoms = seed.build_random_atoms(timeout=60)
```

### Writing Intermediate Files

Keep the .cell file that buildcell uses:

```{code-cell} ipython3
:tags: [hide-output]

import tempfile
import os

# Write seed to temporary file
temp_dir = tempfile.gettempdir()
cell_file = os.path.join(temp_dir, 'temp.cell')

builder = Buildcell(seed)
atoms = builder.generate(timeout=10, write_cell=cell_file)

if atoms is not None and os.path.exists(cell_file):
    print(f"Seed written to: {cell_file}")
    print("You can inspect or modify this file")
```

## Troubleshooting

### Generation Always Fails

If `build_random_atoms()` consistently returns None:

1. **Check constraints are achievable:**
   ```python
   # Too tight - 8 atoms with 3 Å spacing in 3×3×3 cell is impossible
   # seed.gentags.minsep = 3.0

   # More reasonable
   seed.gentags.minsep = 1.5
   ```

2. **Increase cell size:**
   ```python
   # Instead of cell=[2, 2, 2]
   seed = SeedAtoms('C8', cell=[4, 4, 4], pbc=True)
   ```

3. **Increase timeout:**
   ```python
   atoms = seed.build_random_atoms(timeout=60)  # Wait longer
   ```

### Buildcell Not Found

If you get "buildcell not found":
- Verify AIRSS is installed: Run `which buildcell` in terminal
- Check PATH includes AIRSS bin directory
- See [Installation Guide](../getting-started/installation.md)

## Best Practices

1. **Test your seed**: Generate a few structures first to verify constraints work
2. **Reasonable timeout**: Start with 10-30 seconds, adjust based on complexity
3. **Handle None**: Always check if `build_random_atoms()` returns None
4. **Batch with care**: When generating many structures, consider parallelization
5. **Check buildcell output**: If failing repeatedly, inspect error messages

## See Also

- [Creating Seeds](create-seeds.md) - How to design effective seeds
- [Working with RES Files](work-with-res-files.md) - Saving and analyzing results
- [Buildcell Parameters](../explanation/buildcell-parameters.md) - Complete parameter reference

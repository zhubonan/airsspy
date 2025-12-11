# How to Generate Random Structures

Learn different approaches for generating random structures with airsspy.

## Basic Generation

### Using SeedAtoms Directly

The simplest way to generate a structure:

```python
from airsspy import SeedAtoms

seed = SeedAtoms('C4', cell=[3, 3, 3], pbc=True)
seed[0].num = 4
seed.gentags.minsep = 1.5

# Generate one random structure
atoms = seed.build_random_atoms()

if atoms is not None:
    print(f"Generated {len(atoms)} atoms")
    print(f"Volume: {atoms.get_volume():.2f} Å³")
```

### Using Buildcell Class

For more control, use the Buildcell class directly:

```python
from airsspy import SeedAtoms, Buildcell

seed = SeedAtoms('Si8', cell=[4, 4, 4], pbc=True)
seed[0].num = 8

# Create Buildcell instance
builder = Buildcell(seed)

# Generate structure with custom timeout
atoms = builder.generate(timeout=30)
```

## Handling Generation Failures

Buildcell may fail or timeout. Handle these cases:

### With fail_ok (Default)

```python
# Returns None on failure (default behavior)
atoms = seed.build_random_atoms(fail_ok=True)

if atoms is None:
    print("Generation failed or timed out")
else:
    print("Success!")
```

### Raise Exception on Failure

```python
from airsspy.common import BuildcellError

try:
    atoms = seed.build_random_atoms(fail_ok=False)
except BuildcellError:
    print("Buildcell failed to generate structure")
```

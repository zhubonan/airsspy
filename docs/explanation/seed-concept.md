# Understanding Seed Files

Learn what seed files are, how they work, and how to design effective seeds for structure searching.

## What is a Seed?

A seed file is a template that defines the search space for random structure generation. It specifies:

- Which elements to include
- How many atoms of each type
- The initial cell dimensions
- Constraints on the generated structures
- Symmetry requirements

Think of a seed as a recipe: it doesn't specify the exact structure, but rather the rules for generating possible structures.

## Seed File Format

AIRSS uses a CASTEP-style `.cell` file format. Here's a simple example:

```
%BLOCK lattice_cart
4.0 0.0 0.0
0.0 4.0 0.0
0.0 0.0 4.0
%ENDBLOCK lattice_cart

%BLOCK positions_abs
Si 0.0 0.0 0.0 # Si0 % NUM=4
%ENDBLOCK positions_abs

#MINSEP=2.0
#VARVOL=20
#SYMMOPS=2-4
```

The `#` is a comment symbol for CASTEP, but not for `buildcell`. So the actual `builcell` directives always has a `#` prefix. A `##` symbol will comment out the line for build CASTEP and `buildcell`.

### Components

1. **lattice_cart**: Initial cell vectors (will be randomized)
2. **positions_abs**: Atomic positions with generation tags
3. **Generation tags**: Comments starting with `#` (e.g., `#MINSEP`, `#VARVOL`)

## Creating Seeds with airsspy

In airsspy, you create seeds programmatically using the `SeedAtoms` class:

```python
from airsspy import SeedAtoms

# Create a seed
seed = SeedAtoms('Si', cell=[4, 4, 4], pbc=True)

# Set per-atom tags
seed[0].num = 4              # Generate 4 Si atoms

# Set global tags
seed.gentags.minsep = 2.0    # Minimum separation
seed.gentags.varvol = 20     # Volume variation
seed.gentags.symmops = (2, 4)  # 2-4 symmetry operations
```

This is equivalent to the file format shown above, but more convenient and less error-prone.

## How buildcell Uses Seeds

When you generate a structure from a seed:

1. **Reads composition**: Determines which atoms and how many
2. **Calculates target volume**: Based on cell size and number of atoms
3. **Randomizes cell**: Varies volume and cell shape
4. **Places atoms**: Randomly, respecting constraints
5. **Applies symmetry**: If symmetry operations requested
6. **Checks constraints**: Ensures minsep and other requirements met

If constraints can't be satisfied, buildcell will timeout and return failure.

## Global vs Per-Atom Tags

### Global Tags (gentags)

Global tags apply to the entire structure:

```python
seed.gentags.minsep = 2.0        # All atom pairs
seed.gentags.varvol = 20         # Volume variation
seed.gentags.symmops = (2, 4)    # Symmetry operations
seed.gentags.compact = True      # Niggli reduction
```

These control overall structure properties like symmetry, volume, and cell shape.

### Per-Atom Tags

Per-atom tags apply to specific atom types:

```python
# Access individual atoms
si = seed[0]
si.num = 4           # 4 Si atoms
si.posamp = 1.5      # Position randomization
si.fix = True        # Fix position

o = seed[1]
o.num = 8            # 8 O atoms
o.coord = (3, 4)     # 3-4 coordination
```

These control properties of specific atoms within the structure.

## Designing Effective Seeds

### Start Simple

Begin with minimal constraints:

```python
seed = SeedAtoms('C4', cell=[3, 3, 3], pbc=True)
seed[0].num = 4
seed.gentags.minsep = 1.2  # Basic constraint only
```

Add more constraints if needed based on results.

### Use Chemical Intuition

Apply realistic constraints based on chemistry:

```python
# Si-O system with realistic separations
seed = SeedAtoms('SiO2', cell=[5, 5, 5], pbc=True)
seed.gentags.minsep = (1.2, {
    'Si-O': 1.5,   # Si-O bonds typically ~1.6 Å
    'O-O': 2.0,    # O-O should be further apart
    'Si-Si': 2.5   # Si-Si even further
})
```

### Consider Symmetry

If you expect certain symmetry:

```python
# Likely to be cubic
seed.gentags.system = 'cubic'
seed.gentags.symmops = (4, 8)

# Or specific space group number
seed.gentags.symmno = 225  # Fm-3m
```

But don't over-constrain - you might miss the true structure!

### Volume Scaling

The initial cell volume is multiplied by the number of formula units:

```python
seed = SeedAtoms('C', cell=[2, 2, 2], pbc=True)
seed[0].num = 8

# Initial volume: 2³ = 8 Å³
# Target volume: 8 × 8 = 64 Å³
# Actual volume: 64 ± varvol%
```

Set `targvol` to override this behavior:

```python
seed.gentags.targvol = 10.0  # 10 Å³ per formula unit
```

## Common Seed Patterns

### Binary Compound

```python
# AB₂ compound (e.g., SiO₂, TiO₂)
seed = SeedAtoms('SiO2', cell=[5, 5, 5], pbc=True)
seed.gentags.minsep = (1.5, {'Si-O': 1.6, 'O-O': 2.0})
seed.gentags.varvol = 20
seed.gentags.symmops = (2, 4)
```

### Pure Element

```python
# Elemental system (e.g., carbon, silicon)
seed = SeedAtoms('C', cell=[3, 3, 3], pbc=True)
seed[0].num = 8
seed.gentags.minsep = 1.3
seed.gentags.symmops = (4, 12)  # Higher symmetry likely
```

### Complex Composition

```python
# Multi-element system
seed = SeedAtoms('AlNO', cell=[4, 4, 4], pbc=True)

al = seed[0]
al.num = (2, 4)    # 2-4 Al

n = seed[1]
n.num = 2          # Exactly 2 N

o = seed[2]
o.num = (2, 6)     # 2-6 O

seed.gentags.minsep = (1.4, {
    'Al-N': 1.7,
    'Al-O': 1.6,
    'N-O': 1.2
})
```

### Surface Slab

```python
# Slab with vacuum
seed = SeedAtoms('Pt8', cell=[5, 5, 15], pbc=True)
seed[0].num = 8
seed.gentags.minsep = 2.5
seed.gentags.slab = True     # Constrain to slab geometry
seed.gentags.vacuum = 10.0   # 10 Å vacuum
```

## Troubleshooting Seeds

### Generation Always Fails

If `build_random_atoms()` consistently returns `None`:

**Problem**: Constraints too tight
**Solution**: Relax minsep or increase cell size

```python
# Too tight
seed.gentags.minsep = 3.0  # For 8 atoms in 3×3×3 cell

# Better
seed.gentags.minsep = 1.8
```

**Problem**: Cell too small
**Solution**: Increase initial cell dimensions

```python
# Too small
seed = SeedAtoms('C8', cell=[2, 2, 2], pbc=True)

# Better
seed = SeedAtoms('C8', cell=[4, 4, 4], pbc=True)
```

### All Structures Too Similar

If generated structures lack diversity:

**Problem**: Too much symmetry
**Solution**: Reduce symmetry constraints

```python
# Too constrained
seed.gentags.symmops = (8, 12)

# More diverse
seed.gentags.symmops = (1, 4)
```

**Problem**: Not enough volume variation
**Solution**: Increase varvol

```python
# Too rigid
seed.gentags.varvol = 5

# More variation
seed.gentags.varvol = 30
```

### Unrealistic Structures

If generated structures are chemically unrealistic:

**Problem**: Minsep too small
**Solution**: Use element-specific separations

```python
# Generic
seed.gentags.minsep = 1.0

# Element-specific
seed.gentags.minsep = (1.2, {'Si-O': 1.5, 'O-O': 1.8})
```

## Advanced Concepts

### Formula Units

Use `nform` to control the number of formula units:

```python
# SiO₂ with 2-4 formula units = 6-12 atoms total
seed = SeedAtoms('SiO2', cell=[5, 5, 5], pbc=True)
seed.gentags.nform = (2, 4)
```

### Position Amplitudes

Control how much atoms can move from initial positions:

```python
seed = SeedAtoms('AB', cell=[4, 4, 4], pbc=True)
seed[0].posamp = 2.0   # A can move ±2 Å from initial position
seed[1].posamp = 0.5   # B stays closer to initial position
```

### Nested Ranges

Some parameters support nested ranges for sophisticated control:

```python
# MINSEP with ranges
seed.gentags.minsep = ((1.5, 2.0), {'Si-O': (1.6, 1.9)})

# POSAMP nested ranges
seed.gentags.posamp = ((1.0, 2.0), {})
```

## Best Practices

1. **Test your seed**: Generate a few structures to verify constraints work
2. **Inspect output**: Look at generated structures to check they're reasonable
3. **Iterate**: Start loose, tighten constraints based on results
4. **Document**: Keep notes on which seed parameters work well
5. **Save seeds**: Write successful seeds to files for reproducibility

```python
# Save a seed for future use
seed.write_seed('good-seed.cell')
```

## See Also

- [Creating Seeds](../how-to/create-seeds.md) - Practical examples
- [Buildcell Parameters](buildcell-parameters.md) - Complete parameter reference
- [AIRSS Overview](airss-overview.md) - Understanding the AIRSS method

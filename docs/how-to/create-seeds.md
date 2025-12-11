# How to Create Structure Seeds

Learn how to create seed files for different structure search scenarios.

## Basic Seed Creation

### Simple Composition

Create a seed with a fixed composition:

```python
from airsspy import SeedAtoms

# Fixed composition: 6 carbon and 12 oxygen atoms
seed = SeedAtoms('C6O12', cell=[5, 5, 5], pbc=True)
```

### Variable Composition

Use per-atom tags to specify variable composition:

```python
seed = SeedAtoms('C', cell=[5, 5, 5], pbc=True)
carbon = seed[0]
carbon.num = (4, 8)  # 4-8 carbon atoms will be generated
```

## Setting Global Parameters

Global parameters apply to the entire structure:

```python
seed = SeedAtoms('Si8', cell=[5, 5, 5], pbc=True)

# Basic constraints
seed.gentags.minsep = 2.0           # Minimum separation (Å)
seed.gentags.varvol = 20            # Volume variation (%)

# Symmetry control
seed.gentags.symmops = (2, 4)       # 2-4 symmetry operations
seed.gentags.compact = True         # Apply Niggli reduction

# Cell geometry
seed.gentags.targvol = 10.0         # Target volume per formula unit
```

## Per-Atom Parameters

Set parameters for individual atoms:

```python
seed = SeedAtoms('O2', cell=[5, 5, 5], pbc=True)

# Access first oxygen atom
o1 = seed[0]
o1.tagname = 'O1'        # Custom tag name
o1.posamp = 2.0          # Position randomization amplitude
o1.num = 2               # Number of this atom type

# Access second oxygen atom
o2 = seed[1]
o2.tagname = 'O2'
o2.fix = True            # Fix this atom's position
```

## MINSEP Constraints

The `minsep` parameter controls minimum atomic separations.

### Simple MINSEP

Set a global minimum separation for all atom pairs:

```python
seed.gentags.minsep = 1.5  # All pairs must be ≥1.5 Å apart
```

### Species-Specific MINSEP

Set different separations for specific atom pairs:

```python
seed = SeedAtoms('Si2O4', cell=[5, 5, 5], pbc=True)

# Format: (default_value, {pair_dict})
seed.gentags.minsep = (1.5, {'Si-O': 1.8, 'O-O': 2.0})
```

This means:
- Si-O pairs: minimum 1.8 Å
- O-O pairs: minimum 2.0 Å
- All other pairs: minimum 1.5 Å

### Range MINSEP

Use ranges to allow buildcell to vary the constraint:

```python
# Format: (default_range, {pair_ranges})
seed.gentags.minsep = ((1.5, 2.0), {'Si-O': (1.8, 2.2)})
```

Buildcell will randomly select separations within these ranges.

## Advanced Examples

### Multiple Element Types

Create a seed with different atom types and per-atom constraints:

```python
seed = SeedAtoms('AlNO', cell=[4, 4, 4], pbc=True)

# Aluminum
al = seed[0]
al.num = (2, 4)          # 2-4 Al atoms
al.posamp = 1.5          # Position variation

# Nitrogen
n = seed[1]
n.num = 2                # Fixed at 2 N atoms
n.coord = (3, 4)         # 3-4 coordination

# Oxygen
o = seed[2]
o.num = (2, 6)           # 2-6 O atoms

# Set global minsep
seed.gentags.minsep = (1.5, {
    'Al-N': 1.8,
    'Al-O': 1.7,
    'N-O': 1.2
})
```

### Symmetric Structures

Request specific symmetry operations:

```python
seed = SeedAtoms('C8', cell=[3, 3, 3], pbc=True)
seed[0].num = 8

# Request 4-8 symmetry operations
seed.gentags.symmops = (4, 8)

# Or restrict to specific crystal system
seed.gentags.system = 'cubic'

# Only symmorphic space groups
seed.gentags.symmorphic = True
```

### Cell Volume Control

Control the generated cell volume:

```python
seed = SeedAtoms('Si4', cell=[4, 4, 4], pbc=True)
seed[0].num = 4

# Target volume per formula unit
seed.gentags.targvol = 40.0          # 40 Å³ per formula unit

# Allow ±20% variation
seed.gentags.varvol = 20

# Cell shape variation
seed.gentags.cellamp = 0.5           # Randomize cell shape
```

### Position Constraints

Control where atoms can be placed:

```python
seed = SeedAtoms('Fe2', cell=[5, 5, 5], pbc=True)

fe1 = seed[0]
fe1.num = 1
fe1.fix = True          # Keep at original position

fe2 = seed[1]
fe2.num = 1
fe2.posamp = 2.0        # Allow movement within 2 Å
fe2.xamp = 1.0          # X direction: ±1 Å
fe2.yamp = 1.0          # Y direction: ±1 Å
fe2.zamp = 0.5          # Z direction: ±0.5 Å
```

## Exporting Seeds

### Write to File

Save the seed to a .cell file:

```python
seed.write_seed('myseed.cell')
```

### Get as String

View the seed content:

```python
lines = seed.get_cell_inp_lines()
print('\n'.join(lines))
```

Example output:
```
%BLOCK lattice_cart
5.0000000000  0.0000000000  0.0000000000
0.0000000000  5.0000000000  0.0000000000
0.0000000000  0.0000000000  5.0000000000
%ENDBLOCK lattice_cart
%BLOCK positions_abs
Si  0.0000000000 0.0000000000 0.0000000000 # Si0 % NUM=4
%ENDBLOCK positions_abs
#MINSEP=2.0
#VARVOL=20
```

## Complete Example

Here's a complete example for a silicon-oxygen system:

```python
from airsspy import SeedAtoms

# Create seed
seed = SeedAtoms('SiO2', cell=[4, 4, 4], pbc=True)

# Configure silicon
si = seed[0]
si.num = 2              # 2 Si atoms
si.tagname = 'Si_bulk'

# Configure oxygen
o = seed[1]
o.num = 4               # 4 O atoms
o.tagname = 'O_bulk'

# Global parameters
seed.gentags.minsep = (1.5, {'Si-O': 1.6, 'O-O': 2.0})
seed.gentags.varvol = 15
seed.gentags.symmops = (2, 4)
seed.gentags.compact = True

# Save to file
seed.write_seed('SiO2-search.cell')

# Generate a structure
random_structure = seed.build_random_atoms(timeout=15)
```

## Tips and Best Practices

1. **Start simple**: Begin with basic constraints and add complexity as needed
2. **Test constraints**: Generate a few structures to verify constraints work
3. **Adjust minsep**: If generation fails, check if minsep is too restrictive
4. **Volume matters**: Make sure initial cell isn't too small/large for your composition
5. **Use symmetry wisely**: More symmetry operations = more constrained searches

## See Also

- [Generating Structures](generate-structures.md) - Different generation approaches
- [Buildcell Parameters](../explanation/buildcell-parameters.md) - Complete parameter reference
- [API: SeedAtoms](../reference/index.md) - Full API documentation

# Buildcell Parameters Reference

Complete reference for buildcell generation parameters accessible through airsspy.

## Overview

Buildcell parameters control how random structures are generated. In airsspy, these are set through:

- **Global parameters**: `seed.gentags` (applies to entire structure)
- **Per-atom parameters**: `seed[i]` (applies to specific atoms)

## Global Parameters (gentags)

Access via `seed.gentags.parameter_name`.

### Cell Geometry

| Parameter | Type | Description |
|-----------|------|-------------|
| `targvol` | float or tuple | Target volume per formula unit (Ų) |
| `varvol` | float | Volume variation percentage |
| `cellamp` | float | Cell shape variation amplitude |
| `compact` | bool | Apply Niggli reduction to cell |

Example:
```python
seed.gentags.targvol = 40.0   # 40 ų per formula unit
seed.gentags.varvol = 20      # ±20% variation
seed.gentags.cellamp = 0.5    # Randomize cell shape
seed.gentags.compact = True   # Apply Niggli reduction
```

### Symmetry Control

| Parameter | Type | Description |
|-----------|------|-------------|
| `symmops` | int or tuple | Number of symmetry operations |
| `sgrank` | int or tuple | Minimum spacegroup rank |
| `symmorphic` | bool | Restrict to symmorphic groups |
| `system` | str | Crystal system (e.g., 'cubic', 'hexagonal') |
| `symmno` | int | Specific space group number (1-230) |
| `symm` | str | Space group symbol (e.g., 'Fm-3m') |

Example:
```python
seed.gentags.symmops = (2, 4)     # 2-4 symmetry operations
seed.gentags.sgrank = 3           # Minimum rank 3
seed.gentags.symmorphic = True    # Only symmorphic
seed.gentags.system = 'cubic'     # Cubic system
```

### Atomic Separation

| Parameter | Type | Description |
|-----------|------|-------------|
| `minsep` | float, tuple, or nested | Minimum separation constraints |
| `overlap` | float | Maximum allowed atomic overlap |

The `minsep` parameter supports multiple formats:

```python
# Simple: global value for all pairs
seed.gentags.minsep = 1.5

# With species pairs: (default, {pair_dict})
seed.gentags.minsep = (1.5, {'Si-O': 1.8, 'O-O': 2.0})

# With ranges: ((min, max), {pair_ranges})
seed.gentags.minsep = ((1.5, 2.0), {'Si-O': (1.8, 2.2)})
```

### Position Variation

| Parameter | Type | Description |
|-----------|------|-------------|
| `posamp` | float or nested | Position randomization amplitude |
| `minamp` | float or tuple | Minimum positional amplitude |
| `zamp` | float or tuple | Z-direction amplitude |
| `xamp` | float or tuple | X-direction amplitude |
| `yamp` | float or tuple | Y-direction amplitude |
| `angamp` | float or tuple | Angular amplitude |

Example:
```python
seed.gentags.posamp = 2.0    # ±2 Å randomization
seed.gentags.zamp = 0.5      # ±0.5 Å in Z direction
```

### Composition Control

| Parameter | Type | Description |
|-----------|------|-------------|
| `nform` | int or tuple | Number of formula units |
| `natom` | int or tuple | Total number of atoms |
| `species` | str | Species list |
| `formula` | str | Chemical formula |

Example:
```python
seed.gentags.nform = (2, 4)      # 2-4 formula units
seed.gentags.natom = 8           # Exactly 8 atoms
```

### Advanced Options

| Parameter | Type | Description |
|-----------|------|-------------|
| `adjgen` | int | Adjust general positions |
| `breakamp` | float | Symmetry breaking amplitude |
| `celladapt` | bool | Enable cell adaptation |
| `cellcon` | any | Cell constraints |
| `cons` | str | Additional constraints |
| `cylinder` | any | Cylindrical constraint |
| `flip` | bool | Enable mirror reflection |
| `focus` | any | Focus on specific species |
| `maxbangle` | float | Maximum bond angle |
| `maxtime` | float | Maximum time for attempts |
| `minbangle` | float | Minimum bond angle |
| `molecules` | any | Molecular fragments |
| `nocompact` | bool | Disable compaction |
| `nopush` | bool | Disable atomic pushing |
| `octet` | any | Octet rule constraints |
| `permfrac` | float | Permutation fraction |
| `permute` | any | Species to permute |
| `rad` | any | Atomic radii |
| `rash` | any | Relax-and-shake option |
| `rash_angamp` | float | RASH angular amplitude |
| `rash_posamp` | float | RASH position amplitude |
| `remove` | any | Remove species |
| `seed` | int | Random seed |
| `slab` | any | Slab geometry |
| `slack` | float or tuple | Slack factor |
| `species` | str | Species specification |
| `sphere` | any | Spherical constraint |
| `spin` | any | Magnetic spin |
| `supercell` | any | Supercell expansion |
| `surface` | any | Surface constraint |
| `three` | any | Three-body constraints |
| `tight` | bool | Tight packing |
| `vacancies` | any | Vacancy creation |
| `vacuum` | float | Vacuum thickness |
| `width` | float | Width parameter |

## Per-Atom Parameters

Access via `seed[index].parameter_name`.

### Basic Properties

| Parameter | Type | Description |
|-----------|------|-------------|
| `tagname` | str | Custom tag name for this atom |
| `num` | int or tuple | Number of this atom type |
| `spin` | float | Magnetic spin moment |
| `occ` | float | Occupation factor |
| `rad` | float | Atomic radius |
| `vol` | float | Atomic volume |
| `mult` | int | Multiplicity |

Example:
```python
atom = seed[0]
atom.tagname = 'Fe_center'
atom.num = (2, 4)        # 2-4 atoms of this type
atom.spin = 2.5          # Magnetic moment
```

### Position Control

| Parameter | Type | Description |
|-----------|------|-------------|
| `posamp` | float or tuple | Position amplitude |
| `minamp` | float or tuple | Minimum amplitude |
| `zamp` | float or tuple | Z-direction amplitude |
| `xamp` | float or tuple | X-direction amplitude |
| `yamp` | float or tuple | Y-direction amplitude |
| `angamp` | float or tuple | Angular amplitude |

Example:
```python
atom = seed[0]
atom.posamp = 2.0        # General position variation
atom.xamp = 1.0          # X: ±1 Å
atom.yamp = 1.0          # Y: ±1 Å
atom.zamp = 0.5          # Z: ±0.5 Å
```

### Constraints

| Parameter | Type | Description |
|-----------|------|-------------|
| `fix` | bool | Fix position completely |
| `nomove` | bool | Fix during generation only |
| `perm` | bool | Randomly permute species |
| `adatom` | bool | Add after supercell creation |
| `athole` | bool | Place at hole position |

Example:
```python
atom = seed[0]
atom.fix = True          # Fix this atom
atom.nomove = False      # Allow movement
atom.perm = False        # Don't permute
```

### Coordination

| Parameter | Type | Description |
|-----------|------|-------------|
| `coord` | int or tuple | Target coordination number |

Example:
```python
atom = seed[0]
atom.coord = (3, 4)      # 3-4 coordinated
```

## Usage Examples

### Basic Structure Search

```python
from airsspy import SeedAtoms

seed = SeedAtoms('Si4', cell=[4, 4, 4], pbc=True)
seed[0].num = 4

# Global constraints
seed.gentags.minsep = 2.0
seed.gentags.varvol = 20
seed.gentags.symmops = (2, 4)
```

### Multi-Element System

```python
seed = SeedAtoms('SiO2', cell=[5, 5, 5], pbc=True)

# Per-atom control
si = seed[0]
si.num = 2
si.coord = (4, 6)

o = seed[1]
o.num = 4
o.coord = (2, 3)

# Species-specific separations
seed.gentags.minsep = (1.5, {'Si-O': 1.6, 'O-O': 2.0})
```

### High-Symmetry Search

```python
seed = SeedAtoms('C8', cell=[3, 3, 3], pbc=True)
seed[0].num = 8

seed.gentags.system = 'cubic'
seed.gentags.symmops = (8, 12)
seed.gentags.compact = True
```

### Surface Slab

```python
seed = SeedAtoms('Pt8', cell=[5, 5, 15], pbc=True)
seed[0].num = 8

seed.gentags.slab = True
seed.gentags.vacuum = 10.0
seed.gentags.minsep = 2.5
```

## Parameter Ranges

Many parameters accept ranges as tuples:

```python
# Single value
seed.gentags.symmops = 4

# Range (randomly selected)
seed.gentags.symmops = (2, 4)

# Per-atom
atom.num = (4, 8)
atom.coord = (3, 4)
```

## Tips

1. **Start simple**: Use basic constraints first (minsep, varvol)
2. **Test constraints**: Generate structures to verify they work
3. **Adjust iteratively**: Tighten constraints based on results
4. **Use chemistry**: Set realistic separations based on expected bonds
5. **Consider symmetry**: More symmetry = faster DFT but less diversity

## See Also

- [Creating Seeds](../how-to/create-seeds.md) - Practical examples
- [Seed Concept](seed-concept.md) - Understanding seeds
- [Quickstart Tutorial](../getting-started/quickstart.md) - Complete workflow

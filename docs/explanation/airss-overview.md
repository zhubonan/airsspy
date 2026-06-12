# AIRSS Overview

Learn about AIRSS (Ab initio Random Structure Searching) and how airsspy provides a Python interface to this powerful method.

## What is AIRSS?

AIRSS (Ab initio Random Structure Searching) is a structure prediction method that combines random structure generation with ab initio relaxation. It was developed at Cambridge University and is widely used for crystal structure prediction.

### The AIRSS Workflow

1. **Define Search Space**: Create a seed file specifying composition, constraints, and symmetry
2. **Generate Structures**: Use `buildcell` to create random structures within constraints
3. **Relax Structures**: Optimize each structure with ab initio or empirical methods
4. **Analyze Results**: Identify low-energy structures and unique minima

### Key Advantages

- **Simple conceptually**: Easy to understand and implement
- **Flexible**: Works with any relaxation method (DFT, force fields, etc.)
- **Robust**: Finds structures even with minimal prior information
- **Parallelizable**: Each structure is independent
- **Well-tested**: Used in numerous successful structure prediction studies

## What is airsspy?

airsspy is a Python package that makes AIRSS accessible through the Atomic Simulation Environment (ASE). It provides:

- **Modern Python Interface**: Clean, typed API for AIRSS workflows
- **ASE Integration**: Seamless integration with ASE calculators and optimizers
- **Buildcell Wrapper**: Simplified access to AIRSS buildcell functionality
- **RES File Support**: Read and write CASTEP .res output files
- **Pythonic Workflows**: Use AIRSS with Python scripts and Jupyter notebooks

## The buildcell Program

`buildcell` is the core AIRSS executable that generates random structures. It:

- Takes a seed file as input
- Applies random variations (positions, cell, symmetry)
- Enforces constraints (minimum separations, symmetry operations)
- Outputs a randomized structure

airsspy wraps `buildcell` subprocess calls, making it accessible from Python:

```python
from airsspy import SeedAtoms

seed = SeedAtoms('Si4', cell=[4, 4, 4], pbc=True)
seed[0].num = 4
seed.gentags.minsep = 2.0

# This calls buildcell under the hood
atoms = seed.build_random_atoms()
```

## Seed Files

A seed file defines the search template. It specifies:

- **Composition**: Elements and number of atoms
- **Cell**: Initial lattice (will be randomized)
- **Constraints**: Minimum separations, volume ranges
- **Symmetry**: Number and type of symmetry operations

Example seed file:
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
```

In airsspy, you create seeds programmatically:

```python
seed = SeedAtoms('Si', cell=[4, 4, 4], pbc=True)
seed[0].num = 4              # 4 Si atoms
seed.gentags.minsep = 2.0    # Minimum 2 Å separation
seed.gentags.varvol = 20     # ±20% volume variation
```

See [Seed Concept](seed-concept.md) for more details.

## Integration with ASE

airsspy extends ASE's Atoms class, making it compatible with ASE's ecosystem:

### ASE Calculators

Use any ASE calculator for relaxation:

```python
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from ase.calculators.vasp import Vasp

# Example with Lennard-Jones
atoms.set_calculator(LennardJones())
```

### ASE Optimizers

Use ASE optimization algorithms:

```python
from ase.optimize import BFGS, FIRE, BFGSLineSearch
from ase.filters import UnitCellFilter

# Optimize both positions and cell
opt = BFGS(UnitCellFilter(atoms))
opt.run(fmax=0.05)
```

### ASE I/O

Read and write using ASE formats:

```python
from ase.io import read, write

write('structure.xyz', atoms)
write('structure.cif', atoms)
write('POSCAR', atoms)
```

## Typical AIRSS-ASE Workflow

Here's a complete structure search workflow:

```python
from airsspy import SeedAtoms, save_airss_res
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.filters import UnitCellFilter

# 1. Define seed
seed = SeedAtoms('Al8', cell=[3, 3, 3], pbc=True)
seed[0].num = 8
seed.gentags.minsep = 1.5

# 2. Set up calculator
calc = EMT()

# 3. Search loop
results = []
for i in range(20):
    # Generate random structure
    atoms = seed.build_random_atoms()
    if atoms is None:
        continue

    # Relax structure
    atoms.set_calculator(calc)
    opt = BFGS(UnitCellFilter(atoms), logfile=None)
    opt.run(fmax=0.05)

    results.append(atoms)

# 4. Analyze and save
results.sort(key=lambda a: a.get_potential_energy())

for i, atoms in enumerate(results):
    info = {
        'uid': f'Al-{i}',
        'P': 0.0,
        'V': atoms.get_volume(),
        'H': atoms.get_potential_energy(),
        'nat': len(atoms),
        'sym': 'P1'
    }
    save_airss_res(atoms, info, f'Al-{i}.res')

print(f"Lowest energy: {results[0].get_potential_energy():.4f} eV")
```

## Further Reading

- **Official AIRSS Documentation**: [http://www.mtg.msm.cam.ac.uk/Codes/AIRSS](http://www.mtg.msm.cam.ac.uk/Codes/AIRSS)
- **Key Papers**:
  - Pickard & Needs, Phys. Rev. Lett. 97, 045504 (2006)
  - Pickard & Needs, J. Phys.: Condens. Matter 23, 053201 (2011)
- **ASE Documentation**: [https://wiki.fysik.dtu.dk/ase/](https://wiki.fysik.dtu.dk/ase/)

## Next Steps

- Learn more about [Seed Files](seed-concept.md)
- Explore [Buildcell Parameters](buildcell-parameters.md)
- Try the [Quickstart Tutorial](../getting-started/quickstart.md)

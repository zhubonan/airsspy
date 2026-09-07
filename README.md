airsspy
---------

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/zhubonan/airsspy/HEAD)
[![Documentation Status](https://readthedocs.org/projects/airsspy/badge/?version=latest)](https://airsspy.readthedocs.io/en/latest/?badge=latest)

A package to help work with Ab initio Random Structure Searching ([AIRSS](https://www.mtg.msm.cam.ac.uk/Codes/AIRSS))
using the Atomic Simulation Environment ([ASE](https://wiki.fysik.dtu.dk/ase/)).
It supports building search seeds for AIRSS in interactive Python environments,
running local or jobflow-backed searches, ranking results, and converting
between AIRSS `.res` and extxyz files.
One of the important steps for performing a successful search with AIRSS is to have a sensible seed for generating
random structures, which are subsequently relaxed using the method of choice.
In general, AIRSS only relies on a few simple parameters to generate random structures, such as numbers of atoms,
numbers of species, and cell volume.
However, for complicated search involving surfaces and/or interfaces, hand-building seed files becomes a
tedious or impossible job to do.
ASE has a suite of tools for manipulating atomic structures which can be very helpful for building structures,
and here, for building search seeds.

AIRSS is an open source code licensed under GPLv2;
this package does not contain any source code of AIRSS nor links to it.

Documentation
-------------

Full documentation is available at [airsspy.readthedocs.io](https://airsspy.readthedocs.io/).

What this does
--------------
* Allows preparing AIRSS seeds using ASE's `Atoms` interface
* Wraps the external `buildcell` executable
* Reads, writes, ranks, packs, unpacks, and converts AIRSS `.res` files
* Provides the `ap` CLI for local searches, relaxations, single-points,
  jobflow deployment, and database queries

Try interactively
-----------------

Interactive jupyter notebook examples can be found in the `examples` folder.
Click the *binder* badge above to launch these examples in a pre-built environment and try it in your browser!

Dependencies
-----------

* [ase](https://wiki.fysik.dtu.dk/ase/): The atomic simulation environment
* [castepinput](https://gitlab.com/bz1/castepinput): A lightweight writer/reader for the input files of [CASTEP](https://www.castep.org/).

Installation
-----------

The package can be installed from PyPI together with the dependencies:

```
pip install airsspy
```

Alternatively, install directly from the repository:

```
pip install git+https://github.com/zhubonan/airsspy
```

Usage
-----
Assuming you are familiar with `ase`, Python, and have some basic knowledge of AIRSS.
To prepare a seed for generating a *sensible* random structure:

```python
from airsspy import SeedAtoms
seed = SeedAtoms('C6')
seed.build.varvol = 20
seed.build.symmops = (2, 4)

# Can also access per-atom tags/keywords just like in ASE
for i in range(0, 6, 2):
    atom = seed[i]
    atom.tagname = 'CX'
    atom.posamp = 2
```

To write the seed file onto the disk:

```python
seed.write_seed('C6.cell')
# With IPython
# Use the buildcell executable to generate the file
!buildcell < C6.cell > C6-rand.cell
```

To generate a cell we can create a `Buildcell` instance,
which is a helper wrapper around the AIRSS `buildcell` program:

```python
from airsspy import Buildcell
buildcell = Buildcell(seed)
random_atoms = buildcell.generate()
```

A shortcut is also available as a method of `SeedAtoms`:

```python
random_atoms = seed.build_random_atoms()
```

Command-line interface
----------------------

The installed command is `ap`:

```bash
ap --help
ap check airss
ap run search --seed Si --code castep --nmax 100 --pack
ap run relax --cell "*.res" --code eddp --calculator /path/to/model.json \
  --eddp-project /path/to/EDDPotentials.jl
ap rank packed.res -t 20 -de 0.05
ap convert packed.res structures.xyz
```

See the command-line reference in the documentation for the full interface.

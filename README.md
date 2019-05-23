airsspy
---------
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gl/bz1%2Fairsspy/master)

Package to help working with the Ab initio Random Structure Searching ([AIRSS](https://www.mtg.msm.cam.ac.uk/Codes/AIRSS))
using Atomic simulation environment ([ase](https://wiki.fysik.dtu.dk/ase/)).
AIRSS is a open source code licensed under GPLv2, 
this package does not contain any source code of AIRSS nor links to it.


What this does
--------------
* Allow preparing seed for AIRSS using ASE's `atoms` interface
* Allow ase's calculators to be used in AIRSS to do relaxations

Try interactively
-----------------
Interactive jupyter notebook examples can be found in the `examples` folder.
Click the *binder* badge above to launch these examples in a pre-built environment and try it in your browser!

Requirements
-----------
* [ase](https://wiki.fysik.dtu.dk/ase/): The atomic simulation environment
* [castepinput](https://gitlab.com/bz1/castepinput): A light weight writer/reader for the input files of [CASTEP](www.caste.org).

Usage
-----
Assuming you are familiar with `ase`, python and has some basic knowledge of AIRSS.
To prepare a seed for generating a *sensible* random structure:
```
from airsspy import SeedAtoms
seed = SeedAtoms('C6')
seed.buiid.varvol = 20
seed.build.symmops = (2, 4)

# Can also access per `atom` tags/ketwords just like in ASE
for i in range(0, 6, 2):
    atom = seed[i]
    atom.tagname = 'CX'
    atom.posamp = 2
```

To write the seed file onto the disk:
```
atoms.write_seed('C6.cell')
# With IPython
# Use the buildcell executable to generate the file
!buildcell < C6.cell > C6-rand.cell
```

To generate a cell we can create a `Buildcell` instance,
which is helping wrapper to the `buildcell` program of AIRSS:

```
from airsspy import Buildcell
buidcell = Buildcell(seed)
random_atoms = builcell.generate()
```

A shortcut is also available as an method of the `SeedAtoms`:
```
random_atoms = seed.build_random_atoms()
```

Limitations
-----------
Due to the lack of `timeout` argument of `Popen.communicate` in python 2.7,
communication with the `buildcell` is not available. Hence, direct generation and 
retrieval of the random structure are not supported in python. However, it is 
still possible to write the seed out and call the program externally.

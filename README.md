airss-ase
---------

Package to help working with Ab initio Random Structure Searching ([AIRSS](https://www.mtg.msm.cam.ac.uk/Codes/AIRSS))
using Atomic simulation environment ([ase](https://wiki.fysik.dtu.dk/ase/))


What this does
--------------
* Allow preparing seed for AIRSS using ASE's `atoms` interface
* Allow ase's calculators to be used in AIRSS to do relaxations


Requirement
-----------
* [ase](https://wiki.fysik.dtu.dk/ase/): The atomic simulation environment
* [castepinput](https://gitlab.com/bz1/castepinput): A light weight writer/reader for the input files of [CASTEP](www.caste.org).

Usage
-----
Prepare a seed for generation:
```
from airss_ase import TemplateAtoms
seed = TemplateAtoms('C6')
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
which is an wrapp to the `buildcell` program of AIRSS:

```
from airss_ase import Buildcell
buidcell = Buildcell(atoms=seed)
random_atoms = builcell.generate()
```

A shortcut is also avaliable as an method of the `TemplateAtoms`:
```
random_atoms = seed.get_random_atoms()
```

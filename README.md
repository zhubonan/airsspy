airss-ase
---------

Package to help working with Ab initio Random Structure Searching ([AIRSS](https://www.mtg.msm.cam.ac.uk/Codes/AIRSS))
using Atomic simulation environment ([ase](https://wiki.fysik.dtu.dk/ase/))


What this does
--------------
* Allow preparing seed for AIRSS using ASE's `atoms` interface
* Allow ase's calculators to be used in AIRSS to do relaxations


Usage
-----

Prepare a seed for generation
```
from airss_ase import TemplateAtom
seed = TemplateAtom('C6')
seed.buiid.varvol = 20
seed.build.symmops = (2, 4)
seed.build.nato = (2, 4)
```

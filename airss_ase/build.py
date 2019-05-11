"""
Module for building the random cell
"""

from castepinput import CellInput
from .seed import write_buildcell_seed


class Buildcell:
    """
    File based inteface to the buildcell program which is part of
    Ab inito Radnom Structure Searhcing (AIRSS) package
    """

    def __init__(self, atoms):
        """Initialise an Buildcell object"""
        self.atoms = atoms
        self.proc = None
        self.res = None
        self.bc_out = None
        self.bc_err = None

    def generate(self, timeout=10, write_cell=None):
        """Generate a random atom based on a template
        timeout: time to wait for buildcell binary
        write_seed : Name of the output cell to be written"""
        import io
        import ase.io.castep
        import subprocess as sbp
        tmp = io.StringIO()
        write_buildcell_seed(tmp, self.atoms)
        bc = sbp.Popen('buildcell',
                       universal_newlines=True,
                       stdin=sbp.PIPE,
                       stdout=sbp.PIPE,
                       stderr=sbp.PIPE)
        self.proc = bc
        tmp.seek(0)
        self.seed = tmp.read()
        try:
            self.bc_out, self.bc_err = bc.communicate(input=self.seed,
                                                      timeout=timeout)
        except sbp.TimeoutExpired:
            bc.kill()
            self.bc_out, self.bc_err = bc.communicate()
            print('Generation Failed to finished. Output captured')
            return

        tmp.close()
        tmp = io.StringIO()
        tmp.write(self.bc_out)
        tmp.seek(0)
        res = ase.io.castep.read_castep_cell(tmp)

        # Write the output from buildcell
        if write_cell:
            with open(write_cell + '.cell', 'w') as output:
                tmp.seek(0)
                for line in tmp:
                    output.write(line)
        tmp.close()

        # Detach unnecessary calculator
        res.calc = None
        # Store as attribute
        self.res = res
        return res

    def write_seed(self, seedname):
        """Write the seed for buildcell to the disk"""
        with open(seedname + '.cell', 'w') as tmp:
            write_buildcell_seed(tmp, self.atoms)

    def gen_and_view(self, viewer=None, wrap=False, timeout=20):
        from ase.visualize import view
        res = self.generate(timeout=timeout)
        if not res:
            return
        if wrap:
            res.wrap()
        if viewer:
            view(res, viewer=viewer)
        else:
            view(res)

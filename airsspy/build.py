# -*- coding: utf-8 -*-
###########################################################################
# airss-ase                                                               #
# Copyright (C) 2019  Bonan Zhu                                           #
#                                                                         #
# This program is free software; you can redistribute it and/or modify    #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation; either version 2 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# This program is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License along #
# with this program; if not, write to the Free Software Foundation, Inc., #
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             #
###########################################################################
"""
Module for building the random cell
"""

import subprocess as sbp

from castepinput import CellInput
from castepinput.parser import PlainParser
from ase import Atoms

from .common import BuildcellError


class Buildcell:
    """
    File based inteface to the buildcell program which is part of
    Ab inito Radnom Structure Searhcing (AIRSS) package
    """

    def __init__(self, atoms):
        """Initialise an Buildcell object"""
        self.atoms = atoms
        self.proc = None
        self.res_atoms = None
        # Input and output from the buildcell program
        self.bc_out = None
        self.bc_err = None
        self.bc_in = None

    def generate(self, timeout=10, write_cell=None):
        """Generate a random atom based on a template
        timeout: time to wait for buildcell binary
        write_seed : Name of the output cell to be written"""
        bc_proc = sbp.Popen(
            'buildcell',
            universal_newlines=True,
            stdin=sbp.PIPE,
            stdout=sbp.PIPE,
            stderr=sbp.PIPE)
        self.proc = bc_proc
        cell = '\n'.join(self.atoms.get_cell_inp_lines())
        self.bc_in = cell
        try:
            self.bc_out, self.bc_err = bc_proc.communicate(
                input=cell, timeout=timeout)
        except sbp.TimeoutExpired:
            bc_proc.kill()
            self.bc_out, self.bc_err = bc_proc.communicate()
            raise BuildcellError
            return None
        else:
            bc_proc.kill()

        # Write the output from buildcell
        if write_cell:
            with open(write_cell + '.cell', 'w') as output:
                output.write(self.bc_out)

        # Process the result
        outlines = self.bc_out.split('\n')
        parser = PlainParser(outlines)
        parser.parse()
        cellout = CellInput()
        for k, value in parser.get_dict().items():
            cellout.__setitem__(k, value)

        cell = cellout.get_cell()
        elements, positions, _ = cellout.get_positions()
        atoms = Atoms(symbols=elements, cell=cell, positions=positions)
        self.res_atoms = atoms
        return atoms

    def write_seed(self, seedname):
        """Write the seed for buildcell to the disk"""
        self.atoms.write_seed(seedname)

    def gen_and_view(self, viewer=None, wrap=False, timeout=20):
        """Geneate one and view with viewer immediately. Wrap if needed."""
        from ase.visualize import view
        atoms = self.generate(timeout=timeout)
        if not atoms:
            return
        if wrap:
            atoms.wrap()
        if viewer:
            view(atoms, viewer=viewer)
        else:
            view(atoms)

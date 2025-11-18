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
import pytest
from airsspy.seed import SeedAtom, SeedAtoms, SeedAtomTag, BuildcellParam, tuple2range


def test_bc_param():
    bcp = BuildcellParam()
    bcp.fix = True
    bcp.nform = 3
    assert "NFORM" in bcp.to_string()
    bcp.minsep = [2, {"Ce-O": (2, 3)}]
    assert "Ce-O=2-3" in bcp.to_string()


def test_nested_range():
    bcp = BuildcellParam()

    # Test different possible configurations
    bcp.minsep = 2
    assert "MINSEP=2" in bcp.to_string()

    bcp.minsep = (2, 3)
    assert "MINSEP=2-3" in bcp.to_string()

    bcp.minsep = (2, {"Ce-O": 1})
    assert "MINSEP=2 Ce-O=1" in bcp.to_string()

    bcp.minsep = ((2, 3), {"Ce-O": 1})
    assert "MINSEP=2-3 Ce-O=1" in bcp.to_string()

    bcp.minsep = (2, {"Ce-O": (1, 2)})
    assert "MINSEP=2 Ce-O=1-2" in bcp.to_string()


def test_atom_param():
    """
    Test the single atom specification
    """

    param = SeedAtomTag()
    param.num = 3
    param.posamp = (1, 2)
    param.tagname = "O1"
    string = param.to_string()
    assert string.startswith("# O1 %")
    assert "POSAMP=1-2" in string
    assert "NUM=3" in string


def test_template_atom():
    ta = SeedAtom(symbol="C")
    ta.xamp = 0
    ta.tagname = "C1"
    assert ta.to_string() == "# C1 % XAMP=0"


def test_template_atom_from_tmp():
    bc = SeedAtoms(symbols="C2")
    c1 = bc[0]
    c1.posamp = 1
    assert c1.to_string() == "# C0 % POSAMP=1"
    assert bc.arrays["atom_gentags"][0].to_string() == "# C0 % POSAMP=1"

    cell = bc.get_cell_inp()
    assert "positions_abs" in cell


def test_tuple2range():
    """
    Test tuple/number to string conversion
    """
    assert tuple2range(2) == "2"
    assert tuple2range((2, 3)) == "2-3"


def test_missing_keywords():
    """Test newly added missing keywords"""
    bcp = BuildcellParam()

    # Test new global keywords
    bcp.formula = "Si2O4"
    bcp.seed = "test123"
    bcp.vol = (100, 200)
    bcp.nfails = 10
    bcp.hole = 1.5
    bcp.holepos = "0.0 0.0 0.5"
    bcp.shift = "0.1 0.1 0.1"
    bcp.permute = "Si O"

    output = bcp.to_string()
    assert "FORMULA=Si2O4" in output
    assert "SEED=test123" in output
    assert "VOL=100-200" in output
    assert "NFAILS=10" in output
    assert "HOLE=1.5" in output
    assert "HOLEPOS=0.0 0.0 0.5" in output
    assert "SHIFT=0.1 0.1 0.1" in output
    assert "PERMUTE=Si O" in output


def test_new_per_atom_keywords():
    """Test newly added per-atom keywords"""
    param = SeedAtomTag()
    param.tagname = "Si1"
    param.vol = 20.0
    param.mult = 2
    param.spin = 1.5
    param.athole = True

    output = param.to_string()
    assert output.startswith("# Si1 %")
    assert "VOL=20.0" in output
    assert "MULT=2" in output
    assert "SPIN=1.5" in output
    assert "ATHOLE" in output


def test_formula_vs_species_conflict():
    """Test that FORMULA and SPECIES can coexist in airsspy"""
    bcp = BuildcellParam()
    bcp.formula = "C2H4"
    bcp.species = "C H"

    output = bcp.to_string()
    assert "FORMULA=C2H4" in output
    assert "SPECIES=C H" in output


def test_template_atom_write():
    """
    Test writing output from a template atom
    """
    from tempfile import mkdtemp
    import os

    atoms = SeedAtoms("C4")
    atoms.gentags.symmops = (2, 3)
    atoms.gentags.sgrank = 2
    atoms.gentags.minsep = [2, {"C-C": (2, 3)}]
    atoms[0].posamp = 3
    atoms[0].xamp = (2, 3)
    atoms[1].tagname = "CX0"
    buff = "\n".join(atoms.get_cell_inp_lines())

    assert "#SGRANK=2" in buff
    assert "#SYMMOPS=2-3" in buff
    assert "#MINSEP=2 C-C=2-3" in buff
    assert "XAMP=2-3" in buff
    assert "# C0" in buff
    assert "# CX0 " in buff

    # Test actually writing the files
    tmp = mkdtemp()
    fpath = os.path.join(tmp, "test.cell")
    atoms.write_seed(fpath)

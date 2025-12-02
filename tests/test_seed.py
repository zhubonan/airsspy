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


# ============================================================================
# PHASE 1: Comprehensive Tests (added for dataclass modernization)
# ============================================================================


def test_all_buildcell_tag_properties():
    """Test all tag properties in BuildcellParam"""
    bcp = BuildcellParam()

    # Test tag properties (boolean-like)
    tag_properties = [
        "fix",
        "abfix",
        "autoslack",
        "flip",
        "surface",
        "symmorphic",
        "cfix",
        "cluster",
        "tight",
        "compact",
    ]

    for prop_name in tag_properties:
        # Test setting to True
        setattr(bcp, prop_name, True)
        assert getattr(bcp, prop_name) is True, f"{prop_name} should be True"

        # Test setting to False should remove it
        setattr(bcp, prop_name, False)
        assert getattr(bcp, prop_name) is None, f"{prop_name} should be None after False"


def test_all_buildcell_range_properties():
    """Test all range properties in BuildcellParam"""
    bcp = BuildcellParam()

    # Test range properties (number or 2-tuple)
    range_properties = [
        "coord",
        "nform",
        "symmops",
        "minamp",
        "zamp",
        "xamp",
        "yamp",
        "angamp",
        "sgrank",
        "slack",
        "overlap",
        "natom",
        "vol",
        "targvol",
    ]

    for prop_name in range_properties:
        # Test single number
        setattr(bcp, prop_name, 5)
        assert getattr(bcp, prop_name) == 5, f"{prop_name} should be 5"

        # Test tuple
        setattr(bcp, prop_name, (2, 10))
        assert getattr(bcp, prop_name) == (2, 10), f"{prop_name} should be (2, 10)"

        # Test validation - tuple with wrong length should fail
        with pytest.raises(ValueError, match="two element"):
            setattr(bcp, prop_name, (1, 2, 3))

        # Test validation - tuple with non-numeric should fail
        with pytest.raises(ValueError, match="number"):
            setattr(bcp, prop_name, (1, "abc"))


def test_all_buildcell_generic_properties():
    """Test all generic properties in BuildcellParam"""
    bcp = BuildcellParam()

    # Test generic properties (any value)
    generic_properties = {
        "adjgen": "test_value",
        "breakamp": 1.5,
        "celladapt": "value",
        "cellamp": 2.0,
        "cellcon": "constraint",
        "cylinder": 3.0,
        "maxbangle": 180,
        "maxtime": 100,
        "minbangle": 0,
        "focus": "composition",
        "molecules": "H2O",
        "nocompact": True,
        "nopush": True,
        "octet": "value",
        "permfrac": 0.5,
        "permute": "Si O",
        "rad": 1.5,
        "rash": "value",
        "rash_angamp": 1.0,
        "rash_posamp": 0.5,
        "remove": "H",
        "slab": "001",
        "species": "C H O",  # Generic type per user requirement
        "sphere": 5.0,
        "spin": 2.0,
        "supercell": "2x2x2",
        "symm": "P1",
        "symmno": 1,
        "system": "cubic",
        "three": "value",
        "vacancies": "Si",
        "vacuum": 10.0,
        "width": 5.0,
        "varvol": 100.0,
        "cons": "constraint",
        "formula": "Si2O4",
        "seed": "12345",
        "nfails": 10,
        "hole": 1.5,
        "holepos": "0 0 0",
        "shift": "0.1 0.1 0.1",
    }

    for prop_name, test_value in generic_properties.items():
        setattr(bcp, prop_name, test_value)
        assert getattr(bcp, prop_name) == test_value, f"{prop_name} should be {test_value}"


def test_all_buildcell_nested_range_properties():
    """Test all nested range properties in BuildcellParam"""
    bcp = BuildcellParam()

    # Test minsep (nested range)
    bcp.minsep = 2
    assert bcp.minsep == 2

    bcp.minsep = (2, 3)
    assert bcp.minsep == (2, 3)

    bcp.minsep = (2, {"Ce-O": (1, 2)})
    assert bcp.minsep == (2, {"Ce-O": (1, 2)})

    # Test species (nested range) - note: not posamp, which is range
    bcp.species = ((1, 2), {"Si": 3})
    assert bcp.species == ((1, 2), {"Si": 3})


def test_all_seedatomtag_properties():
    """Test all properties in SeedAtomTag"""
    tag = SeedAtomTag()

    # Generic property
    tag.tagname = "C1"
    assert tag.tagname == "C1"

    # Range properties
    tag.posamp = 1.5
    assert tag.posamp == 1.5
    tag.posamp = (1, 2)
    assert tag.posamp == (1, 2)

    tag.minamp = 0.5
    tag.zamp = (0, 1)
    tag.xamp = (0, 1)
    tag.yamp = (0, 1)
    tag.num = 3
    tag.coord = (4, 6)
    tag.angamp = (0, 180)

    # Generic properties
    tag.rad = 1.5
    tag.occ = "1/3"
    tag.vol = 20.0
    tag.mult = 2
    tag.spin = 1.5

    # Tag properties
    tag.adatom = True
    assert tag.adatom is True
    tag.adatom = False
    assert tag.adatom is None

    tag.fix = True
    assert tag.fix is True
    tag.nomove = True
    assert tag.nomove is True
    tag.perm = True
    tag.athole = True


def test_serialization_snapshot_buildcell():
    """Capture exact serialization output for BuildcellParam"""
    bcp = BuildcellParam()
    bcp.nform = 3
    bcp.symmops = (2, 4)
    bcp.fix = True
    bcp.minsep = [2, {"Ce-O": (2, 3)}]
    bcp.formula = "Si2O4"
    bcp.sgrank = 2

    output = bcp.to_string()

    # Verify exact format
    assert "#NFORM=3" in output
    assert "#SYMMOPS=2-4" in output
    assert "#SGRANK=2" in output
    assert "#MINSEP=2 Ce-O=2-3" in output
    assert "#FORMULA=Si2O4" in output
    # FIX should NOT appear in output (special case)
    assert "#FIX" not in output


def test_serialization_snapshot_seedatomtag():
    """Capture exact serialization output for SeedAtomTag"""
    tag = SeedAtomTag()
    tag.tagname = "C1"
    tag.posamp = (1, 2)
    tag.num = 3
    tag.fix = True

    output = tag.to_string()

    # Verify exact format
    assert output.startswith("# C1 %")
    assert "POSAMP=1-2" in output
    assert "NUM=3" in output
    assert "FIX" in output


def test_numpy_array_storage():
    """Verify SeedAtomTag works correctly in numpy arrays"""
    import numpy as np

    # Create SeedAtomTag instances
    tags = [SeedAtomTag() for _ in range(5)]
    tags[0].tagname = "C0"
    tags[0].posamp = (1, 2)

    # Store in numpy array
    arr = np.array(tags, dtype=object)

    # Verify they're the same objects
    assert arr[0] is tags[0]

    # Verify modification works
    arr[1].tagname = "O1"
    assert tags[1].tagname == "O1"

    arr[2].num = 3
    assert tags[2].num == 3

    # Verify serialization
    assert "POSAMP=1-2" in arr[0].to_string()
    assert "O1" in arr[1].to_string()


def test_full_workflow_integration():
    """Integration test with SeedAtoms"""
    atoms = SeedAtoms("C4")

    # Set BuildcellParam properties
    atoms.gentags.symmops = (2, 3)
    atoms.gentags.sgrank = 2
    atoms.gentags.minsep = [2, {"C-C": (2, 3)}]
    atoms.gentags.compact = True

    # Set SeedAtomTag properties
    atoms[0].posamp = 3
    atoms[0].xamp = (2, 3)
    atoms[1].tagname = "CX0"
    atoms[2].fix = True

    # Get output
    lines = atoms.get_cell_inp_lines()
    buff = "\n".join(lines)

    # Verify all properties appear correctly
    assert "#SGRANK=2" in buff
    assert "#SYMMOPS=2-3" in buff
    assert "#MINSEP=2 C-C=2-3" in buff
    assert "XAMP=2-3" in buff
    assert "# C0" in buff
    assert "# CX0 " in buff


def test_property_deletion():
    """Test that properties can be deleted"""
    bcp = BuildcellParam()
    bcp.nform = 3
    assert bcp.nform == 3

    # Delete property
    del bcp.nform
    assert bcp.nform is None

    # Verify it doesn't appear in serialization
    output = bcp.to_string()
    assert "NFORM" not in output


def test_type_registry_tracking():
    """Test that type_registry correctly tracks property types"""
    bcp = BuildcellParam()

    # Set different property types
    bcp.fix = True  # tag
    bcp.formula = "Si2O4"  # generic
    bcp.nform = 3  # range
    bcp.minsep = (2, {})  # nested_range

    # Verify type registry
    assert bcp.type_registry["FIX"] == "tag"
    assert bcp.type_registry["FORMULA"] == "generic"
    assert bcp.type_registry["NFORM"] == "range"
    assert bcp.type_registry["MINSEP"] == "nested_range"


def test_clear_all():
    """Test clear_all method"""
    bcp = BuildcellParam()
    bcp.nform = 3
    bcp.formula = "Si2O4"
    bcp.fix = True

    # Clear all
    bcp.clear_all()

    # Verify all cleared
    assert bcp.nform is None
    assert bcp.formula is None
    assert bcp.fix is None
    assert len(bcp.prop_data) == 0

"""Tests for tools/modcell module."""

import pytest

from airsspy.tools.modcell import modify_cell, replace_block

SAMPLE_CELL = """\
%BLOCK LATTICE_CART
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Si  0.0 0.0 0.0
Si  0.5 0.5 0.5
%ENDBLOCK POSITIONS_FRAC

kpoints_mp_grid : 4 4 4
"""


def test_replace_block_basic():
    """Test replacing a block in cell file lines."""
    lines = SAMPLE_CELL.splitlines()
    new_lattice = ["2.0 0.0 0.0", "0.0 2.0 0.0", "0.0 0.0 2.0"]
    result = replace_block(lines, "LATTICE_CART", "LATTICE_CART", new_lattice)

    assert any("2.0 0.0 0.0" in line for line in result)
    assert "%BLOCK LATTICE_CART" in result
    assert "%ENDBLOCK LATTICE_CART" in result
    # Other content should remain
    assert any("kpoints_mp_grid" in line for line in result)


def test_replace_block_pattern_match():
    """Test that block pattern matching works for LATTICE_(CART|ABC)."""
    lines = SAMPLE_CELL.splitlines()
    new_lattice = ["3.0 0.0 0.0", "0.0 3.0 0.0", "0.0 0.0 3.0"]
    result = replace_block(lines, "LATTICE_CART", "LATTICE_(CART|ABC)", new_lattice)
    assert any("3.0 0.0 0.0" in line for line in result)


def test_replace_block_missing_block():
    """Test that missing block raises RuntimeError when check=True."""
    lines = ["no blocks here"]
    with pytest.raises(RuntimeError, match="Did not find start"):
        replace_block(lines, "NONEXISTENT", "NONEXISTENT", ["data"])


def test_replace_block_no_check():
    """Test that missing block is silently ignored when check=False."""
    lines = ["no blocks here"]
    result = replace_block(lines, "NONEXISTENT", "NONEXISTENT", ["data"], check=False)
    assert result == ["no blocks here"]


def test_modify_cell(tmp_path):
    """Test modifying a cell file with an ASE Atoms object."""
    from ase import Atoms

    # Write base cell file
    cell_file = tmp_path / "base.cell"
    cell_file.write_text(SAMPLE_CELL)

    atoms = Atoms(
        "Ge2",
        positions=[[0.0, 0.0, 0.0], [2.0, 2.0, 2.0]],
        cell=[[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]],
        pbc=True,
    )

    result = modify_cell(str(cell_file), atoms)

    # Check lattice was replaced
    result_str = "\n".join(result)
    assert "4.0000000000" in result_str
    # Check positions were replaced
    assert "Ge" in result_str
    # Check other content preserved
    assert any("kpoints_mp_grid" in line for line in result)

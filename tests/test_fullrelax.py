"""Tests for fullrelax module."""

import pytest

from airsspy.common import RelaxError
from airsspy.fullrelax import (
    FullRelax,
    check_relax_status,
    geom_to_cell,
    parse_geom_text_output,
)

GEOM_FILE_CONTENT = """\
begin header
  CASTEP
end header
        -1.95584852E+02                             <-- E
   5.41699994   0.00000000   0.00000000             <-- h
   0.00000000   5.41699994   0.00000000             <-- h
   0.00000000   0.00000000   5.41699994             <-- h
Si                1  0.0000000000  0.0000000000  0.0000000000    <-- R
Si                2  0.2500000000  0.2500000000  0.2500000000    <-- R

        -1.95600000E+02                             <-- E
   5.42000000   0.00000000   0.00000000             <-- h
   0.00000000   5.42000000   0.00000000             <-- h
   0.00000000   0.00000000   5.42000000             <-- h
Si                1  0.0000000000  0.0000000000  0.0000000000    <-- R
Si                2  0.2500000000  0.2500000000  0.2500000000    <-- R

"""

CASTEP_COMPLETED = """\
Geometry optimization completed
: finished iteration        1
: finished iteration        2
: finished iteration        3
"""

CASTEP_FAILED = """\
Geometry optimization failed
: finished iteration        1
: finished iteration        2
"""


def test_parse_geom_text_output():
    """Test parsing .geom file content."""
    lines = GEOM_FILE_CONTENT.splitlines()
    result = parse_geom_text_output(lines)

    assert "cells" in result
    assert "positions" in result
    assert "forces" in result
    assert "geom_energy" in result
    assert "symbols" in result

    assert result["symbols"] == ["Si", "Si"]
    assert result["cells"].shape[0] == 2  # Two geometry steps
    assert result["positions"].shape[0] == 2
    assert len(result["geom_energy"]) == 2


def test_parse_geom_text_output_energies():
    """Test energy values are converted from Hartree."""
    lines = GEOM_FILE_CONTENT.splitlines()
    result = parse_geom_text_output(lines)
    # First energy: -195.584852 Hartree -> ~-5321 eV
    assert result["geom_energy"][0] < 0
    assert abs(result["geom_energy"][0]) > 100  # Definitely in eV range


def test_parse_geom_empty():
    """Test that empty geom file raises RuntimeError."""
    with pytest.raises(RuntimeError, match="No data found"):
        parse_geom_text_output([])


def test_geom_to_cell(tmp_path):
    """Test converting last geom configuration to cell blocks."""
    geom_file = tmp_path / "test.geom"
    geom_file.write_text(GEOM_FILE_CONTENT)

    cell_block, pos_block = geom_to_cell(str(geom_file))

    assert "%BLOCK LATTICE_CART" in cell_block
    assert "%ENDBLOCK LATTICE_CART" in cell_block
    assert "%BLOCK POSITIONS_ABS" in pos_block
    assert "%ENDBLOCK POSITIONS_ABS" in pos_block
    assert "Si" in pos_block


def test_check_relax_status_completed(tmp_path):
    """Test relax status detection for completed optimisation."""
    castep_file = tmp_path / "test.castep"
    castep_file.write_text(CASTEP_COMPLETED)

    success, count = check_relax_status(str(castep_file))
    assert success is True
    assert count == 3


def test_check_relax_status_failed(tmp_path):
    """Test relax status detection for failed optimisation."""
    castep_file = tmp_path / "test.castep"
    castep_file.write_text(CASTEP_FAILED)

    success, count = check_relax_status(str(castep_file))
    assert success is False
    assert count == 2


def test_fullrelax_init():
    """Test FullRelax initialisation."""
    fr = FullRelax("castep.mpi", "test", maxit=200)
    assert fr.exe == "castep.mpi"
    assert fr.struct_name == "test"
    assert fr.maxit == 200
    assert fr._init_relax == 4  # default initial_cycle
    assert fr.success == 3  # default for non-zero maxit


def test_fullrelax_singlepoint():
    """Test FullRelax initialisation for single-point (maxit=0)."""
    fr = FullRelax("castep.mpi", "test", maxit=0)
    assert fr._init_relax == 0
    assert fr.success == 2


def test_fullrelax_properties():
    """Test FullRelax file name properties."""
    fr = FullRelax("castep.mpi", "Si-test", maxit=100)
    assert fr.dot_castep == "Si-test.castep"
    assert fr.dot_param == "Si-test.param"
    assert fr.dot_cell == "Si-test.cell"
    assert fr.dot_cell_out == "Si-test-out.cell"


def test_relax_error():
    """Test RelaxError exception."""
    assert issubclass(RelaxError, RuntimeError)
    with pytest.raises(RelaxError):
        raise RelaxError("test error")

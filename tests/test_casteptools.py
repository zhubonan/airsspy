"""Tests for casteptools module."""

import io

import pytest

from airsspy.casteptools import (
    CastepManualTimedout,
    CastepRunError,
    CastepSkip,
    castep_finish_ok,
    castep_geom_count,
    get_rand_cell_name,
    gulp_relax_finish_ok,
    parse_dot_castep,
    parse_param,
)

# --- Sample data ---

CASTEP_CONTENT = """\
 Everything done
         BFGS: finished iteration        1 with enthalpy=   -53.2326561 eV
               1      -53.23265609          0.00000000   0.0000000000        0.54    <-- SCF
         BFGS: finished iteration        2 with enthalpy=   -53.2400000 eV
               1      -53.24000000          0.00000000   0.0000000000        1.10    <-- SCF
 Writing model to .castep
         BFGS: finished iteration        3 with enthalpy=   -53.2450000 eV
               1      -53.24500000          0.00000000   0.0000000000        1.70    <-- SCF
 Total time  =     12.34 s
"""

CASTEP_NO_FINISH = """\
         BFGS: finished iteration        1 with enthalpy=   -53.2326561 eV
               1      -53.23265609          0.00000000   0.0000000000        0.54    <-- SCF
"""

PARAM_CONTENT = """\
# This is a comment
task : geometryoptimization
cut_off_energy : 500
xc_functional : PBE
"""

GULP_OUTPUT_WITH_FINAL = """\
  Final Enthalpy = -123.4567890     eV
"""

GULP_OUTPUT_NO_FINAL = """\
  Some output without final enthalpy
"""


# --- Tests ---


def test_parse_param(tmp_path):
    """Test parsing a .param file."""
    param_file = tmp_path / "test.param"
    param_file.write_text(PARAM_CONTENT)
    result = parse_param(str(tmp_path / "test"))

    assert result["task"] == "geometryoptimization"
    assert result["cut_off_energy"] == "500"
    assert result["xc_functional"] == "PBE"
    # Comments should be skipped
    assert "this" not in result


def test_parse_dot_castep_basic():
    """Test parsing geometry convergence from castep output."""
    fb = io.StringIO(CASTEP_CONTENT)
    result = parse_dot_castep(fb)

    assert result["H"] == [-53.2326561, -53.2400000, -53.2450000]
    assert result["iter_num"] == [1, 2, 3]
    assert result["name"] == "BFGS"
    assert result["unit"] == "eV"
    assert len(result["time"]) == 3


def test_parse_dot_castep_aggregate():
    """Test aggregate mode for timer correction."""
    fb = io.StringIO(CASTEP_CONTENT)
    result = parse_dot_castep(fb, aggregate=True)
    assert len(result["time"]) == 3
    # Times should be monotonically increasing
    times = list(result["time"])
    assert times == sorted(times)


def test_parse_dot_castep_save_times():
    """Test that save times are captured."""
    fb = io.StringIO(CASTEP_CONTENT)
    result = parse_dot_castep(fb)
    assert len(result["save_times"]) == 1
    assert len(result["save_iter"]) == 1


def test_castep_finish_ok(tmp_path):
    """Test checking for successful CASTEP completion."""
    castep_file = tmp_path / "test.castep"
    castep_file.write_text(CASTEP_CONTENT)
    assert castep_finish_ok(str(castep_file)) is True

    no_finish = tmp_path / "nofinish.castep"
    no_finish.write_text(CASTEP_NO_FINISH)
    assert castep_finish_ok(str(no_finish)) is False

    assert castep_finish_ok(str(tmp_path / "nonexistent.castep")) is False


def test_castep_geom_count(tmp_path):
    """Test counting geometry iterations."""
    castep_file = tmp_path / "test.castep"
    castep_file.write_text(CASTEP_CONTENT)
    # "finished iteration" appears 3 times, but "starting iteration" appears 0
    count = castep_geom_count(str(castep_file))
    assert count == 0  # pattern matches "starting iteration", not "finished iteration"


def test_castep_geom_count_with_starting(tmp_path):
    """Test counting with 'starting iteration' lines (lowercase match)."""
    lines = [
        "Starting geom iteration\n",  # Capital S - won't match
        "starting iteration\n",        # lowercase - matches
        "starting iteration\n",        # lowercase - matches
    ]
    castep_file = tmp_path / "test.castep"
    castep_file.write_text("".join(lines))
    assert castep_geom_count(str(castep_file)) == 2


def test_get_rand_cell_name():
    """Test random cell name generation."""
    name = get_rand_cell_name("Si")
    assert name.startswith("Si-")
    assert name.endswith(".cell")
    # Should have format Si-YYMMDD-HHMMSS-XXXXXX.cell
    parts = name.replace("Si-", "").replace(".cell", "").split("-")
    assert len(parts) == 3
    assert len(parts[2]) == 6  # uuid part


def test_gulp_relax_finish_ok(tmp_path):
    """Test GULP finish detection."""
    gulp_file = tmp_path / "test.gout"
    gulp_file.write_text(GULP_OUTPUT_WITH_FINAL)
    assert gulp_relax_finish_ok(str(gulp_file)) is True

    no_final = tmp_path / "nofinal.gout"
    no_final.write_text(GULP_OUTPUT_NO_FINAL)
    assert gulp_relax_finish_ok(str(no_final)) is False

    assert gulp_relax_finish_ok(str(tmp_path / "nonexistent.gout")) is False


def test_exception_classes():
    """Test exception classes exist and inherit correctly."""
    assert issubclass(CastepRunError, RuntimeError)
    assert issubclass(CastepSkip, RuntimeError)
    assert issubclass(CastepManualTimedout, RuntimeError)

    with pytest.raises(CastepRunError):
        raise CastepRunError("test error")

    with pytest.raises(CastepSkip):
        raise CastepSkip("skip")

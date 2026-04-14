"""Tests for gulptools module."""


from airsspy.gulptools import Ginfo, check_gulp, geom_opt_progress

GULP_NORMAL_OUTPUT = """\
  Cycle:     1 Energy: -123.4567890  Gnorm:    0.1234  CPU:  1.23
  Cycle:     2 Energy: -124.5678901  Gnorm:    0.0567  CPU:  2.46
  Cycle:     3 Energy: -125.6789012  Gnorm:    0.0012  CPU:  3.69
"""

GULP_OVERFLOW_OUTPUT = """\
  Cycle:     1 Energy: -123.4567890  Gnorm:    0.1234  CPU:  1.23
  Cycle:     2 Energy:************  Gnorm:*********  CPU:  2.46
"""

GULP_DIVERGING_OUTPUT = """\
  Cycle:     1 Energy: -123.4567890  Gnorm:    0.1234  CPU:  1.23
  Cycle:     2 Energy: -124.5678901  Gnorm:    5.0000  CPU:  2.46
  Cycle:     3 Energy: -124.1000000  Gnorm:   15.0000  CPU:  3.69
  Cycle:     4 Energy: -123.8000000  Gnorm:   20.0000  CPU:  4.92
  Cycle:     5 Energy: -123.5000000  Gnorm:   25.0000  CPU:  6.15
  Cycle:     6 Energy: -123.2000000  Gnorm:   30.0000  CPU:  7.38
  Cycle:     7 Energy: -122.9000000  Gnorm:   50.0000  CPU:  8.61
"""


def test_geom_opt_progress_from_lines():
    """Test parsing GULP geometry optimisation output."""
    data = geom_opt_progress(GULP_NORMAL_OUTPUT.splitlines())
    assert len(data) == 3
    assert isinstance(data[0], Ginfo)
    assert data[0].cycle == "1"
    assert data[0].energy == "-123.4567890"
    assert data[0].gnorm == "0.1234"
    assert data[0].cpu == "1.23"


def test_geom_opt_progress_from_file(tmp_path):
    """Test parsing from a file path."""
    path = tmp_path / "test.gout"
    path.write_text(GULP_NORMAL_OUTPUT)
    data = geom_opt_progress(str(path))
    assert len(data) == 3


def test_check_gulp_normal():
    """Test that a healthy GULP run passes the check."""
    assert check_gulp(GULP_NORMAL_OUTPUT.splitlines()) is True


def test_check_gulp_overflow():
    """Test that overflow markers cause check to fail."""
    assert check_gulp(GULP_OVERFLOW_OUTPUT.splitlines()) is False


def test_check_gulp_diverging():
    """Test that diverging Gnorm causes check to fail."""
    assert check_gulp(GULP_DIVERGING_OUTPUT.splitlines()) is False


def test_check_gulp_empty():
    """Test empty output passes (no bad patterns)."""
    assert check_gulp([]) is True
    assert check_gulp(["no geom data here"]) is True

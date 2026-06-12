"""Tests for scftools module."""


import pytest

from airsspy.scftools import SCFInfo

CASTEP_SCF_CONTENT = """\
 Welcome to CASTEP
 Starting electronic iteration
      1     -53.23265609          0.00000000   0.0000000000        0.54    <-- SCF
      2     -53.23500000          0.00000000  -0.0001000000        1.10    <-- SCF
      3     -53.24000000          0.00000000  -0.0000500000        1.70    <-- SCF
 Starting electronic iteration
      1     -53.24500000          0.00000000   0.0000000000        2.30    <-- SCF
      2     -53.24800000          0.00000000  -0.0000300000        2.90    <-- SCF
"""


@pytest.fixture
def castep_file(tmp_path):
    """Create a temporary .castep file with SCF data."""
    path = tmp_path / "test.castep"
    path.write_text(CASTEP_SCF_CONTENT)
    return str(path)


def test_scfinfo_parse(castep_file):
    """Test SCFInfo parses SCF cycles correctly."""
    scf = SCFInfo(castep_file)
    # Two geometry steps, each starting with loop=1
    assert len(scf) == 2


def test_scfinfo_data_fields(castep_file):
    """Test that parsed data has expected fields."""
    scf = SCFInfo(castep_file)
    cycle0 = scf.scf_data[0]
    assert len(cycle0.loops) == 3
    assert len(cycle0.energies) == 3
    assert len(cycle0.fermi_energies) == 3
    assert len(cycle0.gains) == 3
    assert len(cycle0.timers) == 3


def test_scfinfo_cycle_values(castep_file):
    """Test specific values from the first cycle."""
    scf = SCFInfo(castep_file)
    cycle0 = scf.scf_data[0]
    assert cycle0.loops == [1, 2, 3]
    assert cycle0.energies[0] == pytest.approx(-53.23265609)
    assert cycle0.timers[0] == pytest.approx(0.54)


def test_scfinfo_converge_data(castep_file):
    """Test convergence summary data."""
    scf = SCFInfo(castep_file)
    assert "energies" in scf.conv_data
    assert "timings" in scf.conv_data
    assert "durations" in scf.conv_data
    assert "steps" in scf.conv_data
    assert "average_scf_time" in scf.conv_data
    assert len(scf.conv_data["energies"]) == 2  # Two cycles


def test_scfinfo_get_summary(castep_file):
    """Test summary statistics."""
    scf = SCFInfo(castep_file)
    summary = scf.get_summary()
    assert "avg_ionic_time" in summary
    assert "avg_elec_time" in summary
    assert "avg_elec_steps" in summary
    assert "ionic_steps" in summary
    assert "total_time" in summary
    assert summary["ionic_steps"] == 2
    assert summary["avg_elec_steps"] == pytest.approx(2.5)  # (3 + 2) / 2


def test_scfinfo_reload(tmp_path):
    """Test that reload re-reads the file."""
    path = tmp_path / "test.castep"
    path.write_text(CASTEP_SCF_CONTENT)
    scf = SCFInfo(str(path))

    assert len(scf) == 2

    # Modify the file and reload
    path.write_text(CASTEP_SCF_CONTENT + "      1     -60.00000000          0.00000000   0.0000000000        3.50    <-- SCF\n")
    scf.reload()
    assert len(scf) == 3

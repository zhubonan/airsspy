"""Tests for ABACUS output parsing and result composition."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

# --- Sample log outputs ---

LOG_DEVELOP = """\
ITER   ETOT(eV)       EDIFF(eV)
CG1    -1.983062e+02  0.000000e+00
CG2    -1.969518e+02  1.354331e+00
#SCF IS CONVERGED#
#TOTAL ENERGY# -1.971286e+02 eV
TOTAL-PRESSURE: 20.928031 KBAR
Cell volume (Bohr^3) = 320.983
Cell volume (A^3) = 47.5648
STEP OF RELAXATION : 1
Largest force is 0.050000 eV/Angstrom
STEP OF RELAXATION : 2
Largest force is 0.001000 eV/Angstrom
Relaxation is converged!
"""

LOG_LTS = """\
ITER   ETOT(eV)       EDIFF(eV)
CG1    -1.983062e+02  0.000000e+00
charge density convergence is achieved
final etot is -1.971286e+02 eV
#TOTAL-PRESSURE# (EXCLUDE KINETIC PART OF IONS): 2.092803 GPa
Cell volume (A^3) = 47.5648
STEP OF RELAXATION : 1
STEP OF RELAXATION : 2
Relaxation is converged!
"""

LOG_FINAL_ETOT = """\
ITER   ETOT(eV)       EDIFF(eV)
CG1    -1.983062e+02  0.000000e+00
charge density convergence is achieved
!FINAL_ETOT_IS -197.128600000000 eV
#TOTAL-PRESSURE# (EXCLUDE KINETIC PART OF IONS): 0.000001 GPa
Cell volume (A^3) = 47.5648
STEP OF RELAXATION : 1
Relaxation is converged!
"""

LOG_BOHR_VOLUME = """\
!FINAL_ETOT_IS -197.128600000000 eV
TOTAL-PRESSURE: 0.000000 KBAR
Cell volume (Bohr^3) = 320.983
STEP OF RELAXATION : 1
"""

LOG_NOT_CONVERGED = """\
ITER   ETOT(eV)       EDIFF(eV)
CG1    -1.983062e+02  0.000000e+00
!!SCF IS NOT CONVERGED!!
#TOTAL ENERGY# -1.971286e+02 eV
TOTAL-PRESSURE: 0.000000 KBAR
STEP OF RELAXATION : 1
STEP OF RELAXATION : 2
STEP OF RELAXATION : 3
"""

LOG_EMPTY = ""

# --- Sample STRU content ---

STRU_DIRECT = """\
ATOMIC_SPECIES
Si 28.086 Si.UPF

LATTICE_CONSTANT
1.8897259886

LATTICE_VECTORS
5.430000 0.000000 0.000000
0.000000 5.430000 0.000000
0.000000 0.000000 5.430000

ATOMIC_POSITIONS
Direct

Si
0.0
2
0.000000 0.000000 0.000000 1 1 1
0.250000 0.250000 0.250000 1 1 1
"""

STRU_CARTESIAN = """\
ATOMIC_SPECIES
C 12.011 C.UPF

LATTICE_CONSTANT
1.8897259886

LATTICE_VECTORS
2.460000 0.000000 0.000000
0.000000 2.460000 0.000000
0.000000 0.000000 6.700000

ATOMIC_POSITIONS
Cartesian

C
0.0
4
0.000000 0.000000 1.675000 1 1 1
0.000000 0.000000 5.025000 1 1 1
0.615000 0.615000 0.000000 1 1 1
1.845000 1.845000 0.000000 1 1 1
"""


# --- parse_abacus_log tests ---


class TestParseAbacusLogEnergy:
    def test_final_etot_is(self, tmp_path, monkeypatch):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_FINAL_ETOT)
        result = parse_abacus_log(str(log))
        assert result["energy"] == pytest.approx(-197.1286)

    def test_total_energy_develop(self, tmp_path, monkeypatch):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_DEVELOP)
        result = parse_abacus_log(str(log))
        assert result["energy"] == pytest.approx(-197.1286)

    def test_final_etot_lts(self, tmp_path, monkeypatch):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_LTS)
        result = parse_abacus_log(str(log))
        assert result["energy"] == pytest.approx(-197.1286)

    def test_no_energy(self, tmp_path, monkeypatch):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_EMPTY)
        result = parse_abacus_log(str(log))
        assert result["energy"] is None

    def test_missing_file(self):
        from airsspy.abacustools import parse_abacus_log

        result = parse_abacus_log("/nonexistent/path.log")
        assert result["energy"] is None


class TestParseAbacusLogPressure:
    def test_pressure_kbar(self, tmp_path):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_DEVELOP)
        result = parse_abacus_log(str(log))
        assert result["pressure"] == pytest.approx(2.0928031)

    def test_pressure_gpa(self, tmp_path):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_LTS)
        result = parse_abacus_log(str(log))
        assert result["pressure"] == pytest.approx(2.092803)


class TestParseAbacusLogVolume:
    def test_volume_angstrom(self, tmp_path):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_DEVELOP)
        result = parse_abacus_log(str(log))
        assert result["volume"] == pytest.approx(47.5648)

    def test_volume_bohr(self, tmp_path):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_BOHR_VOLUME)
        result = parse_abacus_log(str(log))
        expected = 320.983 * 0.148184743
        assert result["volume"] == pytest.approx(expected, rel=1e-4)


class TestParseAbacusLogConvergence:
    def test_relaxation_converged(self, tmp_path):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_DEVELOP)
        result = parse_abacus_log(str(log))
        assert result["converged"] is True

    def test_not_converged(self, tmp_path):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_NOT_CONVERGED)
        result = parse_abacus_log(str(log))
        assert result["converged"] is False

    def test_scf_converged_develop(self, tmp_path):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_DEVELOP)
        result = parse_abacus_log(str(log))
        assert result["scf_converged"] is True

    def test_scf_converged_lts(self, tmp_path):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_LTS)
        result = parse_abacus_log(str(log))
        assert result["scf_converged"] is True

    def test_scf_not_converged(self, tmp_path):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_NOT_CONVERGED)
        result = parse_abacus_log(str(log))
        assert result["scf_converged"] is False

    def test_ionic_step_count(self, tmp_path):
        from airsspy.abacustools import parse_abacus_log

        log = tmp_path / "log"
        log.write_text(LOG_NOT_CONVERGED)
        result = parse_abacus_log(str(log))
        assert result["n_ionic_steps"] == 3


# --- parse_abacus_stru tests ---


class TestCellToStru:
    def test_positions_abs_writes_cartesian_positions(self):
        from airsspy.abacustools import cell_to_stru

        cell = """\
%BLOCK LATTICE_CART
5.0 0.0 0.0
0.0 5.0 0.0
0.0 0.0 5.0
%ENDBLOCK LATTICE_CART
%BLOCK POSITIONS_ABS
Si 2.5 0.0 0.0
%ENDBLOCK POSITIONS_ABS
%BLOCK SPECIES_POT
Si Si.UPF
%ENDBLOCK SPECIES_POT
"""

        stru = cell_to_stru(cell)

        assert "Si 28.086 Si.UPF" in stru
        assert "ATOMIC_POSITIONS\nCartesian" in stru
        assert "2.5000000000 0.0000000000 0.0000000000 1 1 1" in stru

    def test_positions_frac_takes_precedence_over_abs(self):
        from airsspy.abacustools import cell_to_stru

        cell = """\
%BLOCK LATTICE_CART
5.0 0.0 0.0
0.0 5.0 0.0
0.0 0.0 5.0
%ENDBLOCK LATTICE_CART
%BLOCK POSITIONS_ABS
Si 2.5 0.0 0.0
%ENDBLOCK POSITIONS_ABS
%BLOCK POSITIONS_FRAC
Si 0.25 0.0 0.0
%ENDBLOCK POSITIONS_FRAC
%BLOCK SPECIES_POT
Si Si.UPF
%ENDBLOCK SPECIES_POT
"""

        stru = cell_to_stru(cell)

        assert "ATOMIC_POSITIONS\nDirect" in stru
        assert "0.2500000000 0.0000000000 0.0000000000 1 1 1" in stru


class TestParseAbacusStru:
    def test_stru_direct(self, tmp_path):
        from airsspy.abacustools import parse_abacus_stru

        stru = tmp_path / "STRU"
        stru.write_text(STRU_DIRECT)
        elements, positions, cell = parse_abacus_stru(str(stru))

        assert elements == ["Si", "Si"]
        assert positions.shape == (2, 3)
        assert cell.shape == (3, 3)
        # Lattice constant ~1.8897 converts Bohr to Angstrom
        assert cell[0, 0] == pytest.approx(5.43, rel=1e-3)
        assert cell[1, 1] == pytest.approx(5.43, rel=1e-3)
        assert cell[2, 2] == pytest.approx(5.43, rel=1e-3)

    def test_stru_cartesian(self, tmp_path):
        from airsspy.abacustools import parse_abacus_stru

        stru = tmp_path / "STRU"
        stru.write_text(STRU_CARTESIAN)
        elements, positions, cell = parse_abacus_stru(str(stru))

        assert elements == ["C", "C", "C", "C"]
        assert positions.shape == (4, 3)
        assert cell.shape == (3, 3)
        # Positions should be fractional after conversion
        assert np.all(positions >= -0.01)
        assert np.all(positions <= 1.01)

    def test_stru_cartesian_skew_cell_matches_row_vector_conversion(self, tmp_path):
        from airsspy.abacustools import parse_abacus_stru

        stru = tmp_path / "STRU"
        stru.write_text(
            """\
ATOMIC_SPECIES
Ga 69.723 Ga.UPF

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
3.2459100000 0.0000000000 0.0000000000
-0.2621630475 4.5544108663 0.0000000000
-0.8895542149 -2.3255978880 5.9584787154

ATOMIC_POSITIONS
Cartesian

Ga
0.0
1
0.5392494814 -1.1326931555 3.5297265225 1 1 1
"""
        )

        _, positions, cell = parse_abacus_stru(str(stru))
        expected = np.array([[0.3328221, 0.0537855, 0.5923872]])
        assert np.allclose(positions, expected, atol=1e-7)
        assert np.allclose(positions @ cell, [[0.5392494814, -1.1326931555, 3.5297265225]])

    def test_stru_direct_positions(self, tmp_path):
        from airsspy.abacustools import parse_abacus_stru

        stru = tmp_path / "STRU"
        stru.write_text(STRU_DIRECT)
        elements, positions, cell = parse_abacus_stru(str(stru))

        # Direct positions should be preserved as-is
        assert positions[0, 0] == pytest.approx(0.0)
        assert positions[1, 0] == pytest.approx(0.25)


# --- compose_abacus_task_doc tests ---


@patch("airsspy.restools.save_airss_res")
@patch("airsspy.abacustools.parse_abacus_stru")
@patch("airsspy.abacustools.parse_abacus_log")
@patch("airsspy.abacustools.detect_logfile")
def test_compose_abacus_task_doc(
    mock_detect, mock_parse_log, mock_parse_stru, mock_save, tmp_path, monkeypatch
):
    from airsspy.abacustools import GPA_TO_EV_PER_ANG3, compose_abacus_task_doc

    monkeypatch.chdir(tmp_path)

    mock_detect.return_value = str(tmp_path / "test.abacus" / "OUT.ABACUS" / "running.log")
    mock_parse_log.return_value = {
        "energy": -197.1286,
        "pressure": 2.09,
        "volume": 47.5648,
        "converged": True,
        "scf_converged": True,
        "n_ionic_steps": 2,
    }
    mock_parse_stru.return_value = (
        ["Si", "Si"],
        np.array([[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]),
        np.eye(3) * 5.43,
    )

    # Create workdir with dummy files
    workdir = tmp_path / "test.abacus"
    out_dir = workdir / "OUT.ABACUS"
    out_dir.mkdir(parents=True)
    (out_dir / "STRU_ION_D").write_text("dummy")
    (workdir / "abacus_out").write_text("TOTAL  Time : 120.3")
    (tmp_path / "test.INPUT").write_text("calculation cell-relax\npress1 50\npress2 50\npress3 50\n")

    doc = compose_abacus_task_doc("test")

    assert doc["energy"] == pytest.approx(-197.1286)
    assert doc["pressure"] == pytest.approx(5.0)
    assert doc["natoms"] == 2
    assert doc["total_time"] == pytest.approx(120.3)
    info = mock_save.call_args.args[1]
    assert info["P"] == pytest.approx(5.0)
    assert info["H"] == pytest.approx(-197.1286 + 5.0 * 47.5648 * GPA_TO_EV_PER_ANG3)
    mock_save.assert_called_once()


@patch("airsspy.restools.save_airss_res")
@patch("airsspy.abacustools.parse_abacus_stru")
@patch("airsspy.abacustools.parse_abacus_log")
@patch("airsspy.abacustools.detect_logfile")
def test_compose_abacus_task_doc_falls_back_to_logged_pressure(
    mock_detect, mock_parse_log, mock_parse_stru, mock_save, tmp_path, monkeypatch
):
    from airsspy.abacustools import GPA_TO_EV_PER_ANG3, compose_abacus_task_doc

    monkeypatch.chdir(tmp_path)

    mock_detect.return_value = str(tmp_path / "test.abacus" / "OUT.ABACUS" / "running.log")
    mock_parse_log.return_value = {
        "energy": -197.1286,
        "pressure": 2.09,
        "volume": 47.5648,
        "converged": True,
        "scf_converged": True,
        "n_ionic_steps": 2,
    }
    mock_parse_stru.return_value = (
        ["Si", "Si"],
        np.array([[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]),
        np.eye(3) * 5.43,
    )

    workdir = tmp_path / "test.abacus"
    out_dir = workdir / "OUT.ABACUS"
    out_dir.mkdir(parents=True)
    (out_dir / "STRU_ION_D").write_text("dummy")

    doc = compose_abacus_task_doc("test")

    info = mock_save.call_args.args[1]
    assert doc["pressure"] == pytest.approx(2.09)
    assert info["P"] == pytest.approx(2.09)
    assert info["H"] == pytest.approx(-197.1286 + 2.09 * 47.5648 * GPA_TO_EV_PER_ANG3)


@patch("airsspy.abacustools.detect_logfile")
def test_compose_abacus_task_doc_requires_log(mock_detect, tmp_path, monkeypatch):
    from airsspy.abacustools import compose_abacus_task_doc

    monkeypatch.chdir(tmp_path)
    mock_detect.return_value = None

    with pytest.raises(RuntimeError, match="ABACUS log file not found"):
        compose_abacus_task_doc("test")

    assert not (tmp_path / "test.res").exists()


@patch("airsspy.abacustools.parse_abacus_log")
@patch("airsspy.abacustools.detect_logfile")
def test_compose_abacus_task_doc_rejects_unconverged_scf(
    mock_detect, mock_parse_log, tmp_path, monkeypatch
):
    from airsspy.abacustools import compose_abacus_task_doc

    monkeypatch.chdir(tmp_path)
    mock_detect.return_value = str(tmp_path / "test.abacus" / "OUT.ABACUS" / "running.log")
    mock_parse_log.return_value = {
        "energy": -197.1286,
        "pressure": 2.09,
        "volume": 47.5648,
        "converged": False,
        "scf_converged": False,
        "n_ionic_steps": 2,
    }

    with pytest.raises(RuntimeError, match="ABACUS SCF did not converge"):
        compose_abacus_task_doc("test")

    assert not (tmp_path / "test.res").exists()

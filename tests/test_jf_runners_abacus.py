"""Tests for AirssAbacusRelaxRunner."""

from pathlib import Path
from unittest.mock import patch

import pytest

from airsspy.jf.runners import AirssAbacusRelaxRunner

# Sample .cell content
CELL_CONTENT = """\
%BLOCK LATTICE_CART
  5.430000  0.000000  0.000000
  0.000000  5.430000  0.000000
  0.000000  0.000000  5.430000
%ENDBLOCK LATTICE_CART
%BLOCK POSITIONS_FRAC
Si   0.000000  0.000000  0.000000
Si   0.250000  0.250000  0.250000
%ENDBLOCK POSITIONS_FRAC
%BLOCK SPECIES_POT
Si Si_pbe_gga_6.0au_100Ry.upf
%ENDBLOCK SPECIES_POT
"""

# Sample INPUT content
INPUT_CONTENT = """\
calculation cell-relax
ecutwfc 50
pseudo_dir ../
"""


class TestPrepareInputs:
    def test_writes_cell_and_input(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)

        runner = AirssAbacusRelaxRunner()
        runner.prepare_inputs("Si-001", CELL_CONTENT, INPUT_CONTENT)

        assert Path("Si-001.cell").exists()
        assert Path("Si-001.INPUT").exists()
        assert Path("Si-001.abacus/STRU").exists()
        assert Path("Si-001.abacus/INPUT").exists()

        stru = Path("Si-001.abacus/STRU").read_text()
        assert "ATOMIC_SPECIES" in stru
        assert "Si" in stru
        assert "5.430000" in stru

    def test_cell_to_stru_invalid_lattice_raises(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        bad_cell = "%BLOCK POSITIONS_FRAC\nSi 0 0 0\n%ENDBLOCK POSITIONS_FRAC\n"

        runner = AirssAbacusRelaxRunner()
        with pytest.raises(ValueError, match="Expected 3 lattice vectors"):
            runner.prepare_inputs("Si-001", bad_cell, INPUT_CONTENT)

    def test_writes_external_pressure_in_kbar(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)

        runner = AirssAbacusRelaxRunner(pressure=5.0)
        runner.prepare_inputs("Si-001", CELL_CONTENT, INPUT_CONTENT)

        for input_file in ("Si-001.INPUT", "Si-001.abacus/INPUT"):
            content = Path(input_file).read_text()
            assert "press1 50.0" in content
            assert "press2 50.0" in content
            assert "press3 50.0" in content

    def test_singlepoint_axis_map_rotates_stru_and_kspacing(self, tmp_path, monkeypatch):
        from airsspy.jf.runners import AirssAbacusSinglePointRunner

        monkeypatch.chdir(tmp_path)
        input_content = INPUT_CONTENT + "kspacing 0.1 0.2 1.0\n"

        runner = AirssAbacusSinglePointRunner(cell_axis_map="z:x")
        runner.prepare_inputs("Si-001", CELL_CONTENT, input_content)

        stru = Path("Si-001.abacus/STRU").read_text()
        assert "5.4300000000  0.0000000000  0.0000000000" in stru
        assert "0.2500000000 0.2500000000 0.2500000000 1 1 1" in stru

        abacus_input = Path("Si-001.abacus/INPUT").read_text()
        assert "kspacing 1.0 0.1 0.2" in abacus_input


class TestSetInputParam:
    def test_set_existing_param(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        runner = AirssAbacusRelaxRunner()
        input_path = "test.INPUT"
        Path(input_path).write_text("calculation cell-relax\nrelax_nmax 50\n")
        runner._set_input_param(input_path, "relax_nmax", "3")
        content = Path(input_path).read_text()
        assert "relax_nmax 3" in content
        assert "relax_nmax 50" not in content

    def test_add_new_param(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        runner = AirssAbacusRelaxRunner()
        input_path = "test.INPUT"
        Path(input_path).write_text("calculation cell-relax\n")
        runner._set_input_param(input_path, "press1", "100.0")
        content = Path(input_path).read_text()
        assert "press1 100.0" in content


def _make_runner_result(converged: bool, n_steps: int = 1):
    """Create a result dict for mocked _run_single."""
    return {
        "energy": -197.1286,
        "pressure": 0.0,
        "volume": 160.18,
        "converged": converged,
        "scf_converged": True,
        "n_ionic_steps": n_steps,
    }


def _setup_workdir(tmp_path, struct_name="Si-001"):
    """Create the workdir and INPUT file that run() expects."""
    workdir = tmp_path / f"{struct_name}.abacus"
    workdir.mkdir()
    (workdir / "INPUT").write_text(INPUT_CONTENT)


class TestRunnerConverged:
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._update_cell")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._run_single")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner.prepare_inputs")
    def test_converged_in_phase1(
        self, mock_prepare, mock_run, mock_update, tmp_path, monkeypatch
    ):
        """All 3 rough runs converge — no phase 2 needed."""
        monkeypatch.chdir(tmp_path)
        _setup_workdir(tmp_path)
        mock_run.return_value = _make_runner_result(converged=True, n_steps=5)

        runner = AirssAbacusRelaxRunner(max_iterations=200)
        result = runner.run("Si-001", CELL_CONTENT, INPUT_CONTENT)

        assert result == 0
        assert mock_run.call_count == 3  # 3 rough runs

    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._update_cell")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._run_single")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner.prepare_inputs")
    def test_converged_in_phase2(
        self, mock_prepare, mock_run, mock_update, tmp_path, monkeypatch
    ):
        """Phase 1 not converged, phase 2 converges after 2 runs."""
        monkeypatch.chdir(tmp_path)
        _setup_workdir(tmp_path)
        # Phase 1: not converged; Phase 2: converged
        mock_run.side_effect = [
            _make_runner_result(converged=False, n_steps=3),  # rough 1
            _make_runner_result(converged=False, n_steps=3),  # rough 2
            _make_runner_result(converged=False, n_steps=3),  # rough 3
            _make_runner_result(converged=True, n_steps=10),  # phase 2 #1
            _make_runner_result(converged=True, n_steps=5),   # phase 2 #2
        ]

        runner = AirssAbacusRelaxRunner(max_iterations=200)
        result = runner.run("Si-001", CELL_CONTENT, INPUT_CONTENT)

        assert result == 0
        assert mock_run.call_count == 5  # 3 rough + 2 phase 2

    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._update_cell")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._run_single")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner.prepare_inputs")
    def test_not_converged_max_iterations(
        self, mock_prepare, mock_run, mock_update, tmp_path, monkeypatch
    ):
        """Never converges, stops at max_iterations."""
        monkeypatch.chdir(tmp_path)
        _setup_workdir(tmp_path)
        mock_run.return_value = _make_runner_result(converged=False, n_steps=50)

        runner = AirssAbacusRelaxRunner(max_iterations=200)
        result = runner.run("Si-001", CELL_CONTENT, INPUT_CONTENT)

        assert result == 1
        # Phase 1: 3 * 50 = 150 iter. Phase 2: 50 more = 200. Total 4 runs.
        assert mock_run.call_count == 4


class TestRunnerMaxIterations:
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._update_cell")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._run_single")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner.prepare_inputs")
    def test_max_iterations_exceeded_in_phase1(
        self, mock_prepare, mock_run, mock_update, tmp_path, monkeypatch
    ):
        """Max iterations exceeded during phase 1."""
        monkeypatch.chdir(tmp_path)
        _setup_workdir(tmp_path)
        mock_run.return_value = _make_runner_result(converged=False, n_steps=99)

        runner = AirssAbacusRelaxRunner(max_iterations=200)
        result = runner.run("Si-001", CELL_CONTENT, INPUT_CONTENT)

        assert result == 1
        # 99 * 2 = 198 < 200, 3rd run puts it at 297 > 200 → stop after 3
        assert mock_run.call_count == 3


class TestRunnerMaxFails:
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._run_single")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner.prepare_inputs")
    def test_max_fails(self, mock_prepare, mock_run, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        _setup_workdir(tmp_path)
        mock_run.return_value = None  # crash every time

        runner = AirssAbacusRelaxRunner(max_fails=2)
        result = runner.run("Si-001", CELL_CONTENT, INPUT_CONTENT)

        assert result == 1
        assert mock_run.call_count == 3  # initial + 2 retries

"""Tests for AirssAbacusRelaxRunner."""

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
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
"""

# Sample INPUT content
INPUT_CONTENT = """\
calculation cell-relax
ecutwfc 50
pseudo_dir ../
"""

# Sample STRU output from cell2stru
STRU_OUTPUT = """\
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


class TestPrepareInputs:
    @patch("airsspy.jf.runners.subprocess.run")
    def test_writes_cell_and_input(self, mock_run, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        mock_run.return_value = MagicMock(stdout=STRU_OUTPUT, returncode=0)

        runner = AirssAbacusRelaxRunner()
        runner.prepare_inputs("Si-001", CELL_CONTENT, INPUT_CONTENT)

        assert Path("Si-001.cell").exists()
        assert Path("Si-001.INPUT").exists()
        assert Path("Si-001.abacus/STRU").exists()
        assert Path("Si-001.abacus/INPUT").exists()

    @patch("airsspy.jf.runners.subprocess.run")
    def test_cell2stru_failure_raises(self, mock_run, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        mock_run.return_value = MagicMock(
            stdout="", stderr="cell2stru error", returncode=1
        )

        runner = AirssAbacusRelaxRunner()
        with pytest.raises(RuntimeError, match="cell2stru failed"):
            runner.prepare_inputs("Si-001", CELL_CONTENT, INPUT_CONTENT)


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
    def test_converged_after_two_successes(
        self, mock_prepare, mock_run, mock_update, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        _setup_workdir(tmp_path)
        mock_run.return_value = _make_runner_result(converged=True, n_steps=5)

        runner = AirssAbacusRelaxRunner(max_iterations=200)
        result = runner.run("Si-001", CELL_CONTENT, INPUT_CONTENT)

        assert result == 0
        # Phase 1 (3 rough runs) — all converge so success_counter >= 2
        assert mock_run.call_count == 3

    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._update_cell")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._run_single")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner.prepare_inputs")
    def test_not_converged_alternating(
        self, mock_prepare, mock_run, mock_update, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        _setup_workdir(tmp_path)
        mock_run.return_value = _make_runner_result(converged=False, n_steps=50)

        runner = AirssAbacusRelaxRunner(cycles=3, max_iterations=200)
        result = runner.run("Si-001", CELL_CONTENT, INPUT_CONTENT)

        assert result == 1


class TestRunnerMaxIterations:
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._update_cell")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._run_single")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner.prepare_inputs")
    def test_max_iterations_exceeded(
        self, mock_prepare, mock_run, mock_update, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        _setup_workdir(tmp_path)
        mock_run.return_value = _make_runner_result(converged=False, n_steps=99)

        runner = AirssAbacusRelaxRunner(cycles=10, max_iterations=200)
        result = runner.run("Si-001", CELL_CONTENT, INPUT_CONTENT)

        assert result == 1


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


class TestRunnerSinglePoint:
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._update_cell")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner._run_single")
    @patch("airsspy.jf.runners.AirssAbacusRelaxRunner.prepare_inputs")
    def test_single_point_negative_maxit(
        self, mock_prepare, mock_run, mock_update, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        _setup_workdir(tmp_path)
        mock_run.return_value = _make_runner_result(converged=True)

        runner = AirssAbacusRelaxRunner(max_iterations=-1)
        result = runner.run("Si-001", CELL_CONTENT, INPUT_CONTENT)

        assert result == 0
        # Should only run once for single-point
        assert mock_run.call_count == 1

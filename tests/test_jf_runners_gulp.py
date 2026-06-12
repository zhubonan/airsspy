"""Tests for GULP and pp3 runner classes."""

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from airsspy.jf.runners import AirssGulpRelaxRunner, AirssPp3RelaxRunner


class TestGulpRunnerCmdConstruction:
    """Test GULP runner command construction."""

    def test_default_cmd(self):
        runner = AirssGulpRelaxRunner()
        cmd = runner._get_cmd("test-001")
        assert cmd == ["gulp_relax", "ggulp", "0", "0.0", "test-001"]

    def test_custom_executable(self):
        runner = AirssGulpRelaxRunner(executable="gulp")
        cmd = runner._get_cmd("mystruct")
        assert cmd == ["gulp_relax", "gulp", "0", "0.0", "mystruct"]

    def test_cluster_and_pressure(self):
        runner = AirssGulpRelaxRunner(cluster=True, pressure=5.0)
        cmd = runner._get_cmd("test-001")
        assert cmd == ["gulp_relax", "ggulp", "1", "5.0", "test-001"]

    def test_param_suffix(self):
        runner = AirssGulpRelaxRunner()
        assert runner._param_suffix == ".lib"


class TestPp3RunnerCmdConstruction:
    """Test pp3 runner command construction."""

    def test_default_cmd(self):
        runner = AirssPp3RelaxRunner()
        cmd = runner._get_cmd("test-001")
        assert cmd == ["pp3_relax", "pp3", "test-001"]

    def test_custom_executable(self):
        runner = AirssPp3RelaxRunner(executable="pp3.opt")
        cmd = runner._get_cmd("mystruct")
        assert cmd == ["pp3_relax", "pp3.opt", "mystruct"]

    def test_param_suffix(self):
        runner = AirssPp3RelaxRunner()
        assert runner._param_suffix == ".pp"


class TestInputPreparation:
    """Test that runners write files with correct suffixes."""

    def test_gulp_writes_lib_file(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        runner = AirssGulpRelaxRunner()
        runner._prepare_inputs(
            "test-001",
            "%BLOCK LATTICE_CART\n%ENDBLOCK LATTICE_CART",
            "opt cellonly",
            seed_name="Si",
        )
        assert Path("test-001.cell").exists()
        assert Path("Si.lib").exists()
        assert not Path("test-001.lib").exists()

    def test_gulp_no_rename_when_no_seed(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        runner = AirssGulpRelaxRunner()
        runner._prepare_inputs(
            "test-001",
            "%BLOCK LATTICE_CART\n%ENDBLOCK LATTICE_CART",
            "opt cellonly",
        )
        assert Path("test-001.cell").exists()
        assert Path("test-001.lib").exists()

    def test_pp3_writes_pp_file(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        runner = AirssPp3RelaxRunner()
        runner._prepare_inputs(
            "test-001",
            "%BLOCK LATTICE_CART\n%ENDBLOCK LATTICE_CART",
            "pp3 input content",
        )
        assert Path("test-001.cell").exists()
        assert Path("test-001.pp").exists()


class TestSuccessDetection:
    """Test success checking logic."""

    def test_gulp_success_with_volume(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        runner = AirssGulpRelaxRunner()
        # Create a mock .castep file with "Final Enthalpy"
        Path("test-001.castep").write_text("Final Enthalpy = -123.456 eV\n")
        assert runner._check_success("test-001", "some output with Volume") is True

    def test_gulp_failure_no_volume(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        runner = AirssGulpRelaxRunner()
        Path("test-001.castep").write_text("Final Enthalpy = -123.456 eV\n")
        assert runner._check_success("test-001", "no volume info") is False

    def test_gulp_failure_no_final_enthalpy(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        runner = AirssGulpRelaxRunner()
        Path("test-001.castep").write_text("Incomplete run\n")
        assert runner._check_success("test-001", "some Volume") is False

    def test_pp3_success(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        runner = AirssPp3RelaxRunner()
        Path("test-001.castep").write_text("Final Enthalpy = -789.012 eV\n")
        assert runner._check_success("test-001", "") is True

    def test_pp3_failure(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        runner = AirssPp3RelaxRunner()
        Path("test-001.castep").write_text("Error in calculation\n")
        assert runner._check_success("test-001", "") is False


class TestRetryBehavior:
    """Test retry on timeout."""

    @patch("airsspy.jf.runners.subprocess.run")
    def test_retry_on_timeout(self, mock_run):
        mock_run.side_effect = subprocess.TimeoutExpired("cmd", 60)
        runner = AirssGulpRelaxRunner(max_attempts=2)
        result = runner.run(
            "test-001", "cell content", "lib content", seed_name="Si"
        )
        assert result == 1
        assert mock_run.call_count == 2

    @patch("airsspy.jf.runners.subprocess.run")
    @patch("airsspy.jf.runners.AirssGulpRelaxRunner._check_success")
    @patch("airsspy.jf.runners.AirssGulpRelaxRunner._prepare_inputs")
    def test_success_on_second_attempt(
        self, mock_prepare, mock_success, mock_run
    ):
        mock_run.side_effect = [
            subprocess.TimeoutExpired("cmd", 60),
            MagicMock(stdout="Volume output", returncode=0),
        ]
        mock_success.return_value = True
        runner = AirssGulpRelaxRunner(max_attempts=3)
        result = runner.run(
            "test-001", "cell content", "lib content", seed_name="Si"
        )
        assert result == 0
        assert mock_run.call_count == 2


class TestEndToEnd:
    """Test full run() with mocked subprocess."""

    @patch("airsspy.jf.runners.subprocess.run")
    @patch("airsspy.jf.runners.AirssGulpRelaxRunner._check_success")
    @patch("airsspy.jf.runners.AirssGulpRelaxRunner._prepare_inputs")
    def test_gulp_run_success(self, mock_prepare, mock_success, mock_run):
        mock_run.return_value = MagicMock(stdout="Volume output", returncode=0)
        mock_success.return_value = True
        runner = AirssGulpRelaxRunner()
        result = runner.run("test-001", "cell", "lib")
        assert result == 0

    @patch("airsspy.jf.runners.subprocess.run")
    @patch("airsspy.jf.runners.AirssGulpRelaxRunner._check_success")
    @patch("airsspy.jf.runners.AirssGulpRelaxRunner._prepare_inputs")
    def test_gulp_run_failure(self, mock_prepare, mock_success, mock_run):
        mock_run.return_value = MagicMock(stdout="no volume", returncode=1)
        mock_success.return_value = False
        runner = AirssGulpRelaxRunner(max_attempts=1)
        result = runner.run("test-001", "cell", "lib")
        assert result == 1

    @patch("airsspy.jf.runners.subprocess.run")
    @patch("airsspy.jf.runners.AirssPp3RelaxRunner._check_success")
    @patch("airsspy.jf.runners.AirssPp3RelaxRunner._prepare_inputs")
    def test_pp3_run_success(self, mock_prepare, mock_success, mock_run):
        mock_run.return_value = MagicMock(stdout="", returncode=0)
        mock_success.return_value = True
        runner = AirssPp3RelaxRunner()
        result = runner.run("test-001", "cell", "pp content")
        assert result == 0

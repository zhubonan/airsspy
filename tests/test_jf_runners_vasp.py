"""Tests for VASP runner classes."""

from pathlib import Path
from unittest.mock import MagicMock, patch

from airsspy.jf.runners import AirssVaspRelaxRunner, AirssVaspSinglePointRunner


@patch("airsspy.jf.runners.subprocess.run")
@patch("airsspy.vasptools.structure_from_cell_text")
@patch("airsspy.vasptools.prepare_vasp_inputs")
def test_vasp_relax_runner_prepares_and_runs(
    mock_prepare, mock_structure, mock_run, tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    mock_structure.return_value = object()
    mock_prepare.return_value = {"workdir": Path("Si-001.vasp")}
    mock_run.return_value = MagicMock(returncode=0)

    runner = AirssVaspRelaxRunner(executable="vasp_std", pressure=10.0)
    rc = runner.run("Si-001", "cell text", "ENCUT = 400")

    assert rc == 0
    assert Path("Si-001.cell").read_text() == "cell text"
    assert Path("Si-001.INCAR").read_text() == "ENCUT = 400"
    mock_prepare.assert_called_once()
    assert mock_prepare.call_args.kwargs["mode"] == "relax"
    assert mock_prepare.call_args.kwargs["pressure"] == 10.0
    mock_run.assert_called_once()
    assert mock_run.call_args.kwargs["cwd"] == Path("Si-001.vasp")


@patch("airsspy.vasptools.structure_from_cell_text")
@patch("airsspy.vasptools.prepare_vasp_inputs")
def test_vasp_single_point_uses_sp_mode(mock_prepare, mock_structure, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    mock_structure.return_value = object()
    mock_prepare.return_value = {"workdir": Path("Si-001.vasp")}

    runner = AirssVaspSinglePointRunner()
    runner.prepare_inputs("Si-001", "cell text", "ENCUT = 400")

    assert mock_prepare.call_args.kwargs["mode"] == "sp"

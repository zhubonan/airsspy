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

    def fake_run(*args, **kwargs):
        Path(kwargs["cwd"], "CONTCAR").write_text("relaxed structure")
        return MagicMock(returncode=0)

    mock_run.side_effect = fake_run
    runner = AirssVaspRelaxRunner(
        executable="vasp_std",
        pressure=10.0,
        max_iterations=4,
    )
    runner._read_vasp_status = MagicMock(side_effect=[(True, 1), (True, 1)])
    rc = runner.run("Si-001", "cell text", "ENCUT = 400")

    assert rc == 0
    assert Path("Si-001.cell").read_text() == "cell text"
    assert Path("Si-001.INCAR").read_text() == "ENCUT = 400"
    mock_prepare.assert_called_once()
    assert mock_prepare.call_args.kwargs["mode"] == "relax"
    assert mock_prepare.call_args.kwargs["pressure"] == 10.0
    assert mock_run.call_count == 2
    assert mock_run.call_args.kwargs["cwd"] == Path("Si-001.vasp")
    assert Path("Si-001.vasp/POSCAR").read_text() == "relaxed structure"


@patch("airsspy.jf.runners.subprocess.run")
@patch("airsspy.vasptools.structure_from_cell_text")
@patch("airsspy.vasptools.prepare_vasp_inputs")
def test_vasp_relax_runner_does_not_restart_after_unparsable_run(
    mock_prepare, mock_structure, mock_run, tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    workdir = Path("Si-001.vasp")
    workdir.mkdir()
    Path(workdir, "POSCAR").write_text("previous structure")
    Path(workdir, "CONTCAR").write_text("possibly stale structure")
    mock_structure.return_value = object()
    mock_prepare.return_value = {"workdir": workdir}
    mock_run.return_value = MagicMock(returncode=0)

    runner = AirssVaspRelaxRunner(max_fails=0, max_iterations=4)
    runner._read_vasp_status = MagicMock(return_value=None)
    rc = runner.run("Si-001", "cell text", "ENCUT = 400")

    assert rc == 1
    mock_run.assert_called_once()
    assert Path(workdir, "POSCAR").read_text() == "previous structure"


@patch("airsspy.jf.runners.subprocess.run")
@patch("airsspy.vasptools.structure_from_cell_text")
@patch("airsspy.vasptools.prepare_vasp_inputs")
def test_vasp_single_point_uses_sp_mode(
    mock_prepare, mock_structure, mock_run, tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    mock_structure.return_value = object()
    mock_prepare.return_value = {"workdir": Path("Si-001.vasp")}
    mock_run.return_value = MagicMock(returncode=0)

    runner = AirssVaspSinglePointRunner()
    rc = runner.run("Si-001", "cell text", "ENCUT = 400")

    assert rc == 0
    assert mock_prepare.call_args.kwargs["mode"] == "sp"
    mock_run.assert_called_once()

"""Tests for the native EDDPotentials.jl bridge runners."""

import json
import os
import shutil
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from ase import Atoms
from ase.io import read as ase_read
from click.testing import CliRunner

from airsspy.cli import cli, cmd_run
from airsspy.jf.eddp_runners import (
    AirssEddpRelaxRunner,
    AirssEddpSinglePointRunner,
    _bridge_path,
    compose_eddp_task_doc,
)

OUTPUT_RES = """\
TITL Si-001 0.000 250.000 -10.0000 0.00 0.00 2 (P1) n - 1
CELL 1.0 5.000000 5.000000 10.000000 90.000000 90.000000 90.000000
LATT -1
SFAC Si
Si 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0
Si 1 0.5000000000000 0.5000000000000 0.5000000000000 1.0
END
"""


def _fake_julia_result(command, **kwargs):
    Path(command[-9]).write_text(OUTPUT_RES)
    Path(command[-8]).write_text(
        json.dumps(
            {
                "energy": -10.5,
                # Julia bridge writes the 3 x nat EDDP force matrix by rows.
                "forces": [[0.1, -0.1], [0.2, -0.2], [0.3, -0.3]],
                "stress": [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]],
                "pressure_gpa": -0.320435,
                "converged": True,
                "iterations": 7,
                "fmax": 0.03,
                "smax_gpa": 0.08,
            }
        )
    )
    return SimpleNamespace(returncode=0, stdout="native EDDP ok", stderr="")


def test_native_relax_runner_invokes_bridge_and_writes_extxyz(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    model = tmp_path / "model.json"
    model.write_text("{}")
    project = tmp_path / "EDDPotentials.jl"
    project.mkdir()
    atoms = Atoms(
        "Si2",
        scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        cell=[5, 5, 5],
        pbc=True,
    )

    runner = AirssEddpRelaxRunner(
        model,
        executable="julia --startup-file=no",
        project=project,
        method="fire",
        max_steps=25,
        force_tolerance=0.04,
        stress_tolerance_gpa=0.2,
        pressure=3.0,
        relax_cell=False,
    )
    with patch(
        "airsspy.jf.eddp_runners.subprocess.run", side_effect=_fake_julia_result
    ):
        assert runner.run("Si-001", atoms) == 0

    assert runner.last_command[:2] == ["julia", "--startup-file=no"]
    assert f"--project={project}" in runner.last_command
    assert runner.last_command[-12:] == [
        "relax",
        str(model),
        "Si-001.eddp-input.res",
        "Si-001.eddp-output.res",
        "Si-001.eddp-result.json",
        "Si-001",
        "fire",
        "25",
        "0.04",
        "0.2",
        "3.0",
        "false",
    ]
    result = ase_read("Si-001.extxyz")
    assert result.get_potential_energy() == pytest.approx(-10.5)
    assert result.get_forces() == pytest.approx(
        np.array([[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3]])
    )
    assert result.info["relax_converged"]
    assert result.info["relax_steps"] == 7
    assert "native EDDP ok" in Path("Si-001.eddp.log").read_text()

    doc = compose_eddp_task_doc("Si-001", str(model))
    assert f"EDDP Model {model}" in doc["res_content"]
    assert "EDDP Relax status converged" in doc["res_content"]
    assert "ML Relax" not in doc["res_content"]


def test_native_singlepoint_runner_uses_singlepoint_mode(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    runner = AirssEddpSinglePointRunner(tmp_path / "model.json")
    atoms = Atoms("Si2", positions=[[0, 0, 0], [1, 1, 1]], cell=[5, 5, 5], pbc=True)

    with patch(
        "airsspy.jf.eddp_runners.subprocess.run", side_effect=_fake_julia_result
    ):
        assert runner.run("Si-001", atoms) == 0

    assert runner.last_command[-12] == "singlepoint"
    result = ase_read("Si-001.extxyz")
    assert "relax_converged" not in result.info
    assert "relax_status" not in result.info


def test_native_runner_reports_failed_julia_process(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    runner = AirssEddpSinglePointRunner(tmp_path / "model.json")
    failed = SimpleNamespace(returncode=2, stdout="", stderr="bad artifact")

    Path("Si-001.extxyz").write_text("stale result")
    with patch("airsspy.jf.eddp_runners.subprocess.run", return_value=failed):
        assert runner.run("Si-001", Atoms("Si", cell=[5, 5, 5], pbc=True)) == 1

    assert not Path("Si-001.extxyz").exists()
    assert "bad artifact" in Path("Si-001.eddp.log").read_text()


@pytest.mark.parametrize("method", ["lbfgs", "bfgs"])
def test_native_relax_runner_rejects_unknown_method(method):
    with pytest.raises(ValueError, match="tpsd.*fire"):
        AirssEddpRelaxRunner("model.json", method=method)


def test_cli_factories_create_native_eddp_runners(tmp_path):
    model = tmp_path / "model.json"
    model.write_text("{}")
    project = tmp_path / "EDDPotentials.jl"
    project.mkdir()

    relax = cmd_run._create_runner(
        "eddp",
        "julia",
        44,
        False,
        2.5,
        calculator_spec=str(model),
        fmax=0.03,
        eddp_project=str(project),
        eddp_method="fire",
        eddp_stress_tol=0.15,
        eddp_fixed_cell=True,
    )
    singlepoint = cmd_run._create_sp_runner(
        "eddp",
        "julia",
        calculator_spec=str(model),
        eddp_project=str(project),
        pressure=2.5,
    )

    assert isinstance(relax, AirssEddpRelaxRunner)
    assert relax.method == "fire"
    assert relax.max_steps == 44
    assert relax.force_tolerance == pytest.approx(0.03)
    assert relax.stress_tolerance_gpa == pytest.approx(0.15)
    assert relax.relax_cell is False
    assert isinstance(singlepoint, AirssEddpSinglePointRunner)
    assert singlepoint.pressure == pytest.approx(2.5)


def test_packaged_bridge_uses_native_eddp_apis():
    bridge = _bridge_path()
    content = bridge.read_text()

    assert bridge.is_file()
    for api in (
        "load_calculator",
        "eddp_energy",
        "eddp_forces",
        "eddp_stress",
        "multirelax!",
    ):
        assert api in content


@pytest.mark.e2e
def test_real_julia_eddp_singlepoint_when_configured(tmp_path, monkeypatch):
    """Exercise the real Julia API when CI supplies an EDDP model artifact."""
    julia = shutil.which("julia")
    model = os.environ.get("AIRSSPY_TEST_EDDP_MODEL")
    project = os.environ.get("AIRSSPY_TEST_EDDP_PROJECT")
    input_res = os.environ.get("AIRSSPY_TEST_EDDP_RES")
    if julia is None or not model or not input_res:
        pytest.skip(
            "requires julia, AIRSSPY_TEST_EDDP_MODEL, and AIRSSPY_TEST_EDDP_RES"
        )

    from airsspy.restools import read_res_atoms

    _, atoms = read_res_atoms(Path(input_res).read_text().splitlines())
    monkeypatch.chdir(tmp_path)
    runner = AirssEddpSinglePointRunner(
        model,
        executable=julia,
        project=project,
    )
    assert runner.run("Si-001", atoms) == 0
    result = ase_read("Si-001.extxyz")
    assert np.isfinite(result.get_potential_energy())
    assert result.get_forces().shape == (len(atoms), 3)

    compose_eddp_task_doc("Si-001", model)
    from airsspy.convert import _parse_res_forces
    from airsspy.restools import read_res_atoms

    res_lines = Path("Si-001.res").read_text().splitlines()
    _, res_atoms = read_res_atoms(res_lines)
    res_forces = _parse_res_forces(res_lines)
    assert len(res_atoms) == len(atoms)
    assert res_forces is not None
    assert len(res_forces) == len(atoms)


@pytest.mark.e2e
def test_real_julia_eddp_relax_when_configured(tmp_path, monkeypatch):
    """Exercise native variable-cell relaxation on an already relaxed fixture."""
    julia = shutil.which("julia")
    model = os.environ.get("AIRSSPY_TEST_EDDP_MODEL")
    project = os.environ.get("AIRSSPY_TEST_EDDP_PROJECT")
    input_res = os.environ.get("AIRSSPY_TEST_EDDP_RELAXED_RES")
    if julia is None or not model or not input_res:
        pytest.skip(
            "requires julia, AIRSSPY_TEST_EDDP_MODEL, and AIRSSPY_TEST_EDDP_RELAXED_RES"
        )

    from airsspy.restools import read_res_atoms

    _, atoms = read_res_atoms(Path(input_res).read_text().splitlines())
    monkeypatch.chdir(tmp_path)
    runner = AirssEddpRelaxRunner(
        model,
        executable=julia,
        project=project,
        method="fire",
        max_steps=10,
        force_tolerance=0.01,
        stress_tolerance_gpa=0.1,
        pressure=0.01,
    )
    assert runner.run("LiTaOCl4-relax", atoms) == 0
    assert runner.last_result is not None
    assert runner.last_result["converged"] is True
    assert runner.last_result["fmax"] <= 0.01
    assert runner.last_result["smax_gpa"] <= 0.1


def test_run_search_routes_native_eddp_options():
    click_runner = CliRunner()
    with click_runner.isolated_filesystem():
        model = Path("model.json")
        model.write_text("{}")
        project = Path("EDDPotentials.jl")
        project.mkdir()
        Path("Si.cell").write_text("#SPECIES=Si\n#NATOM=1\n")
        Path("Si-001.cell").write_text("cell content\n")
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        build_result = {
            "struct_name": "Si-001",
            "seed_name": "Si",
            "struct_content": "cell content\n",
        }
        with patch("airsspy.jf.runners.run_buildcell", return_value=build_result):
            with patch(
                "airsspy.cli.cmd_run._create_runner", return_value=fake_runner
            ) as create:
                with patch("airsspy.cli.cmd_run._collect_result") as collect:
                    result = click_runner.invoke(
                        cli,
                        [
                            "run",
                            "search",
                            "--seed",
                            "Si",
                            "--nmax",
                            "1",
                            "--code",
                            "eddp",
                            "--calculator",
                            str(model),
                            "--eddp-project",
                            str(project),
                            "--eddp-method",
                            "fire",
                            "--eddp-fixed-cell",
                        ],
                    )

    assert result.exit_code == 0, result.output
    assert create.call_args.kwargs["eddp_method"] == "fire"
    assert create.call_args.kwargs["eddp_fixed_cell"] is True
    fake_runner.run.assert_called_once_with("Si-001", "cell content\n")
    collect.assert_called_once()


@pytest.mark.parametrize("command", ["relax", "sp"])
def test_run_existing_res_routes_native_eddp(command):
    click_runner = CliRunner()
    with click_runner.isolated_filesystem():
        Path("model.json").write_text("{}")
        Path("Si-001.res").write_text(OUTPUT_RES)
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        factory = (
            "airsspy.cli.cmd_run._create_task_runner"
            if command == "relax"
            else "airsspy.cli.cmd_run._create_sp_runner"
        )
        with patch(factory, return_value=fake_runner) as create:
            with patch("airsspy.cli.cmd_run._collect_result") as collect:
                result = click_runner.invoke(
                    cli,
                    [
                        "run",
                        command,
                        "--cell",
                        "*.res",
                        "--code",
                        "eddp",
                        "--calculator",
                        "model.json",
                    ],
                )

    assert result.exit_code == 0, result.output
    assert create.call_args.kwargs["calculator_spec"].endswith("model.json")
    fake_runner.run.assert_called_once()
    assert isinstance(fake_runner.run.call_args.args[1], Atoms)
    collect.assert_called_once()


def test_run_crud_routes_native_eddp_without_root_input_files():
    click_runner = CliRunner()
    with click_runner.isolated_filesystem():
        Path("model.json").write_text("{}")
        Path("hopper").mkdir()
        Path("hopper/Si-001.res").write_text(OUTPUT_RES)
        fake_runner = MagicMock()
        fake_runner.run.return_value = 0
        with patch(
            "airsspy.cli.cmd_run._create_task_runner", return_value=fake_runner
        ) as create:
            with patch("airsspy.cli.cmd_run._collect_result") as collect:
                result = click_runner.invoke(
                    cli,
                    [
                        "run",
                        "crud",
                        "--code",
                        "eddp",
                        "--calculator",
                        "model.json",
                        "--keep",
                    ],
                )

    assert result.exit_code == 0, result.output
    assert create.call_args.kwargs["singlepoint"] is False
    fake_runner.run.assert_called_once()
    assert isinstance(fake_runner.run.call_args.args[1], Atoms)
    collect.assert_called_once()


def test_crud_cleanup_removes_eddp_intermediate_res_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path("Si-001.res").write_text(OUTPUT_RES)
    Path("Si-001.eddp-input.res").write_text(OUTPUT_RES)
    Path("Si-001.eddp-output.res").write_text(OUTPUT_RES)
    Path("Si-001.eddp-result.json").write_text("{}")

    cmd_run._cleanup_crud_artifacts("Si-001")

    assert Path("Si-001.res").is_file()
    assert not Path("Si-001.eddp-input.res").exists()
    assert not Path("Si-001.eddp-output.res").exists()
    assert not Path("Si-001.eddp-result.json").exists()

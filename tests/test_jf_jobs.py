"""Tests for jobflow Makers (AirssSearchMaker, AirssRelaxMaker, AirssValidateMaker)."""

import subprocess
from unittest.mock import MagicMock, patch

import pytest

jobflow = pytest.importorskip("jobflow")

from jobflow import run_locally  # noqa: E402
from castepinput.inputs import CellInput, ParamInput  # noqa: E402
from pymatgen.core import Lattice, Structure  # noqa: E402

from airsspy.jf.documents import RelaxOutcome  # noqa: E402
from airsspy.jf.jobs import (  # noqa: E402
    AirssRelaxMaker,
    AirssSearchMaker,
    AirssValidateMaker,
)


def _make_task_doc(**overrides):
    """Create a default compose_task_doc return dict with optional overrides."""
    defaults = {
        "structure": None,
        "volume": 100.0,
        "reduced_formula": "Si",
        "formula": "Si4",
        "natoms": 4,
        "energy": -10.0,
        "energy_per_atom": -2.5,
        "pressure": 0.1,
        "spin": 0.0,
        "mod_spin": 0.0,
        "symmetry": "(Fd-3m)",
        "res_content": "TITL Si\nCELL ...",
        "parallel_efficiency": 0.8,
        "total_time": 100.0,
    }
    defaults.update(overrides)
    return defaults


def _make_si_structure():
    """Create a minimal Si structure for testing."""
    return Structure(Lattice.cubic(5.0), ["Si"], [[0, 0, 0]])


def _make_paraminput():
    """Create a real ParamInput for testing."""
    p = ParamInput()
    p["task"] = "geometryoptimization"
    return p


def _make_cellinput():
    """Create a real CellInput for testing."""
    return CellInput()


# --- AirssSearchMaker tests ---


def test_search_maker_castep_success():
    maker = AirssSearchMaker(n_structures=1, code="castep")
    job = maker.make(
        seed_name="Si",
        seed_content="seed cell content",
        paraminput=_make_paraminput(),
        project_name="test_project",
    )

    with (
        patch("airsspy.jf.jobs.run_buildcell") as mock_buildcell,
        patch("airsspy.jf.jobs.AirssCastepRelaxRunner") as mock_runner_cls,
        patch("airsspy.jf.jobs.compose_task_doc") as mock_compose,
        patch("castepinput.inputs.CellInput") as mock_cellinput_cls,
    ):
        mock_buildcell.return_value = {
            "struct_name": "Si-001",
            "seed_name": "Si",
            "struct_content": "%BLOCK LATTICE_CART\n1.0 0 0\n0 1.0 0\n0 0 1.0\n%ENDBLOCK LATTICE_CART",
            "seed_hash": "abc123",
        }
        mock_runner = MagicMock()
        mock_runner.run.return_value = 0
        mock_runner_cls.return_value = mock_runner
        mock_compose.return_value = _make_task_doc()
        mock_cellinput_cls.from_file.return_value = MagicMock()

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.project_name == "test_project"
    assert output.seed_name == "Si"
    assert output.n_structures == 1
    assert output.n_finished == 1
    assert output.n_errored == 0
    assert len(output.results) == 1
    assert output.results[0].relax_status == RelaxOutcome.FINISHED


def test_search_maker_buildcell_timeout():
    maker = AirssSearchMaker(n_structures=1)
    job = maker.make(
        seed_name="Si",
        seed_content="seed cell content",
        paraminput=_make_paraminput(),
        project_name="test",
    )

    with patch("airsspy.jf.jobs.run_buildcell") as mock_buildcell:
        mock_buildcell.return_value = None
        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_structures == 0
    assert len(output.results) == 0


def test_search_maker_relax_error():
    maker = AirssSearchMaker(n_structures=1)
    job = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=_make_paraminput(),
        project_name="test",
    )

    with (
        patch("airsspy.jf.jobs.run_buildcell") as mock_buildcell,
        patch("airsspy.jf.jobs.AirssCastepRelaxRunner") as mock_runner_cls,
        patch("airsspy.jf.jobs.compose_task_doc") as mock_compose,
        patch("castepinput.inputs.CellInput") as mock_cellinput_cls,
    ):
        mock_buildcell.return_value = {
            "struct_name": "Si-001",
            "seed_name": "Si",
            "struct_content": "cell",
        }
        mock_runner = MagicMock()
        mock_runner.run.return_value = 1
        mock_runner_cls.return_value = mock_runner
        mock_compose.return_value = _make_task_doc()
        mock_cellinput_cls.from_file.return_value = MagicMock()

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_errored == 1
    assert output.n_finished == 0
    assert output.results[0].relax_status == RelaxOutcome.ERRORED


def test_search_maker_gulp_success():
    maker = AirssSearchMaker(n_structures=1, code="gulp")
    job = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=_make_paraminput(),
        project_name="test",
    )

    with (
        patch("airsspy.jf.jobs.run_buildcell") as mock_buildcell,
        patch("airsspy.jf.jobs.AirssGulpRelaxRunner") as mock_runner_cls,
        patch("airsspy.jf.jobs.compose_task_doc") as mock_compose,
        patch("pathlib.Path.read_text", return_value="cell content"),
    ):
        mock_buildcell.return_value = {
            "struct_name": "Si-001",
            "seed_name": "Si",
            "struct_content": "cell",
        }
        mock_runner = MagicMock()
        mock_runner.run.return_value = 0
        mock_runner_cls.return_value = mock_runner
        mock_compose.return_value = _make_task_doc()

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_finished == 1


def test_search_maker_pp3_success():
    maker = AirssSearchMaker(n_structures=1, code="pp3")
    job = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=_make_paraminput(),
        project_name="test",
    )

    with (
        patch("airsspy.jf.jobs.run_buildcell") as mock_buildcell,
        patch("airsspy.jf.jobs.AirssPp3RelaxRunner") as mock_runner_cls,
        patch("airsspy.jf.jobs.compose_task_doc") as mock_compose,
        patch("pathlib.Path.read_text", return_value="cell content"),
    ):
        mock_buildcell.return_value = {
            "struct_name": "Si-001",
            "seed_name": "Si",
            "struct_content": "cell",
        }
        mock_runner = MagicMock()
        mock_runner.run.return_value = 0
        mock_runner_cls.return_value = mock_runner
        mock_compose.return_value = _make_task_doc()

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_finished == 1


def test_search_maker_abacus_success():
    maker = AirssSearchMaker(n_structures=1, code="abacus")
    job = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=_make_paraminput(),
        project_name="test",
    )

    with (
        patch("airsspy.jf.jobs.run_buildcell") as mock_buildcell,
        patch("airsspy.jf.jobs.AirssAbacusRelaxRunner") as mock_runner_cls,
        patch("airsspy.abacustools.compose_abacus_task_doc") as mock_compose,
        patch("pathlib.Path.read_text", return_value="cell content"),
    ):
        mock_buildcell.return_value = {
            "struct_name": "Si-001",
            "seed_name": "Si",
            "struct_content": "cell",
        }
        mock_runner = MagicMock()
        mock_runner.run.return_value = 0
        mock_runner_cls.return_value = mock_runner
        mock_compose.return_value = _make_task_doc()

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_finished == 1


def test_search_maker_invalid_code():
    maker = AirssSearchMaker(n_structures=1, code="vasp")
    job = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=_make_paraminput(),
        project_name="test",
    )

    with (
        patch("airsspy.jf.jobs.run_buildcell") as mock_buildcell,
    ):
        mock_buildcell.return_value = {
            "struct_name": "Si-001",
            "seed_name": "Si",
            "struct_content": "cell",
        }

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_failed == 1
    assert output.results[0].relax_status == RelaxOutcome.FAILED
    assert "Unknown code" in output.results[0].error_message


def test_search_maker_stop_if_all_errored():
    maker = AirssSearchMaker(n_structures=1, stop_if_not_converged=True)
    job = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=_make_paraminput(),
        project_name="test",
    )

    with (
        patch("airsspy.jf.jobs.run_buildcell") as mock_buildcell,
        patch("airsspy.jf.jobs.AirssCastepRelaxRunner") as mock_runner_cls,
        patch("airsspy.jf.jobs.compose_task_doc") as mock_compose,
        patch("castepinput.inputs.CellInput") as mock_cellinput_cls,
    ):
        mock_buildcell.return_value = {
            "struct_name": "Si-001",
            "seed_name": "Si",
            "struct_content": "cell",
        }
        mock_runner = MagicMock()
        mock_runner.run.return_value = 1
        mock_runner_cls.return_value = mock_runner
        mock_compose.return_value = _make_task_doc()
        mock_cellinput_cls.from_file.return_value = MagicMock()

        responses = run_locally(job, ensure_success=True)

    response = responses[job.uuid][1]
    assert response.stop_children is True


def test_search_maker_multiple_structures():
    """Test that multiple structures are iterated correctly."""
    maker = AirssSearchMaker(n_structures=3)
    job = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=_make_paraminput(),
        project_name="test",
    )

    with (
        patch("airsspy.jf.jobs.run_buildcell") as mock_buildcell,
        patch("airsspy.jf.jobs.AirssCastepRelaxRunner") as mock_runner_cls,
        patch("airsspy.jf.jobs.compose_task_doc") as mock_compose,
        patch("castepinput.inputs.CellInput") as mock_cellinput_cls,
    ):
        mock_buildcell.return_value = {
            "struct_name": "Si-001",
            "seed_name": "Si",
            "struct_content": "cell",
        }
        mock_runner = MagicMock()
        mock_runner.run.return_value = 0
        mock_runner_cls.return_value = mock_runner
        mock_compose.return_value = _make_task_doc()
        mock_cellinput_cls.from_file.return_value = MagicMock()

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_structures == 3
    assert mock_buildcell.call_count == 3


# --- AirssRelaxMaker tests ---


def test_relax_maker_castep_success():
    maker = AirssRelaxMaker(code="castep")
    job = maker.make(
        structures=[_make_si_structure()],
        struct_names=["Si-001"],
        cellinputs=[_make_cellinput()],
        paraminput=_make_paraminput(),
        project_name="test",
        seed_name="Si",
    )

    with (
        patch("airsspy.jf.jobs.AirssCastepRelaxRunner") as mock_runner_cls,
        patch("airsspy.jf.jobs.compose_task_doc") as mock_compose,
    ):
        mock_runner = MagicMock()
        mock_runner.run.return_value = 0
        mock_runner_cls.return_value = mock_runner
        mock_compose.return_value = _make_task_doc()

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_finished == 1
    assert output.job_type == "relax"


def test_relax_maker_gulp_success():
    maker = AirssRelaxMaker(code="gulp", executable="ggulp")
    job = maker.make(
        structures=[_make_si_structure()],
        struct_names=["Si-001"],
        cellinputs=[_make_cellinput()],
        paraminput=_make_paraminput(),
        project_name="test",
        seed_name="Si",
    )

    with (
        patch("airsspy.jf.jobs.AirssGulpRelaxRunner") as mock_runner_cls,
        patch("airsspy.jf.jobs.compose_task_doc") as mock_compose,
    ):
        mock_runner = MagicMock()
        mock_runner.run.return_value = 0
        mock_runner_cls.return_value = mock_runner
        mock_compose.return_value = _make_task_doc()

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_finished == 1


def test_relax_maker_abacus_success():
    maker = AirssRelaxMaker(code="abacus", executable="abacus")
    job = maker.make(
        structures=[_make_si_structure()],
        struct_names=["Si-001"],
        cellinputs=[_make_cellinput()],
        paraminput=_make_paraminput(),
        project_name="test",
        seed_name="Si",
    )

    with (
        patch("airsspy.jf.jobs.AirssAbacusRelaxRunner") as mock_runner_cls,
        patch("airsspy.abacustools.compose_abacus_task_doc") as mock_compose,
    ):
        mock_runner = MagicMock()
        mock_runner.run.return_value = 0
        mock_runner_cls.return_value = mock_runner
        mock_compose.return_value = _make_task_doc()

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_finished == 1


def test_relax_maker_invalid_code():
    maker = AirssRelaxMaker(code="vasp")
    job = maker.make(
        structures=[_make_si_structure()],
        struct_names=["Si-001"],
        cellinputs=[_make_cellinput()],
        paraminput=_make_paraminput(),
        project_name="test",
        seed_name="Si",
    )

    responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_failed == 1
    assert output.results[0].relax_status == RelaxOutcome.FAILED
    assert "Unknown code" in output.results[0].error_message


# --- AirssValidateMaker tests ---


def test_validate_all_found():
    maker = AirssValidateMaker()
    job = maker.make()

    with patch("airsspy.jf.jobs.subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        responses = run_locally(job, ensure_success=True)

    response = responses[job.uuid][1]
    assert response.output is None


def test_validate_missing_exe():
    maker = AirssValidateMaker()
    job = maker.make()

    with patch("airsspy.jf.jobs.subprocess.run") as mock_run:
        mock_run.side_effect = subprocess.CalledProcessError(1, "which")
        responses = run_locally(job)

    response = responses[job.uuid][1]
    assert response.stop_jobflow is True


def test_validate_additional_exes():
    maker = AirssValidateMaker(additional_exes=("gulp_relax",))
    job = maker.make()

    with patch("airsspy.jf.jobs.subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        responses = run_locally(job, ensure_success=True)

    response = responses[job.uuid][1]
    assert response.output is None
    # buildcell + castep_relax + castep2res + gulp_relax = 4 calls
    assert mock_run.call_count == 4


# --- Crash resilience tests ---


def test_search_maker_exception_continues_loop():
    """Runner crash on one structure should not prevent others from running."""
    maker = AirssSearchMaker(n_structures=3)
    job = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=_make_paraminput(),
        project_name="test",
    )

    with (
        patch("airsspy.jf.jobs.run_buildcell") as mock_buildcell,
        patch("airsspy.jf.jobs.AirssCastepRelaxRunner") as mock_runner_cls,
        patch("airsspy.jf.jobs.compose_task_doc") as mock_compose,
        patch("castepinput.inputs.CellInput") as mock_cellinput_cls,
    ):
        mock_buildcell.return_value = {
            "struct_name": "Si-001",
            "seed_name": "Si",
            "struct_content": "cell",
        }
        mock_runner = MagicMock()
        # Second call raises, first and third succeed
        mock_runner.run.side_effect = [0, RuntimeError("segfault"), 0]
        mock_runner_cls.return_value = mock_runner
        mock_compose.return_value = _make_task_doc()
        mock_cellinput_cls.from_file.return_value = MagicMock()

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_structures == 3
    assert output.n_finished == 2
    assert output.n_failed == 1
    assert mock_buildcell.call_count == 3
    # Second result is the failed one
    assert output.results[1].relax_status == RelaxOutcome.FAILED
    assert "segfault" in output.results[1].error_message
    assert output.results[0].relax_status == RelaxOutcome.FINISHED
    assert output.results[2].relax_status == RelaxOutcome.FINISHED


def test_relax_maker_exception_continues_loop():
    """Runner crash on one structure should not prevent others from running."""
    maker = AirssRelaxMaker(code="castep")
    job = maker.make(
        structures=[_make_si_structure(), _make_si_structure(), _make_si_structure()],
        struct_names=["Si-001", "Si-002", "Si-003"],
        cellinputs=[_make_cellinput(), _make_cellinput(), _make_cellinput()],
        paraminput=_make_paraminput(),
        project_name="test",
        seed_name="Si",
    )

    with (
        patch("airsspy.jf.jobs.AirssCastepRelaxRunner") as mock_runner_cls,
        patch("airsspy.jf.jobs.compose_task_doc") as mock_compose,
    ):
        mock_runner = MagicMock()
        mock_runner.run.side_effect = [0, RuntimeError("disk full"), 0]
        mock_runner_cls.return_value = mock_runner
        mock_compose.return_value = _make_task_doc()

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_structures == 3
    assert output.n_finished == 2
    assert output.n_failed == 1
    assert output.results[1].relax_status == RelaxOutcome.FAILED
    assert "disk full" in output.results[1].error_message


def test_search_maker_buildcell_exception_continues():
    """Exception from run_buildcell itself should be caught."""
    maker = AirssSearchMaker(n_structures=2)
    job = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=_make_paraminput(),
        project_name="test",
    )

    with patch("airsspy.jf.jobs.run_buildcell") as mock_buildcell:
        mock_buildcell.side_effect = [RuntimeError("buildcell not found"), None]

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_failed == 1
    assert output.n_structures == 1
    assert output.results[0].relax_status == RelaxOutcome.FAILED
    assert "buildcell not found" in output.results[0].error_message


def test_search_maker_all_exceptions():
    """All structures failing with exceptions should still produce a valid doc."""
    maker = AirssSearchMaker(n_structures=2, stop_if_not_converged=True)
    job = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=_make_paraminput(),
        project_name="test",
    )

    with (
        patch("airsspy.jf.jobs.run_buildcell") as mock_buildcell,
    ):
        mock_buildcell.side_effect = RuntimeError("total failure")

        responses = run_locally(job, ensure_success=True)

    output = responses[job.uuid][1].output
    assert output.n_structures == 2
    assert output.n_failed == 2
    assert output.n_finished == 0
    # stop_if_not_converged should NOT trigger for FAILED (only ERRORED)
    assert responses[job.uuid][1].stop_children is False

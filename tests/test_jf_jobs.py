"""Tests for jobflow Makers (AirssSearchMaker, AirssRelaxMaker, AirssValidateMaker)."""

from unittest.mock import MagicMock, patch

import pytest

jobflow = pytest.importorskip("jobflow")

from airsspy.jf.documents import RelaxOutcome  # noqa: E402
from airsspy.jf.jobs import (  # noqa: E402
    AirssRelaxMaker,
    AirssSearchMaker,
    AirssValidateMaker,
)


def _make_task_doc(**overrides):
    """Create a default compose_task_doc return dict with optional overrides."""
    defaults = {
        "structure": MagicMock(volume=100.0, reduced_formula="Si"),
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


@pytest.fixture
def mock_paraminput():
    m = MagicMock()
    m.get_string.return_value = "task: geometryoptimization"
    return m


# --- AirssSearchMaker tests ---


@patch("airsspy.jf.jobs.compose_task_doc")
@patch("airsspy.jf.jobs.AirssCastepRelaxRunner")
@patch("airsspy.jf.jobs.run_buildcell")
def test_search_maker_castep_success(mock_buildcell, mock_runner_cls, mock_compose):
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

    maker = AirssSearchMaker(n_structures=1, code="castep")
    response = maker.make(
        seed_name="Si",
        seed_content="seed cell content",
        paraminput=MagicMock(),
        project_name="test_project",
    )

    assert response.output.project_name == "test_project"
    assert response.output.seed_name == "Si"
    assert response.output.n_structures == 1
    assert response.output.n_finished == 1
    assert response.output.n_errored == 0
    assert len(response.output.results) == 1
    assert response.output.results[0].relax_status == RelaxOutcome.FINISHED


@patch("airsspy.jf.jobs.compose_task_doc")
@patch("airsspy.jf.jobs.AirssCastepRelaxRunner")
@patch("airsspy.jf.jobs.run_buildcell")
def test_search_maker_buildcell_timeout(mock_buildcell, mock_runner_cls, mock_compose):
    mock_buildcell.return_value = None

    maker = AirssSearchMaker(n_structures=1)
    response = maker.make(
        seed_name="Si",
        seed_content="seed cell content",
        paraminput=MagicMock(),
        project_name="test",
    )

    assert response.output.n_structures == 0
    assert len(response.output.results) == 0


@patch("airsspy.jf.jobs.compose_task_doc")
@patch("airsspy.jf.jobs.AirssCastepRelaxRunner")
@patch("airsspy.jf.jobs.run_buildcell")
def test_search_maker_relax_error(mock_buildcell, mock_runner_cls, mock_compose):
    mock_buildcell.return_value = {
        "struct_name": "Si-001",
        "seed_name": "Si",
        "struct_content": "cell",
    }
    mock_runner = MagicMock()
    mock_runner.run.return_value = 1
    mock_runner_cls.return_value = mock_runner
    mock_compose.return_value = _make_task_doc()

    maker = AirssSearchMaker(n_structures=1)
    response = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=MagicMock(),
        project_name="test",
    )

    assert response.output.n_errored == 1
    assert response.output.n_finished == 0
    assert response.output.results[0].relax_status == RelaxOutcome.ERRORED


@patch("airsspy.jf.jobs.compose_task_doc")
@patch("airsspy.jf.jobs.AirssGulpRelaxRunner")
@patch("airsspy.jf.jobs.run_buildcell")
def test_search_maker_gulp_success(mock_buildcell, mock_runner_cls, mock_compose):
    mock_buildcell.return_value = {
        "struct_name": "Si-001",
        "seed_name": "Si",
        "struct_content": "cell",
    }
    mock_runner = MagicMock()
    mock_runner.run.return_value = 0
    mock_runner_cls.return_value = mock_runner
    mock_compose.return_value = _make_task_doc()

    maker = AirssSearchMaker(n_structures=1, code="gulp")
    response = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=MagicMock(),
        project_name="test",
    )

    assert response.output.n_finished == 1


@patch("airsspy.jf.jobs.compose_task_doc")
@patch("airsspy.jf.jobs.AirssPp3RelaxRunner")
@patch("airsspy.jf.jobs.run_buildcell")
def test_search_maker_pp3_success(mock_buildcell, mock_runner_cls, mock_compose):
    mock_buildcell.return_value = {
        "struct_name": "Si-001",
        "seed_name": "Si",
        "struct_content": "cell",
    }
    mock_runner = MagicMock()
    mock_runner.run.return_value = 0
    mock_runner_cls.return_value = mock_runner
    mock_compose.return_value = _make_task_doc()

    maker = AirssSearchMaker(n_structures=1, code="pp3")
    response = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=MagicMock(),
        project_name="test",
    )

    assert response.output.n_finished == 1


@patch("airsspy.jf.jobs.compose_task_doc")
@patch("airsspy.jf.jobs.run_buildcell")
def test_search_maker_invalid_code(mock_buildcell, mock_compose):
    mock_buildcell.return_value = {
        "struct_name": "Si-001",
        "seed_name": "Si",
        "struct_content": "cell",
    }
    mock_compose.return_value = _make_task_doc()

    maker = AirssSearchMaker(n_structures=1, code="vasp")
    with pytest.raises(ValueError, match="Unknown code: vasp"):
        maker.make(
            seed_name="Si",
            seed_content="seed",
            paraminput=MagicMock(),
            project_name="test",
        )


@patch("airsspy.jf.jobs.compose_task_doc")
@patch("airsspy.jf.jobs.AirssCastepRelaxRunner")
@patch("airsspy.jf.jobs.run_buildcell")
def test_search_maker_stop_if_all_errored(mock_buildcell, mock_runner_cls, mock_compose):
    mock_buildcell.return_value = {
        "struct_name": "Si-001",
        "seed_name": "Si",
        "struct_content": "cell",
    }
    mock_runner = MagicMock()
    mock_runner.run.return_value = 1
    mock_runner_cls.return_value = mock_runner
    mock_compose.return_value = _make_task_doc()

    maker = AirssSearchMaker(n_structures=1, stop_if_not_converged=True)
    response = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=MagicMock(),
        project_name="test",
    )

    assert response.stop_children is True


@patch("airsspy.jf.jobs.compose_task_doc")
@patch("airsspy.jf.jobs.AirssCastepRelaxRunner")
@patch("airsspy.jf.jobs.run_buildcell")
def test_search_maker_multiple_structures(mock_buildcell, mock_runner_cls, mock_compose):
    """Test that multiple structures are iterated correctly."""
    mock_buildcell.return_value = {
        "struct_name": "Si-001",
        "seed_name": "Si",
        "struct_content": "cell",
    }
    mock_runner = MagicMock()
    mock_runner.run.return_value = 0
    mock_runner_cls.return_value = mock_runner
    mock_compose.return_value = _make_task_doc()

    maker = AirssSearchMaker(n_structures=3)
    response = maker.make(
        seed_name="Si",
        seed_content="seed",
        paraminput=MagicMock(),
        project_name="test",
    )

    assert response.output.n_structures == 3
    assert mock_buildcell.call_count == 3


# --- AirssRelaxMaker tests ---


@patch("airsspy.jf.jobs.compose_task_doc")
@patch("airsspy.jf.jobs.AirssCastepRelaxRunner")
def test_relax_maker_castep_success(mock_runner_cls, mock_compose):
    mock_runner = MagicMock()
    mock_runner.run.return_value = 0
    mock_runner_cls.return_value = mock_runner
    mock_compose.return_value = _make_task_doc()

    structure = MagicMock()
    cellinput = MagicMock()

    maker = AirssRelaxMaker(code="castep")
    response = maker.make(
        structures=[structure],
        struct_names=["Si-001"],
        cellinputs=[cellinput],
        paraminput=MagicMock(),
        project_name="test",
        seed_name="Si",
    )

    assert response.output.n_finished == 1
    assert response.output.job_type == "relax"
    cellinput.set_positions.assert_called_once()
    cellinput.set_cell.assert_called_once()


@patch("airsspy.jf.jobs.compose_task_doc")
@patch("airsspy.jf.jobs.AirssGulpRelaxRunner")
def test_relax_maker_gulp_success(mock_runner_cls, mock_compose):
    mock_runner = MagicMock()
    mock_runner.run.return_value = 0
    mock_runner_cls.return_value = mock_runner
    mock_compose.return_value = _make_task_doc()

    maker = AirssRelaxMaker(code="gulp", executable="ggulp")
    response = maker.make(
        structures=[MagicMock()],
        struct_names=["Si-001"],
        cellinputs=[MagicMock()],
        paraminput=MagicMock(),
        project_name="test",
        seed_name="Si",
    )

    assert response.output.n_finished == 1


@patch("airsspy.jf.jobs.compose_task_doc")
def test_relax_maker_invalid_code(mock_compose):
    mock_compose.return_value = _make_task_doc()

    maker = AirssRelaxMaker(code="vasp")
    with pytest.raises(ValueError, match="Unknown code: vasp"):
        maker.make(
            structures=[MagicMock()],
            struct_names=["Si-001"],
            cellinputs=[MagicMock()],
            paraminput=MagicMock(),
            project_name="test",
            seed_name="Si",
        )


# --- AirssValidateMaker tests ---


@patch("airsspy.jf.jobs.subprocess.run")
def test_validate_all_found(mock_run):
    mock_run.return_value = MagicMock(returncode=0)

    maker = AirssValidateMaker()
    response = maker.make()

    assert response is None


@patch("airsspy.jf.jobs.subprocess.run")
def test_validate_missing_exe(mock_run):
    mock_run.side_effect = FileNotFoundError("not found")

    maker = AirssValidateMaker()
    response = maker.make()

    assert response.stop_jobflow is True


@patch("airsspy.jf.jobs.subprocess.run")
def test_validate_additional_exes(mock_run):
    mock_run.return_value = MagicMock(returncode=0)

    maker = AirssValidateMaker(additional_exes=("gulp_relax",))
    response = maker.make()

    assert response is None
    # buildcell + castep_relax + castep2res + gulp_relax = 4 calls
    assert mock_run.call_count == 4

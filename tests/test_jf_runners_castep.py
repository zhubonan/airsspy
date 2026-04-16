"""Tests for AirssCastepRelaxRunner and compose_task_doc."""

from pathlib import Path
from unittest.mock import MagicMock, patch

# --- compose_task_doc tests ---


CASTEP_OUTPUT_COMPLETE = """\
Starting spin               0.0 s
Initial spin               1.0 s
NB est. 0K energy         -10.123456 eV
Final energy, E             -10.123456 eV
Final free energy (E-TS)   -10.456789 eV
Pressure: 0.123 GPa
Overall parallel efficiency   85%
Total time            = 120.3 s
 Spin density    0.00000 mu_B
 Total  Charge(e)   Spin(hbar/2)
     Si        0.000000        0.000
     Si        0.000000        0.000
|Spin density   0.00000 mu_B
Geometry optimization completed
Finished iteration 50
"""


CASTEP_OUTPUT_FAILED = """\
Starting spin               0.0 s
NB est. 0K energy         -10.123456 eV
Final energy, E             -10.123456 eV
Geometry optimization failed
Finished iteration 50
"""


CASTEP_OUTPUT_NO_ENERGY = """\
Starting spin               0.0 s
Geometry optimization completed
Finished iteration 50
Total time            = 120.3 s
"""


CASTEP_OUTPUT_WITH_SPIN = """\
NB est. 0K energy         -10.123456 eV
Pressure: 0.1 GPa
 Total  Charge(e)   Spin(hbar/2)
     Si        0.000000        0.012
     Si        0.000000       -0.012
 Length (A)
"""


CELL_OUT_CONTENT = """\
%BLOCK LATTICE_CART
  5.430000  0.000000  0.000000
  0.000000  5.430000  0.000000
  0.000000  0.000000  5.430000
%ENDBLOCK LATTICE_CART
%BLOCK POSITIONS_ABS
Si   0.000000  0.000000  0.000000
Si   2.715000  2.715000  2.715000
%ENDBLOCK POSITIONS_ABS
"""


@patch("airsspy.restools.save_airss_res")
@patch("castepinput.inputs.CellInput")
def test_compose_task_doc_basic(mock_cellinput, mock_save):
    """Test basic compose_task_doc with mocked CellInput (no real files)."""
    from airsspy.jf.runners import compose_task_doc

    # Setup mock for CellInput
    mock_ci = MagicMock()
    mock_ci.get_positions.return_value = (["Si", "Si"], [[0, 0, 0], [2.715, 2.715, 2.715]], [None, None])
    mock_ci.get_cell.return_value = [[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]]
    mock_cellinput.from_file.return_value = mock_ci

    doc = compose_task_doc("nonexistent")

    # No .castep file so energy/pressure stay None
    assert doc["energy"] is None
    assert doc["pressure"] is None
    assert doc["natoms"] == 2
    mock_save.assert_called_once()


@patch("airsspy.restools.save_airss_res")
@patch("castepinput.inputs.CellInput")
def test_compose_task_doc_with_real_files(mock_cellinput, mock_save, tmp_path, monkeypatch):
    """Test compose_task_doc reads real files from disk."""
    from airsspy.jf.runners import compose_task_doc

    monkeypatch.chdir(tmp_path)
    Path("Si-001.castep").write_text(CASTEP_OUTPUT_COMPLETE)
    Path("Si-001-out.cell").write_text(CELL_OUT_CONTENT)

    mock_ci = MagicMock()
    mock_ci.get_positions.return_value = (["Si", "Si"], [[0, 0, 0], [2.715, 2.715, 2.715]], [None, None])
    mock_ci.get_cell.return_value = [[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]]
    mock_cellinput.from_file.return_value = mock_ci

    doc = compose_task_doc("Si-001")

    assert doc["energy"] == -10.123456
    assert doc["pressure"] == 0.123
    assert doc["parallel_efficiency"] == 0.85
    assert doc["total_time"] == 120.3
    assert doc["natoms"] == 2
    mock_save.assert_called_once()


@patch("airsspy.restools.save_airss_res")
@patch("castepinput.inputs.CellInput")
def test_compose_task_doc_fallback_cell(mock_cellinput, mock_save, tmp_path, monkeypatch):
    """Test compose_task_doc falls back to .cell when -out.cell is missing."""
    from airsspy.jf.runners import compose_task_doc

    monkeypatch.chdir(tmp_path)
    Path("Si-001.castep").write_text(CASTEP_OUTPUT_COMPLETE)
    Path("Si-001.cell").write_text(CELL_OUT_CONTENT)

    mock_ci = MagicMock()
    mock_ci.get_positions.return_value = (["Si", "Si"], [[0, 0, 0], [2.715, 2.715, 2.715]], [None, None])
    mock_ci.get_cell.return_value = [[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]]
    mock_cellinput.from_file.return_value = mock_ci

    doc = compose_task_doc("Si-001")

    assert doc["energy"] == -10.123456
    # No -out.cell on disk, so from_file called once with .cell
    mock_cellinput.from_file.assert_called_once_with("Si-001.cell")


@patch("airsspy.restools.save_airss_res")
@patch("castepinput.inputs.CellInput")
def test_compose_task_doc_missing_castep(mock_cellinput, mock_save, tmp_path, monkeypatch):
    """Test compose_task_doc handles missing .castep file."""
    from airsspy.jf.runners import compose_task_doc

    monkeypatch.chdir(tmp_path)
    Path("Si-001-out.cell").write_text(CELL_OUT_CONTENT)

    mock_ci = MagicMock()
    mock_ci.get_positions.return_value = (["Si", "Si"], [[0, 0, 0], [2.715, 2.715, 2.715]], [None, None])
    mock_ci.get_cell.return_value = [[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]]
    mock_cellinput.from_file.return_value = mock_ci

    # Should not raise, just return None for missing fields
    doc = compose_task_doc("Si-001")
    assert doc["energy"] is None
    assert doc["pressure"] is None


@patch("airsspy.restools.save_airss_res")
@patch("castepinput.inputs.CellInput")
def test_compose_task_doc_with_spin(mock_cellinput, mock_save, tmp_path, monkeypatch):
    """Test compose_task_doc extracts spin moments."""
    from airsspy.jf.runners import compose_task_doc

    monkeypatch.chdir(tmp_path)
    Path("Si-001.castep").write_text(CASTEP_OUTPUT_WITH_SPIN)
    Path("Si-001-out.cell").write_text(CELL_OUT_CONTENT)

    mock_ci = MagicMock()
    mock_ci.get_positions.return_value = (["Si", "Si"], [[0, 0, 0], [2.715, 2.715, 2.715]], [None, None])
    mock_ci.get_cell.return_value = [[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]]
    mock_cellinput.from_file.return_value = mock_ci

    doc = compose_task_doc("Si-001")
    assert doc["spin"] == 0.0


# --- AirssCastepRelaxRunner tests ---


@patch("castepinput.inputs.CellInput")
@patch("airsspy.jf.runners.subprocess.run")
def test_castep_runner_converged(mock_run, mock_cellinput, tmp_path, monkeypatch):
    """Test runner returns 0 after two consecutive completed optimizations."""
    from castepinput.inputs import ParamInput

    from airsspy.jf.runners import AirssCastepRelaxRunner

    monkeypatch.chdir(tmp_path)

    castep_completed = "Geometry optimization completed\nFinished iteration 50\nTotal time 10.0 s\n"

    mock_run.return_value = MagicMock(returncode=0)
    mock_ci = MagicMock()
    mock_ci.get_cell.return_value = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    mock_ci.get_positions.return_value = (["Si"], [[0, 0, 0]], [None])
    mock_cellinput.from_file.return_value = mock_ci

    def write_castep_side_effect(*args, **kwargs):
        with open("Si-001.castep", "w") as f:
            f.write(castep_completed)
        return MagicMock(returncode=0)

    mock_run.side_effect = write_castep_side_effect

    param = ParamInput()
    runner = AirssCastepRelaxRunner(executable="castep", cycles=3, max_fails=2, max_iterations=200)
    result = runner.run("Si-001", "cell content", param)

    assert result == 0


@patch("castepinput.inputs.CellInput")
@patch("airsspy.jf.runners.subprocess.run")
def test_castep_runner_max_iterations(mock_run, mock_cellinput, tmp_path, monkeypatch):
    """Test runner returns 1 when max iterations exceeded."""
    from castepinput.inputs import ParamInput

    from airsspy.jf.runners import AirssCastepRelaxRunner

    monkeypatch.chdir(tmp_path)
    mock_run.return_value = MagicMock(returncode=0)
    mock_ci = MagicMock()
    mock_ci.get_cell.return_value = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    mock_ci.get_positions.return_value = (["Si"], [[0, 0, 0]], [None])
    mock_cellinput.from_file.return_value = mock_ci

    call_count = [0]

    def write_castep_side_effect(*args, **kwargs):
        call_count[0] += 1
        with open("Si-001.castep", "w") as f:
            if call_count[0] < 3:
                f.write("Geometry optimization\nFinished iteration 99\n")
            else:
                f.write("Geometry optimization completed\nFinished iteration 150\nTotal time 10.0 s\n")
        return MagicMock(returncode=0)

    mock_run.side_effect = write_castep_side_effect

    param = ParamInput()
    runner = AirssCastepRelaxRunner(executable="castep", cycles=5, max_iterations=200)
    result = runner.run("Si-001", "cell content", param)

    assert result == 1


@patch("castepinput.inputs.CellInput")
@patch("airsspy.jf.runners.subprocess.run")
def test_castep_runner_max_fails(mock_run, mock_cellinput, tmp_path, monkeypatch):
    """Test runner returns 1 after too many consecutive failures."""
    from castepinput.inputs import ParamInput

    from airsspy.jf.runners import AirssCastepRelaxRunner

    monkeypatch.chdir(tmp_path)
    mock_run.return_value = MagicMock(returncode=1)
    mock_ci = MagicMock()
    mock_ci.get_cell.return_value = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    mock_ci.get_positions.return_value = (["Si"], [[0, 0, 0]], [None])
    mock_cellinput.from_file.return_value = mock_ci

    param = ParamInput()
    runner = AirssCastepRelaxRunner(executable="castep", cycles=5, max_fails=2)
    result = runner.run("Si-001", "cell content", param)

    assert result == 1
    assert mock_run.call_count == 3  # initial + 2 retries

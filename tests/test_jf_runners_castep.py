"""Tests for AirssCastepRelaxRunner and compose_task_doc."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

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


def test_run_buildcell_writes_original_seed_when_transformed(tmp_path, monkeypatch):
    """Formula sampling should not mutate the source seed file."""
    from airsspy.jf.runners import run_buildcell

    monkeypatch.chdir(tmp_path)
    proc = MagicMock()
    proc.communicate.return_value = ("generated cell", "")

    with (
        patch("airsspy.jf.runners.subprocess.Popen", return_value=proc),
        patch("airsspy.casteptools.get_rand_cell_name", return_value="Si-001.cell"),
    ):
        result = run_buildcell(
            "Si",
            "original seed",
            seed_text_transform=lambda text: text.replace("original", "sampled"),
        )

    assert result["struct_name"] == "Si-001"
    assert Path("Si.cell").read_text() == "original seed"
    assert Path("Si-001.cell").read_text() == "generated cell"
    proc.communicate.assert_called_once_with("sampled seed", timeout=30)


@patch("airsspy.restools.save_airss_res")
@patch("castepinput.inputs.CellInput")
def test_compose_task_doc_basic(mock_cellinput, mock_save):
    """Test basic compose_task_doc with mocked CellInput (no real files)."""
    from airsspy.jf.runners import compose_task_doc

    # Setup mock for CellInput
    mock_ci = MagicMock()
    mock_ci.get_positions.return_value = (
        ["Si", "Si"],
        [[0, 0, 0], [2.715, 2.715, 2.715]],
        [None, None],
    )
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
def test_compose_task_doc_with_real_files(
    mock_cellinput, mock_save, tmp_path, monkeypatch
):
    """Test compose_task_doc reads real files from disk."""
    from airsspy.jf.runners import compose_task_doc

    monkeypatch.chdir(tmp_path)
    Path("Si-001.castep").write_text(CASTEP_OUTPUT_COMPLETE)
    Path("Si-001-out.cell").write_text(CELL_OUT_CONTENT)

    mock_ci = MagicMock()
    mock_ci.get_positions.return_value = (
        ["Si", "Si"],
        [[0, 0, 0], [2.715, 2.715, 2.715]],
        [None, None],
    )
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
def test_compose_task_doc_fallback_cell(
    mock_cellinput, mock_save, tmp_path, monkeypatch
):
    """Test compose_task_doc falls back to .cell when -out.cell is missing."""
    from airsspy.jf.runners import compose_task_doc

    monkeypatch.chdir(tmp_path)
    Path("Si-001.castep").write_text(CASTEP_OUTPUT_COMPLETE)
    Path("Si-001.cell").write_text(CELL_OUT_CONTENT)

    mock_ci = MagicMock()
    mock_ci.get_positions.return_value = (
        ["Si", "Si"],
        [[0, 0, 0], [2.715, 2.715, 2.715]],
        [None, None],
    )
    mock_ci.get_cell.return_value = [[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]]
    mock_cellinput.from_file.return_value = mock_ci

    doc = compose_task_doc("Si-001")

    assert doc["energy"] == -10.123456
    # No -out.cell on disk, so from_file called once with .cell
    mock_cellinput.from_file.assert_called_once_with("Si-001.cell")


@patch("airsspy.restools.save_airss_res")
@patch("castepinput.inputs.CellInput")
def test_compose_task_doc_missing_castep(
    mock_cellinput, mock_save, tmp_path, monkeypatch
):
    """Test compose_task_doc handles missing .castep file."""
    from airsspy.jf.runners import compose_task_doc

    monkeypatch.chdir(tmp_path)
    Path("Si-001-out.cell").write_text(CELL_OUT_CONTENT)

    mock_ci = MagicMock()
    mock_ci.get_positions.return_value = (
        ["Si", "Si"],
        [[0, 0, 0], [2.715, 2.715, 2.715]],
        [None, None],
    )
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
    mock_ci.get_positions.return_value = (
        ["Si", "Si"],
        [[0, 0, 0], [2.715, 2.715, 2.715]],
        [None, None],
    )
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

    castep_completed = (
        "Geometry optimization completed\nFinished iteration 50\nTotal time 10.0 s\n"
    )

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
    runner = AirssCastepRelaxRunner(
        executable="castep", max_fails=2, max_iterations=200
    )
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
                f.write(
                    "Geometry optimization completed\nFinished iteration 150\nTotal time 10.0 s\n"
                )
        return MagicMock(returncode=0)

    mock_run.side_effect = write_castep_side_effect

    param = ParamInput()
    runner = AirssCastepRelaxRunner(executable="castep", max_iterations=200)
    result = runner.run("Si-001", "cell content", param)

    assert result == 1


@patch("castepinput.inputs.CellInput")
@patch("airsspy.jf.runners.subprocess.run")
def test_castep_runner_counts_lowercase_lbfgs_iterations(
    mock_run, mock_cellinput, tmp_path, monkeypatch
):
    """Lowercase CASTEP LBFGS iteration lines should advance max_iterations."""
    from castepinput.inputs import ParamInput

    from airsspy.jf.runners import AirssCastepRelaxRunner

    monkeypatch.chdir(tmp_path)
    mock_run.return_value = MagicMock(returncode=0)
    mock_ci = MagicMock()
    mock_ci.get_cell.return_value = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    mock_ci.get_positions.return_value = (["Si"], [[0, 0, 0]], [None])
    mock_cellinput.from_file.return_value = mock_ci

    def write_castep_side_effect(*args, **kwargs):
        with open("Si-001.castep", "w") as f:
            f.write(
                "LBFGS: finished iteration     0 with enthalpy= -1.0E+000 eV\n"
                "LBFGS: finished iteration     1 with enthalpy= -2.0E+000 eV\n"
                "LBFGS: finished iteration     2 with enthalpy= -3.0E+000 eV\n"
                "LBFGS: WARNING - Geometry optimization failed to converge "
                "after          2 steps\n"
            )
        return MagicMock(returncode=0)

    mock_run.side_effect = write_castep_side_effect

    param = ParamInput()
    runner = AirssCastepRelaxRunner(executable="castep", max_iterations=6)
    result = runner.run("Si-001", "cell content", param)

    assert result == 1
    assert mock_run.call_count == 3


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
    runner = AirssCastepRelaxRunner(executable="castep", max_fails=2)
    result = runner.run("Si-001", "cell content", param)

    assert result == 1
    assert mock_run.call_count == 3  # initial + 2 retries


# --- Restart copy tests: -out.cell -> .cell block update ---

INITIAL_CELL = """\
%BLOCK LATTICE_CART
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_ABS
Si 0.0 0.0 0.0
%ENDBLOCK POSITIONS_ABS

kpoints_mp_grid : 4 4 4
"""

UPDATED_OUT_CELL = """\
%BLOCK LATTICE_CART
5.0 0.0 0.0
0.0 5.0 0.0
0.0 0.0 5.0
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_ABS
Si 0.0 0.0 0.0
Si 2.5 2.5 2.5
%ENDBLOCK POSITIONS_ABS
"""

OUT_CELL_WITH_LATTICE_ABC = """\
%BLOCK LATTICE_ABC
5.0 5.0 5.0
90.0 90.0 90.0
%ENDBLOCK LATTICE_ABC

%BLOCK POSITIONS_ABS
Si 0.0 0.0 0.0
Si 2.5 2.5 2.5
%ENDBLOCK POSITIONS_ABS
"""

OUT_CELL_WITH_POSITIONS_FRAC = """\
%BLOCK LATTICE_CART
5.0 0.0 0.0
0.0 5.0 0.0
0.0 0.0 5.0
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Si 0.0 0.0 0.0
Si 0.5 0.5 0.5
%ENDBLOCK POSITIONS_FRAC
"""

OUT_CELL_MISSING_LATTICE = """\
%BLOCK POSITIONS_ABS
Si 0.0 0.0 0.0
%ENDBLOCK POSITIONS_ABS
"""

OUT_CELL_MISSING_POSITIONS = """\
%BLOCK LATTICE_CART
5.0 0.0 0.0
0.0 5.0 0.0
0.0 0.0 5.0
%ENDBLOCK LATTICE_CART
"""

CASTEP_CONVERGED = (
    "Geometry optimization completed\nFinished iteration 50\nTotal time 10.0 s\n"
)


@patch("airsspy.jf.runners.subprocess.run")
def test_restart_copy_basic(mock_run, tmp_path, monkeypatch):
    """Restart copy updates lattice and positions from -out.cell."""
    from castepinput.inputs import ParamInput

    from airsspy.jf.runners import AirssCastepRelaxRunner

    monkeypatch.chdir(tmp_path)
    Path("Si-001.cell").write_text(INITIAL_CELL)
    Path("Si-001-out.cell").write_text(UPDATED_OUT_CELL)

    def side_effect(*args, **kwargs):
        Path("Si-001.castep").write_text(CASTEP_CONVERGED)
        return MagicMock(returncode=0)

    mock_run.side_effect = side_effect

    runner = AirssCastepRelaxRunner(executable="castep", max_iterations=200)
    runner.run("Si-001", INITIAL_CELL, ParamInput())

    result = Path("Si-001.cell").read_text()
    assert "5.0 0.0 0.0" in result
    assert "Si 2.5 2.5 2.5" in result
    assert "kpoints_mp_grid" in result
    assert result.count("%BLOCK lattice") == 1
    assert result.count("%BLOCK positions") == 1


@patch("airsspy.jf.runners.subprocess.run")
def test_restart_preserves_other_content(mock_run, tmp_path, monkeypatch):
    """Non-structural content from original .cell is preserved."""
    from castepinput.inputs import ParamInput

    from airsspy.jf.runners import AirssCastepRelaxRunner

    cell_with_extra = INITIAL_CELL + "\nspecies_pot : Si POT\nsymmetry_tol : 0.01\n"
    extra_out = UPDATED_OUT_CELL + "\nsymmetry_tol : 0.001\n"

    monkeypatch.chdir(tmp_path)
    Path("Si-001.cell").write_text(cell_with_extra)
    Path("Si-001-out.cell").write_text(extra_out)

    def side_effect(*args, **kwargs):
        Path("Si-001.castep").write_text(CASTEP_CONVERGED)
        return MagicMock(returncode=0)

    mock_run.side_effect = side_effect

    runner = AirssCastepRelaxRunner(executable="castep", max_iterations=200)
    runner.run("Si-001", cell_with_extra, ParamInput())

    result = Path("Si-001.cell").read_text()
    assert "species_pot" in result
    assert "Si POT" in result
    assert "kpoints_mp_grid" in result
    # lattice/positions from -out.cell are present
    assert "5.0 0.0 0.0" in result
    assert "Si 2.5 2.5 2.5" in result


@patch("airsspy.jf.runners.subprocess.run")
def test_restart_lattice_abc_replaces_cart(mock_run, tmp_path, monkeypatch):
    """Output LATTICE_ABC replaces input LATTICE_CART."""
    from castepinput.inputs import ParamInput

    from airsspy.jf.runners import AirssCastepRelaxRunner

    monkeypatch.chdir(tmp_path)
    Path("Si-001.cell").write_text(INITIAL_CELL)
    Path("Si-001-out.cell").write_text(OUT_CELL_WITH_LATTICE_ABC)

    def side_effect(*args, **kwargs):
        Path("Si-001.castep").write_text(CASTEP_CONVERGED)
        return MagicMock(returncode=0)

    mock_run.side_effect = side_effect

    runner = AirssCastepRelaxRunner(executable="castep", max_iterations=200)
    runner.run("Si-001", INITIAL_CELL, ParamInput())

    result = Path("Si-001.cell").read_text()
    assert "%BLOCK lattice_abc" in result
    assert "lattice_cart" not in result
    assert "5.0 5.0 5.0" in result
    assert "kpoints_mp_grid" in result


@patch("airsspy.jf.runners.subprocess.run")
def test_restart_lattice_cart_replaces_abc(mock_run, tmp_path, monkeypatch):
    """Output LATTICE_CART replaces input LATTICE_ABC."""
    from castepinput.inputs import ParamInput

    from airsspy.jf.runners import AirssCastepRelaxRunner

    cell_with_abc = """\
%BLOCK LATTICE_ABC
1.0 1.0 1.0
90.0 90.0 90.0
%ENDBLOCK LATTICE_ABC
%BLOCK POSITIONS_ABS
Si 0.0 0.0 0.0
%ENDBLOCK POSITIONS_ABS
"""
    monkeypatch.chdir(tmp_path)
    Path("Si-001.cell").write_text(cell_with_abc)
    Path("Si-001-out.cell").write_text(UPDATED_OUT_CELL)

    def side_effect(*args, **kwargs):
        Path("Si-001.castep").write_text(CASTEP_CONVERGED)
        return MagicMock(returncode=0)

    mock_run.side_effect = side_effect

    runner = AirssCastepRelaxRunner(executable="castep", max_iterations=200)
    runner.run("Si-001", cell_with_abc, ParamInput())

    result = Path("Si-001.cell").read_text()
    assert "%BLOCK lattice_cart" in result
    assert "lattice_abc" not in result


@patch("airsspy.jf.runners.subprocess.run")
def test_restart_positions_frac_replaces_abs(mock_run, tmp_path, monkeypatch):
    """Output POSITIONS_FRAC replaces input POSITIONS_ABS."""
    from castepinput.inputs import ParamInput

    from airsspy.jf.runners import AirssCastepRelaxRunner

    monkeypatch.chdir(tmp_path)
    Path("Si-001.cell").write_text(INITIAL_CELL)
    Path("Si-001-out.cell").write_text(OUT_CELL_WITH_POSITIONS_FRAC)

    def side_effect(*args, **kwargs):
        Path("Si-001.castep").write_text(CASTEP_CONVERGED)
        return MagicMock(returncode=0)

    mock_run.side_effect = side_effect

    runner = AirssCastepRelaxRunner(executable="castep", max_iterations=200)
    runner.run("Si-001", INITIAL_CELL, ParamInput())

    result = Path("Si-001.cell").read_text()
    assert "%BLOCK positions_frac" in result
    assert "positions_abs" not in result
    assert "Si 0.5 0.5 0.5" in result
    assert "kpoints_mp_grid" in result


@patch("airsspy.jf.runners.subprocess.run")
def test_restart_missing_lattice_raises(mock_run, tmp_path, monkeypatch):
    """Restart raises RuntimeError when -out.cell has no lattice block."""
    from castepinput.inputs import ParamInput

    from airsspy.jf.runners import AirssCastepRelaxRunner

    monkeypatch.chdir(tmp_path)
    Path("Si-001.cell").write_text(INITIAL_CELL)
    Path("Si-001-out.cell").write_text(OUT_CELL_MISSING_LATTICE)

    def side_effect(*args, **kwargs):
        Path("Si-001.castep").write_text(CASTEP_CONVERGED)
        return MagicMock(returncode=0)

    mock_run.side_effect = side_effect

    runner = AirssCastepRelaxRunner(executable="castep", max_iterations=200)
    with pytest.raises(RuntimeError, match="No lattice block"):
        runner.run("Si-001", INITIAL_CELL, ParamInput())


@patch("airsspy.jf.runners.subprocess.run")
def test_restart_missing_positions_raises(mock_run, tmp_path, monkeypatch):
    """Restart raises RuntimeError when -out.cell has no positions block."""
    from castepinput.inputs import ParamInput

    from airsspy.jf.runners import AirssCastepRelaxRunner

    monkeypatch.chdir(tmp_path)
    Path("Si-001.cell").write_text(INITIAL_CELL)
    Path("Si-001-out.cell").write_text(OUT_CELL_MISSING_POSITIONS)

    def side_effect(*args, **kwargs):
        Path("Si-001.castep").write_text(CASTEP_CONVERGED)
        return MagicMock(returncode=0)

    mock_run.side_effect = side_effect

    runner = AirssCastepRelaxRunner(executable="castep", max_iterations=200)
    with pytest.raises(RuntimeError, match="No positions block"):
        runner.run("Si-001", INITIAL_CELL, ParamInput())

"""Tests for ML interatomic potential runners."""

import sys
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

# ---------------------------------------------------------------------------
# _resolve_calculator tests
# ---------------------------------------------------------------------------


class TestResolveCalculator:
    """Tests for _resolve_calculator()."""

    def test_colon_syntax(self):
        """module:Class format."""
        from airsspy.jf.ml_runners import _resolve_calculator

        with patch("importlib.import_module") as mock_import:
            mock_cls = MagicMock(return_value="calculator_instance")
            mock_mod = MagicMock()
            mock_mod.SomeCalc = mock_cls
            mock_import.return_value = mock_mod

            result = _resolve_calculator("my.module:SomeCalc")
            mock_import.assert_called_with("my.module")
            mock_cls.assert_called_once_with()
            assert result == "calculator_instance"

    def test_colon_syntax_with_model(self):
        """module:Class@model format."""
        from airsspy.jf.ml_runners import _resolve_calculator

        with patch("importlib.import_module") as mock_import:
            mock_cls = MagicMock(return_value="calculator_instance")
            mock_mod = MagicMock()
            mock_mod.SomeCalc = mock_cls
            mock_import.return_value = mock_mod

            result = _resolve_calculator("my.module:SomeCalc@medium")
            mock_cls.assert_called_once_with("medium")
            assert result == "calculator_instance"

    def test_colon_syntax_with_path_model(self):
        """module:Class@/path/to/model format."""
        from airsspy.jf.ml_runners import _resolve_calculator

        with patch("importlib.import_module") as mock_import:
            mock_cls = MagicMock(return_value="calculator_instance")
            mock_mod = MagicMock()
            mock_mod.SomeCalc = mock_cls
            mock_import.return_value = mock_mod

            result = _resolve_calculator("my.module:SomeCalc@/path/to/model.pt")
            mock_cls.assert_called_once_with("/path/to/model.pt")
            assert result == "calculator_instance"

    def test_dot_syntax(self):
        """module.Class format."""
        from airsspy.jf.ml_runners import _resolve_calculator

        with patch("importlib.import_module") as mock_import:
            mock_cls = MagicMock(return_value="calculator_instance")
            mock_mod = MagicMock()
            mock_mod.SomeCalc = mock_cls
            mock_import.return_value = mock_mod

            result = _resolve_calculator("my.module.SomeCalc")
            mock_import.assert_called_with("my.module")
            mock_cls.assert_called_once_with()
            assert result == "calculator_instance"

    def test_kwargs_passed(self):
        """Extra kwargs forwarded to constructor."""
        from airsspy.jf.ml_runners import _resolve_calculator

        with patch("importlib.import_module") as mock_import:
            mock_cls = MagicMock(return_value="calculator_instance")
            mock_mod = MagicMock()
            mock_mod.SomeCalc = mock_cls
            mock_import.return_value = mock_mod

            result = _resolve_calculator("my.module:SomeCalc", device="cpu")
            mock_cls.assert_called_once_with(device="cpu")
            assert result == "calculator_instance"

    def test_model_and_kwargs(self):
        """Both @model and kwargs."""
        from airsspy.jf.ml_runners import _resolve_calculator

        with patch("importlib.import_module") as mock_import:
            mock_cls = MagicMock(return_value="calculator_instance")
            mock_mod = MagicMock()
            mock_mod.SomeCalc = mock_cls
            mock_import.return_value = mock_mod

            result = _resolve_calculator("my.module:SomeCalc@medium", device="cpu")
            mock_cls.assert_called_once_with("medium", device="cpu")
            assert result == "calculator_instance"

    def test_invalid_spec_raises(self):
        from airsspy.jf.ml_runners import _resolve_calculator

        with pytest.raises(ValueError, match="Cannot parse"):
            _resolve_calculator("nocolonordot")


# ---------------------------------------------------------------------------
# AirssMlSinglePointRunner tests
# ---------------------------------------------------------------------------


class TestMlSinglePointRunner:
    """Tests for AirssMlSinglePointRunner."""

    def _make_cell_file(self, tmp_path, name="test-001"):
        """Create a minimal .cell file."""
        cell_content = (
            "%BLOCK LATTICE_CART\n"
            "  5.0  0.0  0.0\n"
            "  0.0  5.0  0.0\n"
            "  0.0  0.0  5.0\n"
            "%ENDBLOCK LATTICE_CART\n"
            "%BLOCK POSITIONS_FRAC\n"
            "Si  0.0  0.0  0.0\n"
            "Si  0.25  0.25  0.25\n"
            "%ENDBLOCK POSITIONS_FRAC\n"
        )
        cell_file = tmp_path / (name + ".cell")
        cell_file.write_text(cell_content)
        return str(name), cell_content

    @patch("airsspy.jf.ml_runners.ase_write")
    @patch("airsspy.jf.ml_runners._cell_to_atoms")
    @patch("airsspy.jf.ml_runners._resolve_calculator")
    def test_sp_success(self, mock_resolve, mock_cell_to, mock_write, tmp_path):
        from airsspy.jf.ml_runners import AirssMlSinglePointRunner

        atoms = Atoms(
            "Si2",
            positions=[[0, 0, 0], [1.25, 1.25, 1.25]],
            cell=[[5, 0, 0], [0, 5, 0], [0, 0, 5]],
            pbc=True,
        )

        # Mock the calculator as a dict-like that ASE can use
        # Patch get_potential_energy etc on the atoms level
        atoms.get_potential_energy = MagicMock(return_value=-10.5)
        atoms.get_forces = MagicMock(return_value=np.zeros((2, 3)))
        atoms.get_stress = MagicMock(return_value=np.zeros(6))

        mock_cell_to.return_value = atoms
        mock_resolve.return_value = MagicMock()

        name, cell_content = self._make_cell_file(tmp_path)

        import os

        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            runner = AirssMlSinglePointRunner("fake:Calc@medium")
            rc = runner.run(name, cell_content)
            assert rc == 0
            assert atoms.info["extern_pressure"] == pytest.approx(0.0)
            mock_write.assert_called()
        finally:
            os.chdir(orig)

    @patch("airsspy.jf.ml_runners._resolve_calculator")
    def test_sp_failure(self, mock_resolve, tmp_path):
        from airsspy.jf.ml_runners import AirssMlSinglePointRunner

        mock_resolve.side_effect = ImportError("no such module")

        name, cell_content = self._make_cell_file(tmp_path)

        import os

        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            # Mock _cell_to_atoms to avoid ASE trying to read .cell file
            with patch("airsspy.jf.ml_runners._cell_to_atoms") as mock_cta:
                atoms = Atoms(
                    "Si2",
                    positions=[[0, 0, 0], [1.25, 1.25, 1.25]],
                    cell=[[5, 0, 0], [0, 5, 0], [0, 0, 5]],
                    pbc=True,
                )
                mock_cta.return_value = atoms

                runner = AirssMlSinglePointRunner("fake:Calc")
                rc = runner.run(name, cell_content)
                assert rc == 1
        finally:
            os.chdir(orig)


# ---------------------------------------------------------------------------
# AirssMlRelaxRunner tests
# ---------------------------------------------------------------------------


class TestMlRelaxRunner:
    """Tests for AirssMlRelaxRunner."""

    def _make_cell_file(self, tmp_path, name="test-001"):
        """Create a minimal .cell file."""
        cell_content = (
            "%BLOCK LATTICE_CART\n"
            "  5.0  0.0  0.0\n"
            "  0.0  5.0  0.0\n"
            "  0.0  0.0  5.0\n"
            "%ENDBLOCK LATTICE_CART\n"
            "%BLOCK POSITIONS_FRAC\n"
            "Si  0.0  0.0  0.0\n"
            "Si  0.25  0.25  0.25\n"
            "%ENDBLOCK POSITIONS_FRAC\n"
        )
        return str(name), cell_content

    @patch("airsspy.jf.ml_runners.ase_write")
    @patch("airsspy.jf.ml_runners._cell_to_atoms")
    @patch("airsspy.jf.ml_runners._resolve_calculator")
    def test_relax_success(self, mock_resolve, mock_cell_to, mock_write, tmp_path):
        from airsspy.jf.ml_runners import AirssMlRelaxRunner

        atoms = Atoms(
            "Si2",
            positions=[[0, 0, 0], [1.25, 1.25, 1.25]],
            cell=[[5, 0, 0], [0, 5, 0], [0, 0, 5]],
            pbc=True,
        )
        atoms.get_potential_energy = MagicMock(return_value=-10.5)
        atoms.get_forces = MagicMock(return_value=np.zeros((2, 3)))
        atoms.get_stress = MagicMock(return_value=np.zeros(6))

        mock_cell_to.return_value = atoms
        mock_resolve.return_value = MagicMock()

        name, cell_content = self._make_cell_file(tmp_path)

        import os

        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            with patch("ase.optimize.FIRE") as mock_fire_cls:
                mock_dyn = MagicMock()
                mock_fire_cls.return_value = mock_dyn

                runner = AirssMlRelaxRunner("fake:Calc@medium")
                rc = runner.run(name, cell_content)
                assert rc == 0  # converged
                mock_dyn.run.assert_called_once()
        finally:
            os.chdir(orig)

    @patch("airsspy.jf.ml_runners.ase_write")
    @patch("airsspy.jf.ml_runners._cell_to_atoms")
    @patch("airsspy.jf.ml_runners._resolve_calculator")
    def test_relax_with_pressure(
        self, mock_resolve, mock_cell_to, mock_write, tmp_path
    ):
        from airsspy.jf.ml_runners import AirssMlRelaxRunner

        atoms = Atoms(
            "Si2",
            positions=[[0, 0, 0], [1.25, 1.25, 1.25]],
            cell=[[5, 0, 0], [0, 5, 0], [0, 0, 5]],
            pbc=True,
        )
        atoms.get_potential_energy = MagicMock(return_value=-10.5)
        atoms.get_forces = MagicMock(return_value=np.zeros((2, 3)))
        atoms.get_stress = MagicMock(return_value=np.zeros(6))

        mock_cell_to.return_value = atoms
        mock_resolve.return_value = MagicMock()

        name, cell_content = self._make_cell_file(tmp_path)

        import os

        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            with patch("ase.optimize.FIRE") as mock_fire_cls:
                with patch("ase.filters.ExpCellFilter") as mock_ecf:
                    mock_dyn = MagicMock()
                    mock_fire_cls.return_value = mock_dyn

                    runner = AirssMlRelaxRunner("fake:Calc@medium", pressure=10.0)
                    rc = runner.run(name, cell_content)
                    assert rc == 0
                    mock_ecf.assert_called_once()
        finally:
            os.chdir(orig)

    @patch("airsspy.jf.ml_runners.ase_write")
    @patch("airsspy.jf.ml_runners._cell_to_atoms")
    @patch("airsspy.jf.ml_runners._resolve_calculator")
    def test_relax_max_steps_is_collected_with_metadata(
        self, mock_resolve, mock_cell_to, mock_write, tmp_path
    ):
        from airsspy.jf.ml_runners import AirssMlRelaxRunner

        atoms = Atoms(
            "Si2",
            positions=[[0, 0, 0], [1.25, 1.25, 1.25]],
            cell=[[5, 0, 0], [0, 5, 0], [0, 0, 5]],
            pbc=True,
        )
        atoms.get_potential_energy = MagicMock(return_value=-10.5)
        atoms.get_forces = MagicMock(return_value=np.zeros((2, 3)))
        atoms.get_stress = MagicMock(return_value=np.zeros(6))

        mock_cell_to.return_value = atoms
        mock_resolve.return_value = MagicMock()

        name, cell_content = self._make_cell_file(tmp_path)

        import os

        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            with patch("ase.optimize.FIRE") as mock_fire_cls:
                mock_dyn = MagicMock()
                mock_dyn.run.return_value = False
                mock_dyn.nsteps = 7
                mock_fire_cls.return_value = mock_dyn

                runner = AirssMlRelaxRunner("fake:Calc@medium", max_steps=7)
                rc = runner.run(name, cell_content)

                assert rc == 0
                assert atoms.info["relax_converged"] is False
                assert atoms.info["relax_status"] == "max_steps"
                assert atoms.info["relax_steps"] == 7
        finally:
            os.chdir(orig)


def test_torchsim_relax_converts_pressure_gpa_to_ev_ang3(monkeypatch, tmp_path):
    """TorchSim batch relax keeps CLI pressure in GPa at the wrapper boundary."""
    from airsspy.jf import ml_runners
    from airsspy.jf.ml_runners import EV_PER_ANG3_TO_GPA

    captured = {}

    class FakeDevice:
        def __init__(self, value):
            self.type = value

        def __str__(self):
            return self.type

    fake_torch = SimpleNamespace(
        device=FakeDevice,
        cuda=SimpleNamespace(is_available=lambda: False),
        float32="float32",
        float64="float64",
    )

    fake_ts = SimpleNamespace(
        Optimizer=SimpleNamespace(
            fire="fire",
            lbfgs="lbfgs",
            bfgs="bfgs",
            gradient_descent="gradient_descent",
        ),
        CellFilter=SimpleNamespace(frechet="frechet", unit="unit"),
        io=SimpleNamespace(
            atoms_to_state=lambda atoms, device, dtype: {"atoms": atoms},
            state_to_atoms=lambda state: state["atoms"],
        ),
        optimize=lambda **kwargs: captured.setdefault("kwargs", kwargs)["system"],
    )

    monkeypatch.setitem(sys.modules, "torch", fake_torch)
    monkeypatch.setitem(sys.modules, "torch_sim", fake_ts)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(ml_runners, "_load_torchsim_model", lambda *a, **k: MagicMock())
    monkeypatch.setattr(
        ml_runners,
        "_cell_to_atoms",
        lambda path: Atoms("Si", positions=[[0, 0, 0]], cell=[3, 3, 3], pbc=True),
    )
    written = {}
    monkeypatch.setattr(
        ml_runners,
        "ase_write",
        lambda filename, atoms_obj, format: written.setdefault(filename, atoms_obj),
    )

    rc = ml_runners._torchsim_relax_batch(
        "mace:medium-mpa-0",
        ["Si-001"],
        ["cell"],
        max_steps=1,
        scalar_pressure=10.0,
    )

    assert rc == {"Si-001": 0}
    assert captured["kwargs"]["init_kwargs"]["scalar_pressure"] == pytest.approx(
        10.0 / EV_PER_ANG3_TO_GPA
    )
    assert written["Si-001.extxyz"].info["extern_pressure"] == pytest.approx(10.0)


def test_torchsim_relax_atoms_input_skips_cell_parser(monkeypatch, tmp_path):
    """TorchSim batch relax accepts Atoms without temporary .cell parsing."""
    from airsspy.jf import ml_runners

    class FakeDevice:
        def __init__(self, value):
            self.type = value

    fake_torch = SimpleNamespace(
        device=FakeDevice,
        cuda=SimpleNamespace(is_available=lambda: False),
        float32="float32",
        float64="float64",
    )
    fake_ts = SimpleNamespace(
        Optimizer=SimpleNamespace(
            fire="fire",
            lbfgs="lbfgs",
            bfgs="bfgs",
            gradient_descent="gradient_descent",
        ),
        CellFilter=SimpleNamespace(frechet="frechet", unit="unit"),
        io=SimpleNamespace(
            atoms_to_state=lambda atoms, device, dtype: {"atoms": atoms},
            state_to_atoms=lambda state: state["atoms"],
        ),
        optimize=lambda **kwargs: kwargs["system"],
    )

    monkeypatch.setitem(sys.modules, "torch", fake_torch)
    monkeypatch.setitem(sys.modules, "torch_sim", fake_ts)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(ml_runners, "_load_torchsim_model", lambda *a, **k: MagicMock())
    monkeypatch.setattr(
        ml_runners,
        "_cell_content_to_atoms",
        lambda content: (_ for _ in ()).throw(AssertionError("cell parser called")),
    )
    monkeypatch.setattr(ml_runners, "ase_write", lambda *a, **k: None)

    atoms = Atoms("Si", positions=[[0, 0, 0]], cell=[3, 3, 3], pbc=True)
    rc = ml_runners._torchsim_relax_batch(
        "mace:medium-mpa-0",
        ["Si-001"],
        [atoms],
        max_steps=1,
    )

    assert rc == {"Si-001": 0}


def test_torchsim_static_batch_preserves_per_structure_forces(monkeypatch, tmp_path):
    """TorchSim static returns per-system forces, not one concatenated force array."""
    from airsspy.jf import ml_runners

    class FakeDevice:
        def __init__(self, value):
            self.type = value

    fake_torch = SimpleNamespace(
        device=FakeDevice,
        cuda=SimpleNamespace(is_available=lambda: False),
        float32="float32",
        float64="float64",
    )
    atoms = [
        Atoms("Si", positions=[[0, 0, 0]], cell=[3, 3, 3], pbc=True),
        Atoms("Si2", positions=[[0, 0, 0], [1, 1, 1]], cell=[4, 4, 4], pbc=True),
    ]
    props_list = [
        {
            "potential_energy": np.array([1.0]),
            "forces": np.array([[[1.0, 2.0, 3.0]]]),
            "stress": np.array([np.eye(3)]),
        },
        {
            "potential_energy": np.array([2.0]),
            "forces": np.array([[[4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]]),
            "stress": np.array([np.eye(3) * 2.0]),
        },
    ]
    written = {}
    fake_ts = SimpleNamespace(
        io=SimpleNamespace(
            atoms_to_state=lambda atoms_in, device, dtype: {"atoms": atoms_in},
            state_to_atoms=lambda state: state["atoms"],
        ),
        static=lambda system, model: props_list,
    )

    monkeypatch.setitem(sys.modules, "torch", fake_torch)
    monkeypatch.setitem(sys.modules, "torch_sim", fake_ts)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(ml_runners, "_load_torchsim_model", lambda *a, **k: MagicMock())
    monkeypatch.setattr(
        ml_runners,
        "ase_write",
        lambda filename, atoms_obj, format: written.setdefault(filename, atoms_obj),
    )

    rc = ml_runners._torchsim_static_batch(
        "mace:medium-mpa-0",
        ["Si-001", "Si2-001"],
        atoms,
        scalar_pressure=10.0,
    )

    assert rc == {"Si-001": 0, "Si2-001": 0}
    assert np.allclose(
        written["Si2-001.extxyz"].calc.results["forces"],
        [[4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
    )
    assert np.allclose(
        written["Si2-001.extxyz"].calc.results["stress"],
        np.eye(3) * 2.0,
    )
    assert written["Si-001.extxyz"].info["extern_pressure"] == pytest.approx(10.0)
    assert written["Si2-001.extxyz"].info["extern_pressure"] == pytest.approx(10.0)


# ---------------------------------------------------------------------------
# compose_ml_task_doc tests
# ---------------------------------------------------------------------------


class TestComposeMlTaskDoc:
    """Tests for compose_ml_task_doc()."""

    def test_compose_from_extxyz(self, tmp_path):
        from airsspy.jf.ml_runners import compose_ml_task_doc

        # Create an Atoms with SinglePointCalculator
        atoms = Atoms(
            "Si2",
            positions=[[0, 0, 0], [1.25, 1.25, 1.25]],
            cell=[[5, 0, 0], [0, 5, 0], [0, 0, 5]],
            pbc=True,
        )
        energy = -10.5
        forces = np.array([[0.01, 0.02, 0.03], [-0.01, -0.02, -0.03]])
        stress = np.array([-0.001, -0.001, -0.001, 0.0, 0.0, 0.0])
        atoms.calc = SinglePointCalculator(
            atoms, energy=energy, forces=forces, stress=stress
        )

        import os

        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            from ase.io import write as ase_write

            name = "test-ml-001"
            ase_write(name + ".extxyz", atoms, format="extxyz")

            result = compose_ml_task_doc(name, calculator_spec="fake:Calc@medium")

            assert result["energy"] == -10.5
            assert result["energy_per_atom"] == pytest.approx(-5.25)
            assert result["natoms"] == 2
            assert result["label"] == name
            assert result["volume"] == pytest.approx(125.0)
            assert "Si" in result["formula"]
            assert (tmp_path / (name + ".res")).is_file()
        finally:
            os.chdir(orig)

    @patch("airsspy.restools.save_airss_res")
    def test_compose_prefers_external_pressure_for_enthalpy(self, mock_save, tmp_path):
        from airsspy.jf.ml_runners import EV_PER_ANG3_TO_GPA, compose_ml_task_doc

        atoms = Atoms(
            "Si2",
            positions=[[0, 0, 0], [1.25, 1.25, 1.25]],
            cell=[[5, 0, 0], [0, 5, 0], [0, 0, 5]],
            pbc=True,
        )
        atoms.info["extern_pressure"] = 10.0
        atoms.calc = SinglePointCalculator(
            atoms,
            energy=-10.5,
            forces=np.zeros((2, 3)),
            stress=np.zeros(6),
        )

        import os

        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            from ase.io import write as ase_write

            name = "test-ml-pressure"
            ase_write(name + ".extxyz", atoms, format="extxyz")

            compose_ml_task_doc(name, calculator_spec="fake:Calc@medium")

            info = mock_save.call_args.args[1]
            assert info["P"] == pytest.approx(10.0)
            assert info["H"] == pytest.approx(-10.5 + 10.0 * 125.0 / EV_PER_ANG3_TO_GPA)
        finally:
            os.chdir(orig)

    def test_compose_includes_relax_metadata(self, tmp_path):
        from airsspy.jf.ml_runners import compose_ml_task_doc

        atoms = Atoms(
            "Si2",
            positions=[[0, 0, 0], [1.25, 1.25, 1.25]],
            cell=[[5, 0, 0], [0, 5, 0], [0, 0, 5]],
            pbc=True,
        )
        atoms.info["relax_converged"] = False
        atoms.info["relax_status"] = "max_steps"
        atoms.info["relax_steps"] = 12
        atoms.calc = SinglePointCalculator(
            atoms,
            energy=-10.5,
            forces=np.zeros((2, 3)),
            stress=np.zeros(6),
        )

        import os

        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            from ase.io import write as ase_write

            name = "test-ml-meta"
            ase_write(name + ".extxyz", atoms, format="extxyz")

            result = compose_ml_task_doc(name, calculator_spec="fake:Calc@medium")

            assert result["relax_converged"] is False
            assert result["relax_status"] == "max_steps"
            assert result["relax_steps"] == 12
            assert "ML Relax status max_steps" in "\n".join(result["rem_lines"])
        finally:
            os.chdir(orig)


# ---------------------------------------------------------------------------
# AirssAbacusSinglePointRunner tests
# ---------------------------------------------------------------------------


class TestAbacusSinglePointRunner:
    """Tests for AirssAbacusSinglePointRunner."""

    def test_prepare_inputs_sets_scf(self, tmp_path):
        from airsspy.jf.runners import AirssAbacusSinglePointRunner

        cell_content = (
            "%BLOCK LATTICE_CART\n"
            "  5.0  0.0  0.0\n"
            "  0.0  5.0  0.0\n"
            "  0.0  0.0  5.0\n"
            "%ENDBLOCK LATTICE_CART\n"
            "%BLOCK POSITIONS_FRAC\n"
            "Si  0.0  0.0  0.0\n"
            "%ENDBLOCK POSITIONS_FRAC\n"
        )
        input_content = "calculation cell-relax\necutwfc 50\n"

        import os

        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            name = "test-abacus-sp-001"
            runner = AirssAbacusSinglePointRunner()
            runner.prepare_inputs(name, cell_content, input_content)

            # Check that INPUT file has calculation scf
            input_path = tmp_path / (name + ".INPUT")
            assert input_path.is_file()
            content = input_path.read_text()
            assert "calculation scf" in content
            assert "cell-relax" not in content
            assert "ecutwfc 50" in content
            assert "press1 0.0" in content
            assert "press2 0.0" in content
            assert "press3 0.0" in content

            # Check that ABACUS workdir and STRU exist
            assert (tmp_path / (name + ".abacus")).is_dir()
            assert (tmp_path / (name + ".abacus") / "STRU").is_file()
            assert (tmp_path / (name + ".abacus") / "INPUT").is_file()
        finally:
            os.chdir(orig)

    def test_prepare_inputs_adds_scf_when_missing(self, tmp_path):
        from airsspy.jf.runners import AirssAbacusSinglePointRunner

        cell_content = (
            "%BLOCK LATTICE_CART\n"
            "  5.0  0.0  0.0\n"
            "  0.0  5.0  0.0\n"
            "  0.0  0.0  5.0\n"
            "%ENDBLOCK LATTICE_CART\n"
            "%BLOCK POSITIONS_FRAC\n"
            "Si  0.0  0.0  0.0\n"
            "%ENDBLOCK POSITIONS_FRAC\n"
        )
        input_content = "ecutwfc 50\n"

        import os

        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            name = "test-abacus-sp-002"
            runner = AirssAbacusSinglePointRunner()
            runner.prepare_inputs(name, cell_content, input_content)

            input_path = tmp_path / (name + ".INPUT")
            content = input_path.read_text()
            assert "calculation scf" in content
        finally:
            os.chdir(orig)

    def test_prepare_inputs_writes_external_pressure_in_kbar(self, tmp_path):
        from airsspy.jf.runners import AirssAbacusSinglePointRunner

        cell_content = (
            "%BLOCK LATTICE_CART\n"
            "  5.0  0.0  0.0\n"
            "  0.0  5.0  0.0\n"
            "  0.0  0.0  5.0\n"
            "%ENDBLOCK LATTICE_CART\n"
            "%BLOCK POSITIONS_FRAC\n"
            "Si  0.0  0.0  0.0\n"
            "%ENDBLOCK POSITIONS_FRAC\n"
        )
        input_content = "calculation cell-relax\necutwfc 50\npress1 1.0\n"

        import os

        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            name = "test-abacus-sp-pressure"
            runner = AirssAbacusSinglePointRunner(pressure=5.0)
            runner.prepare_inputs(name, cell_content, input_content)

            content = (tmp_path / (name + ".INPUT")).read_text()
            assert "calculation scf" in content
            assert "press1 50.0" in content
            assert "press2 50.0" in content
            assert "press3 50.0" in content
            assert "press1 1.0" not in content
        finally:
            os.chdir(orig)


# ---------------------------------------------------------------------------
# CLI integration tests
# ---------------------------------------------------------------------------


class TestCliSp:
    """Tests for the `ap run sp` CLI command."""

    def test_sp_requires_calculator_for_ml(self):
        from click.testing import CliRunner

        from airsspy.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["run", "sp", "--cell", "*.cell", "--code", "ml"])
        assert result.exit_code != 0
        assert "--calculator is required" in result.output

    def test_sp_invalid_code(self):
        from click.testing import CliRunner

        from airsspy.cli.main import cli

        runner = CliRunner()
        # GULP should not be a valid choice for sp
        result = runner.invoke(cli, ["run", "sp", "--cell", "*.cell", "--code", "gulp"])
        assert result.exit_code != 0


class TestCliRelaxMl:
    """Tests for `ap run relax --code ml` CLI options."""

    def test_relax_requires_calculator_for_ml(self):
        from click.testing import CliRunner

        from airsspy.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(
            cli, ["run", "relax", "--cell", "*.cell", "--code", "ml"]
        )
        assert result.exit_code != 0
        assert "--calculator is required" in result.output

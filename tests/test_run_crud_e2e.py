"""End-to-end tests for ``ap run crud`` with real backends.

These tests are skipped by default through the existing ``--run-e2e`` gate.
External executables are configured with environment variables so the harness
can run on developer workstations or CI runners without hardcoded local paths.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from airsspy.cli.main import cli

CASTEP_EXE_ENV = "AIRSSPY_E2E_CASTEP_EXE"
ABACUS_EXE_ENV = "AIRSSPY_E2E_ABACUS_EXE"
ABACUS_PSEUDO_ENV = "AIRSSPY_E2E_ABACUS_PSEUDO"
VASP_EXE_ENV = "AIRSSPY_E2E_VASP_EXE"
VASP_POTCAR_DIR_ENV = "AIRSSPY_E2E_VASP_POTCAR_DIR"
MPINP_ENV = "AIRSSPY_E2E_MPINP"
ML_CALCULATOR_ENV = "AIRSSPY_E2E_ML_CALCULATOR"


def _resolve_executable(env_name: str) -> str:
    """Return an executable from *env_name* or skip the test."""
    value = os.environ.get(env_name)
    if not value:
        pytest.skip(f"{env_name} is not set")

    exe = Path(value).expanduser()
    if exe.is_absolute() or os.sep in value:
        if not exe.is_file():
            pytest.skip(f"{env_name} does not point to a file: {value}")
        if not os.access(exe, os.X_OK):
            pytest.skip(f"{env_name} is not executable: {value}")
        return str(exe)

    if shutil.which(value) is None:
        pytest.skip(f"{env_name} executable is not on PATH: {value}")
    return value


def _resolve_vasp_executable() -> str:
    """Return the configured VASP executable, defaulting to ``vasp_std``."""
    value = os.environ.get(VASP_EXE_ENV, "vasp_std")
    exe = Path(value).expanduser()
    if exe.is_absolute() or os.sep in value:
        if not exe.is_file():
            pytest.skip(f"{VASP_EXE_ENV} does not point to a file: {value}")
        if not os.access(exe, os.X_OK):
            pytest.skip(f"{VASP_EXE_ENV} is not executable: {value}")
        return str(exe)
    if shutil.which(value) is None:
        pytest.skip(f"VASP executable is not on PATH: {value}")
    return value


def _resolve_vasp_potcar_dir() -> str:
    """Return a POTCAR root without requiring POTCAR files in the repo."""
    from airsspy.vasptools import _potcar_root

    root = _potcar_root(os.environ.get(VASP_POTCAR_DIR_ENV))
    if root is not None:
        path = root.expanduser()
        if path.is_dir():
            return str(path)
        pytest.skip(f"Configured VASP POTCAR root does not point to a directory: {root}")
    pytest.skip(
        f"{VASP_POTCAR_DIR_ENV}, AIRSSPY_POTCAR_DIR, PMG_VASP_PSP_DIR, "
        "VASP_PSP_DIR, or pymatgen .pmgrc must be set for VASP e2e tests"
    )


def _mpinp_args() -> list[str]:
    """Return optional CRUD MPI arguments from the e2e environment."""
    value = os.environ.get(MPINP_ENV)
    if value is None:
        return []
    try:
        int(value)
    except ValueError:
        pytest.skip(f"{MPINP_ENV} must be an integer, got {value!r}")
    return ["--mpinp", value]


def _write_hopper_res(root: Path, label: str, element: str) -> None:
    """Create a single-structure RES job in ``hopper/``."""
    root.joinpath("hopper").mkdir(exist_ok=True)
    root.joinpath("hopper", f"{label}.res").write_text(
        f"TITL {label} 0.000 125.000 -1.0000 0.00 0.00 1 (P1) n - 1\n"
        "CELL 1.0 5.000000 5.000000 5.000000 90.000000 90.000000 90.000000\n"
        "LATT -1\n"
        f"SFAC {element}\n"
        f"{element} 1 0.0000000000000 0.0000000000000 0.0000000000000 1.0\n"
        "END\n"
    )


def _write_root_cell(
    root: Path,
    element: str,
    species_pot: str | None = None,
) -> None:
    """Create the root cell file used to reconstruct claimed RES jobs."""
    lines = [
        "%BLOCK LATTICE_CART",
        "5.0 0.0 0.0",
        "0.0 5.0 0.0",
        "0.0 0.0 5.0",
        "%ENDBLOCK LATTICE_CART",
        "",
        "%BLOCK POSITIONS_FRAC",
        f"{element} 0.0 0.0 0.0",
        "%ENDBLOCK POSITIONS_FRAC",
        "",
    ]
    if species_pot is not None:
        lines.extend(
            [
                "%BLOCK SPECIES_POT",
                f"{element} {species_pot}",
                "%ENDBLOCK SPECIES_POT",
                "",
            ]
        )
    root.joinpath(f"{element}.cell").write_text("\n".join(lines))


def _assert_crud_success(root: Path, label: str) -> str:
    """Assert common CRUD success artifacts and return output RES text."""
    assert not root.joinpath("hopper", f"{label}.res").exists()

    res_path = root / "good_castep" / f"{label}.res"
    assert res_path.is_file()
    res_text = res_path.read_text()
    assert f"TITL {label}" in res_text
    assert "REM" in res_text
    return res_text


def _write_vasp_inputs(root: Path, input_set: str, *, nsw: int | None = None) -> None:
    """Create small VASP INCAR/KPOINTS inputs for Si e2e smoke tests."""
    incar_lines = [
        f"AIRSSPY_VASP_INPUT_SET = {input_set}",
        "ENCUT = 200",
        "EDIFF = 1E-4",
        "NELM = 20",
        "ISMEAR = 0",
        "SIGMA = 0.05",
        "LWAVE = .FALSE.",
        "LCHARG = .FALSE.",
    ]
    if nsw is not None:
        incar_lines.append(f"NSW = {nsw}")
    root.joinpath("Si.INCAR").write_text("\n".join(incar_lines) + "\n")
    root.joinpath("Si.KPOINTS").write_text(
        "Automatic mesh\n"
        "0\n"
        "Gamma\n"
        "1 1 1\n"
        "0 0 0\n"
    )


@pytest.fixture(autouse=True)
def _e2e_single_thread(monkeypatch):
    """Keep tiny real-executable smoke tests from oversubscribing CPUs."""
    monkeypatch.setenv("OMP_NUM_THREADS", "1")


@pytest.fixture
def crud_workdir(tmp_path, monkeypatch):
    """Use one pytest tmp directory as both cwd and CRUD workdir."""
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.mark.e2e
def test_run_crud_castep_real_executable_smoke(crud_workdir):
    """CASTEP CRUD claims one RES job, runs CASTEP SP, and collects it."""
    exe = _resolve_executable(CASTEP_EXE_ENV)
    pseudo_src = Path(__file__).parent / "data" / "Si_QC5_PBE_OTF.usp"
    shutil.copy2(pseudo_src, crud_workdir / pseudo_src.name)

    _write_root_cell(crud_workdir, "Si", pseudo_src.name)
    _write_hopper_res(crud_workdir, "Si-001", "Si")
    crud_workdir.joinpath("Si.param").write_text(
        "task : geometryoptimization\n"
        "xc_functional : PBE\n"
        "cut_off_energy : 120 eV\n"
        "max_scf_cycles : 30\n"
        "geom_max_iter : 2\n"
        "finite_basis_corr : 0\n"
        "fixed_npw : true\n"
        "opt_strategy : speed\n"
        "write_checkpoint : none\n"
        "write_otfg : false\n"
        "write_bib : false\n"
        "geom_method : LBFGS\n"
    )

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "crud",
            "--code",
            "castep",
            "--exe",
            exe,
            "--max-iterations",
            "6",
            "--singlepoint",
            "--keep",
            *_mpinp_args(),
        ],
    )

    assert result.exit_code == 0, result.output
    res_text = _assert_crud_success(crud_workdir, "Si-001")
    assert "PBE" in res_text
    assert crud_workdir.joinpath("good_castep", "Si-001.castep").is_file()


@pytest.mark.e2e
def test_run_crud_abacus_real_executable_smoke(crud_workdir):
    """ABACUS CRUD claims one RES job, runs ABACUS SP, and collects it."""
    exe = _resolve_executable(ABACUS_EXE_ENV)
    pseudo_value = os.environ.get(ABACUS_PSEUDO_ENV)
    if not pseudo_value:
        pytest.skip(f"{ABACUS_PSEUDO_ENV} is not set")
    pseudo_src = Path(pseudo_value).expanduser()
    if not pseudo_src.is_file():
        pytest.skip(f"{ABACUS_PSEUDO_ENV} does not point to a file: {pseudo_value}")

    shutil.copy2(pseudo_src, crud_workdir / pseudo_src.name)
    _write_root_cell(crud_workdir, "Si", pseudo_src.name)
    _write_hopper_res(crud_workdir, "Si-001", "Si")
    crud_workdir.joinpath("Si.INPUT").write_text(
        "INPUT_PARAMETERS\n"
        "pseudo_dir ../\n"
        "basis_type pw\n"
        "ecutwfc 20\n"
        "scf_thr 1e-5\n"
        "scf_nmax 50\n"
        "calculation relax\n"
        "kspacing 1.0\n"
        "relax_nmax 2\n"
        "cal_stress 1\n"
    )

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "crud",
            "--code",
            "abacus",
            "--exe",
            exe,
            "--max-iterations",
            "8",
            "--singlepoint",
            "--keep",
            *_mpinp_args(),
        ],
    )

    assert result.exit_code == 0, result.output
    res_text = _assert_crud_success(crud_workdir, "Si-001")
    assert "Basis pw" in res_text
    assert crud_workdir.joinpath("good_castep", "Si-001.abacus").is_dir()


@pytest.mark.e2e
def test_run_crud_vasp_real_executable_smoke(crud_workdir):
    """VASP CRUD claims one RES job, runs VASP SP, and collects it."""
    exe = _resolve_vasp_executable()
    potcar_dir = _resolve_vasp_potcar_dir()

    _write_root_cell(crud_workdir, "Si")
    _write_hopper_res(crud_workdir, "Si-001", "Si")
    _write_vasp_inputs(crud_workdir, "MPStaticSet")

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "crud",
            "--code",
            "vasp",
            "--exe",
            exe,
            "--potcar-dir",
            potcar_dir,
            "--potcar-map",
            "Si=Si",
            "--max-iterations",
            "1",
            "--singlepoint",
            "--keep",
            *_mpinp_args(),
        ],
    )

    assert result.exit_code == 0, result.output
    res_text = _assert_crud_success(crud_workdir, "Si-001")
    assert "VASP input set MPStaticSet" in res_text
    assert "POTCAR Si Si sha256=" in res_text
    assert crud_workdir.joinpath("good_castep", "Si-001.vasp", "OUTCAR").is_file()


@pytest.mark.e2e
def test_run_crud_vasp_relax_real_executable_smoke(crud_workdir):
    """VASP CRUD claims one RES job, runs a short relaxation, and collects it."""
    exe = _resolve_vasp_executable()
    potcar_dir = _resolve_vasp_potcar_dir()

    _write_root_cell(crud_workdir, "Si")
    _write_hopper_res(crud_workdir, "Si-001", "Si")
    _write_vasp_inputs(crud_workdir, "MPRelaxSet", nsw=1)

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "crud",
            "--code",
            "vasp",
            "--exe",
            exe,
            "--potcar-dir",
            potcar_dir,
            "--potcar-map",
            "Si=Si",
            "--max-iterations",
            "1",
            "--keep",
            *_mpinp_args(),
        ],
    )

    assert result.exit_code == 0, result.output
    res_text = _assert_crud_success(crud_workdir, "Si-001")
    assert "VASP input set MPRelaxSet" in res_text
    incar_text = crud_workdir.joinpath("good_castep", "Si-001.vasp", "INCAR").read_text()
    assert "NSW = 1" in incar_text
    assert crud_workdir.joinpath("good_castep", "Si-001.vasp", "OUTCAR").is_file()


@pytest.mark.e2e
def test_run_relax_vasp_real_executable_smoke(crud_workdir):
    """VASP run-relax accepts RES input and runs a short relaxation."""
    exe = _resolve_vasp_executable()
    potcar_dir = _resolve_vasp_potcar_dir()

    _write_root_cell(crud_workdir, "Si")
    _write_hopper_res(crud_workdir, "Si-001", "Si")
    shutil.move(
        crud_workdir / "hopper" / "Si-001.res",
        crud_workdir / "Si-001.res",
    )
    _write_vasp_inputs(crud_workdir, "MPRelaxSet", nsw=1)

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "relax",
            "--cell",
            "Si-001.res",
            "--code",
            "vasp",
            "--exe",
            exe,
            "--potcar-dir",
            potcar_dir,
            "--potcar-map",
            "Si=Si",
            "--max-iterations",
            "1",
            "--keep",
            *_mpinp_args(),
        ],
    )

    assert result.exit_code == 0, result.output
    res_text = crud_workdir.joinpath("Si-001.res").read_text()
    assert "VASP input set MPRelaxSet" in res_text
    incar_text = crud_workdir.joinpath("Si-001.vasp", "INCAR").read_text()
    assert "NSW = 1" in incar_text
    assert crud_workdir.joinpath("Si-001.vasp", "OUTCAR").is_file()


@pytest.mark.e2e
def test_run_crud_ml_ase_calculator_smoke(crud_workdir):
    """ML CRUD uses a real ASE calculator implementation and collects output."""
    calculator_spec = os.environ.get(
        ML_CALCULATOR_ENV,
        "ase.calculators.emt:EMT",
    )
    element = "Si" if ML_CALCULATOR_ENV in os.environ else "Al"
    label = f"{element}-001"

    _write_root_cell(crud_workdir, element)
    _write_hopper_res(crud_workdir, label, element)

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "crud",
            "--code",
            "ml",
            "--calculator",
            f"ase:{calculator_spec}",
            "--keep",
        ],
    )

    assert result.exit_code == 0, result.output
    res_text = _assert_crud_success(crud_workdir, label)
    assert f"ML Calculator ase:{calculator_spec}" in res_text
    assert crud_workdir.joinpath("good_castep", f"{label}.extxyz").is_file()

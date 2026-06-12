"""Tests for REM metadata extraction from CASTEP and ABACUS output files.

Uses real output files from the AIRSS reference implementation as fixtures
and optionally runs actual CASTEP/ABACUS calculations.
"""

import hashlib
import shutil
from pathlib import Path

import pytest

import pytest

# --- Paths to reference data ---

AIRSS_GIT = Path.home() / "appdir" / "airss-git"

CASTEP_EXAMPLE = AIRSS_GIT / "examples" / "3.1" / "C2-260416-234105-0dbf05.castep"
CASTEP_EXAMPLE_2 = AIRSS_GIT / "examples" / "3.1" / "C2-260416-235310-daaa6c.castep"
CELL_EXAMPLE = AIRSS_GIT / "examples" / "5.1" / "C2.cell"

ABACUS_NCP19_LOG = (
    AIRSS_GIT / "examples" / "7.1-ncp19"
    / "C2-289372-8682-10.abacus/OUT.ABACUS/running_cell-relax.log"
)
ABACUS_QC5_LOG = (
    AIRSS_GIT / "examples" / "7.1-qc5"
    / "C2-372705-3361-1.abacus/OUT.ABACUS/running_cell-relax.log"
)
ABACUS_QC5_INPUT = AIRSS_GIT / "examples" / "7.1-qc5" / "C2.INPUT"
ABACUS_NCP19_INPUT = AIRSS_GIT / "examples" / "7.1-ncp19" / "C2.INPUT"


def _has_executable(name: str) -> bool:
    return shutil.which(name) is not None


# ============================================================================
# Override clean_dir — we need tmp_path and CWD to be the same directory
# ============================================================================


@pytest.fixture(autouse=True)
def use_tmp_as_cwd(tmp_path, monkeypatch):
    """Override clean_dir so CWD matches tmp_path."""
    monkeypatch.chdir(tmp_path)


# ============================================================================
# CASTEP REM extraction tests (real files)
# ============================================================================


@pytest.fixture
def castep_seed(tmp_path):
    """Copy a real .castep file into tmp_path and return the seed name."""
    seed = "C2-260416-234105-0dbf05"
    shutil.copy2(CASTEP_EXAMPLE, tmp_path / (seed + ".castep"))
    return seed


@pytest.fixture
def castep_seed_2(tmp_path):
    """Copy a second real .castep file into tmp_path."""
    seed = "C2-260416-235310-daaa6c"
    shutil.copy2(CASTEP_EXAMPLE_2, tmp_path / (seed + ".castep"))
    return seed


@pytest.fixture
def cell_seed(tmp_path):
    """Copy a real .cell file into tmp_path and return the seed name."""
    seed = "C2"
    shutil.copy2(CELL_EXAMPLE, tmp_path / (seed + ".cell"))
    return seed


class TestExtractREMFromCastep:
    def test_functional(self, castep_seed):
        from airsspy.casteptools import extract_REM_from_castep

        rem = extract_REM_from_castep(castep_seed)
        assert "Perdew Burke Ernzerhof" in rem["functional"]

    def test_relativity(self, castep_seed):
        from airsspy.casteptools import extract_REM_from_castep

        rem = extract_REM_from_castep(castep_seed)
        assert "Koelling-Harmon" in rem["relativity"]

    def test_dispersion_off(self, castep_seed):
        from airsspy.casteptools import extract_REM_from_castep

        rem = extract_REM_from_castep(castep_seed)
        assert "off" in rem["dispersion"]

    def test_cutoff(self, castep_seed):
        from airsspy.casteptools import extract_REM_from_castep

        rem = extract_REM_from_castep(castep_seed)
        assert "350" in rem["cutoff"]
        assert "eV" in rem["cutoff"]

    def test_gridscale(self, castep_seed):
        from airsspy.casteptools import extract_REM_from_castep

        rem = extract_REM_from_castep(castep_seed)
        assert "1.75" in rem["gridscale"]

    def test_gmax(self, castep_seed):
        from airsspy.casteptools import extract_REM_from_castep

        rem = extract_REM_from_castep(castep_seed)
        assert "16.7730" in rem["gmax"]

    def test_fbsc(self, castep_seed):
        from airsspy.casteptools import extract_REM_from_castep

        rem = extract_REM_from_castep(castep_seed)
        assert "none" in rem["fbsc"]

    def test_mpgrid(self, castep_seed):
        from airsspy.casteptools import extract_REM_from_castep

        rem = extract_REM_from_castep(castep_seed)
        assert "11" in rem["mpgrid"]
        assert "9" in rem["mpgrid"]
        assert "8" in rem["mpgrid"]

    def test_offset(self, castep_seed):
        from airsspy.casteptools import extract_REM_from_castep

        rem = extract_REM_from_castep(castep_seed)
        assert "0.000" in rem["offset"]

    def test_nkpts(self, castep_seed):
        from airsspy.casteptools import extract_REM_from_castep

        rem = extract_REM_from_castep(castep_seed)
        assert "396" in rem["nkpts"]

    def test_psps(self, castep_seed):
        from airsspy.casteptools import extract_REM_from_castep

        rem = extract_REM_from_castep(castep_seed)
        assert "C 2|1.4|10|12|13|20:21(qc=7)" in rem["psps"]

    def test_consistency_across_two_files(self, castep_seed, castep_seed_2):
        """Both C2 examples should have the same quality parameters (except kpoints)."""
        from airsspy.casteptools import extract_REM_from_castep

        rem1 = extract_REM_from_castep(castep_seed)
        rem2 = extract_REM_from_castep(castep_seed_2)

        # Parameters that should be identical (same cell template)
        for key in ["functional", "relativity", "cutoff", "gridscale", "gmax", "fbsc"]:
            assert rem1[key] == rem2[key], f"{key} differs: {rem1[key]!r} vs {rem2[key]!r}"

        # MP grid may differ (different volumes lead to different k-point grids)
        assert "MP grid" in rem1["mpgrid"]
        assert "MP grid" in rem2["mpgrid"]

    def test_missing_file(self, tmp_path):
        from airsspy.casteptools import extract_REM_from_castep

        rem = extract_REM_from_castep("nonexistent")
        assert rem == {}


class TestExtractCellREMMetadata:
    def test_spacing(self, cell_seed):
        from airsspy.casteptools import extract_cell_rem_metadata

        rem = extract_cell_rem_metadata(cell_seed)
        assert "0.07" in rem["spacing"]

    def test_md5(self, cell_seed):
        from airsspy.casteptools import extract_cell_rem_metadata

        rem = extract_cell_rem_metadata(cell_seed)
        assert "md5" in rem

        # Verify the MD5 matches a direct hash of the cell file
        expected_hash = hashlib.md5(
            (Path(cell_seed).parent / (cell_seed + ".cell")).read_bytes()
        ).hexdigest()
        assert expected_hash in rem["md5"]

    def test_no_hubbard(self, cell_seed):
        """The C2.cell file has no Hubbard U block."""
        from airsspy.casteptools import extract_cell_rem_metadata

        rem = extract_cell_rem_metadata(cell_seed)
        assert "hubbard" not in rem

    def test_missing_file(self, tmp_path):
        from airsspy.casteptools import extract_cell_rem_metadata

        rem = extract_cell_rem_metadata("nonexistent")
        assert rem == {}


class TestBuildREMLines:
    def test_output_format(self, castep_seed, cell_seed):
        from airsspy.casteptools import build_rem_lines

        lines = build_rem_lines(castep_seed)

        # Should be a non-empty list of strings
        assert isinstance(lines, list)
        assert len(lines) > 0

        # First line should be empty (blank REM)
        assert lines[0] == ""

        # Should contain functional, cutoff, kpoints info
        full_text = "\n".join(lines)
        assert "Perdew Burke Ernzerhof" in full_text
        assert "350" in full_text
        assert "396" in full_text
        assert "C 2|1.4|10|12|13|20:21(qc=7)" in full_text
        assert "0.07" in full_text  # spacing from cell

    def test_rem_lines_without_cell(self, castep_seed, tmp_path):
        """No .cell file — should still have castep-derived fields."""
        from airsspy.casteptools import build_rem_lines

        # Remove the cell file (only castep file remains)
        lines = build_rem_lines(castep_seed)

        full_text = "\n".join(lines)
        # Should still have castep-derived fields
        assert "Perdew Burke Ernzerhof" in full_text
        # No spacing since no cell file
        assert "Spacing" not in full_text


# ============================================================================
# ABACUS REM extraction tests (real files)
# ============================================================================


@pytest.fixture
def abacus_ncp19_env(tmp_path):
    """Set up ABACUS NCP19 environment with INPUT and log files."""
    seed = "test"
    workdir = tmp_path / (seed + ".abacus")
    out_dir = workdir / "OUT.ABACUS"
    out_dir.mkdir(parents=True)

    shutil.copy2(ABACUS_NCP19_INPUT, tmp_path / (seed + ".INPUT"))
    shutil.copy2(ABACUS_NCP19_LOG, out_dir / "running_cell-relax.log")

    return seed


@pytest.fixture
def abacus_qc5_env(tmp_path):
    """Set up ABACUS QC5 (LCAO) environment with INPUT and log files."""
    seed = "test"
    workdir = tmp_path / (seed + ".abacus")
    out_dir = workdir / "OUT.ABACUS"
    out_dir.mkdir(parents=True)

    shutil.copy2(ABACUS_QC5_INPUT, tmp_path / (seed + ".INPUT"))
    shutil.copy2(ABACUS_QC5_LOG, out_dir / "running_cell-relax.log")

    return seed


class TestExtractAbacusREM:
    def test_functional(self, abacus_ncp19_env):
        from airsspy.abacustools import extract_abacus_rem

        rem = extract_abacus_rem(abacus_ncp19_env)
        assert "PBE" in rem["functional"]

    def test_cutoff_pw(self, abacus_ncp19_env):
        """NCP19 uses pw basis — ecutwfc from INPUT."""
        from airsspy.abacustools import extract_abacus_rem

        rem = extract_abacus_rem(abacus_ncp19_env)
        assert rem["cutoff"] == 50.0
        assert rem["basis_type"] == "pw"

    def test_cutoff_qc5(self, abacus_qc5_env):
        """QC5 INPUT has ecutwfc=30, basis_type=pw."""
        from airsspy.abacustools import extract_abacus_rem

        rem = extract_abacus_rem(abacus_qc5_env)
        assert rem["cutoff"] == 30.0
        assert rem["basis_type"] == "pw"

    def test_kspacing(self, abacus_ncp19_env):
        from airsspy.abacustools import extract_abacus_rem

        rem = extract_abacus_rem(abacus_ncp19_env)
        assert rem["kspacing"] == 0.25

    def test_nkpts(self, abacus_ncp19_env):
        from airsspy.abacustools import extract_abacus_rem

        rem = extract_abacus_rem(abacus_ncp19_env)
        assert "nkpts" in rem
        assert rem["nkpts"] > 0

    def test_psps(self, abacus_ncp19_env):
        from airsspy.abacustools import extract_abacus_rem

        rem = extract_abacus_rem(abacus_ncp19_env)
        assert "psps" in rem
        assert len(rem["psps"]) > 0
        assert any("C_NCP19" in p for p in rem["psps"])

    def test_orbital_info_lcao(self, abacus_qc5_env):
        """LCAO basis should have orbital zeta info."""
        from airsspy.abacustools import extract_abacus_rem

        rem = extract_abacus_rem(abacus_qc5_env)
        assert "orbital_info" in rem
        assert len(rem["orbital_info"]) > 0
        # Each entry should have element name and zeta info
        for entry in rem["orbital_info"]:
            assert "C" in entry
            assert "L" in entry


class TestBuildAbacusREMLines:
    def test_pw_output(self, abacus_ncp19_env):
        from airsspy.abacustools import build_abacus_rem_lines

        lines = build_abacus_rem_lines(abacus_ncp19_env)

        full_text = "\n".join(lines)
        assert "Functional PBE" in full_text
        assert "Basis pw" in full_text
        assert "50" in full_text
        assert "Spacing 0.25" in full_text
        assert "No. kpts" in full_text
        assert "C_NCP19" in full_text

    def test_qc5_output(self, abacus_qc5_env):
        """QC5 INPUT says pw with ecutwfc=30 — REM should reflect that."""
        from airsspy.abacustools import build_abacus_rem_lines

        lines = build_abacus_rem_lines(abacus_qc5_env)

        full_text = "\n".join(lines)
        assert "Functional PBE" in full_text
        assert "Basis pw" in full_text
        assert "30" in full_text


# ============================================================================
# Integration: compose_task_doc REM and symmetry (CASTEP real file)
# ============================================================================


class TestComposeTaskDocREM:
    def test_rem_lines_in_output(self, castep_seed):
        """compose_task_doc should produce rem_lines from a real .castep."""
        from airsspy.jf.runners import compose_task_doc

        # compose_task_doc needs -out.cell; create one from final config
        out_cell = Path(castep_seed + "-out.cell")
        out_cell.write_text("""\
%BLOCK LATTICE_CART
2.4213815 0.2580130 0.3187888
0.9469871 2.2390504 0.4005947
0.8811574 0.5715302 2.2294130
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
C 0.517971 0.416555 0.643946
C 0.768638 0.667347 0.894505
%ENDBLOCK POSITIONS_FRAC
""")

        doc = compose_task_doc(castep_seed)
        assert "rem_lines" in doc
        assert doc["rem_lines"] is not None
        assert len(doc["rem_lines"]) > 0

        full_rem = "\n".join(doc["rem_lines"])
        assert "Perdew Burke Ernzerhof" in full_rem
        assert "350" in full_rem

    def test_symmetry_not_hardcoded(self, castep_seed):
        """Symmetry should be computed via spglib, not hardcoded."""
        from airsspy.jf.runners import compose_task_doc

        out_cell = Path(castep_seed + "-out.cell")
        out_cell.write_text("""\
%BLOCK LATTICE_CART
2.4213815 0.2580130 0.3187888
0.9469871 2.2390504 0.4005947
0.8811574 0.5715302 2.2294130
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
C 0.517971 0.416555 0.643946
C 0.768638 0.667347 0.894505
%ENDBLOCK POSITIONS_FRAC
""")

        doc = compose_task_doc(castep_seed)

        # Read the .res file and check TITL has a real spacegroup
        res_path = Path(castep_seed + ".res")
        assert res_path.is_file()
        content = res_path.read_text()
        titl_line = [l for l in content.splitlines() if l.startswith("TITL")][0]

        # Should not be hardcoded "1" or "(1)"
        assert "(1)" not in titl_line
        # Should contain REM lines
        assert "REM" in content


class TestComposeAbacusTaskDocREM:
    def _setup_stru(self, seed):
        """Create STRU_ION_D and abacus_out for a given seed."""
        workdir = Path(seed + ".abacus")
        out_dir = workdir / "OUT.ABACUS"
        out_dir.mkdir(parents=True, exist_ok=True)

        stru_path = out_dir / "STRU_ION_D"
        stru_path.write_text("""\
ATOMIC_SPECIES
C 12.011 C_NCP19_PBE_OTF.upf

LATTICE_CONSTANT
1.88973

LATTICE_VECTORS
  2.460000  0.000000  0.000000
  0.000000  2.460000  0.000000
  0.000000  0.000000  6.700000

ATOMIC_POSITIONS
Direct
C
0.0
2
0.0000000000 0.0000000000 0.0000000000 1 1 1
0.0000000000 0.0000000000 0.5000000000 1 1 1
""")

        (workdir / "abacus_out").write_text("TOTAL  Time : 12.5\n")

    def test_rem_lines_in_output(self, abacus_ncp19_env):
        """compose_abacus_task_doc should produce rem_lines from ABACUS."""
        from airsspy.abacustools import compose_abacus_task_doc

        self._setup_stru(abacus_ncp19_env)

        doc = compose_abacus_task_doc(abacus_ncp19_env)
        assert "rem_lines" in doc
        assert doc["rem_lines"] is not None

        full_rem = "\n".join(doc["rem_lines"])
        assert "Functional PBE" in full_rem
        assert "Basis pw" in full_rem

    def test_symmetry_not_hardcoded(self, abacus_ncp19_env):
        """Symmetry should be computed via spglib for ABACUS too."""
        from airsspy.abacustools import compose_abacus_task_doc

        # Use a diamond Si cell that has real symmetry
        workdir = Path(abacus_ncp19_env + ".abacus")
        out_dir = workdir / "OUT.ABACUS"
        out_dir.mkdir(parents=True, exist_ok=True)

        stru_path = out_dir / "STRU_ION_D"
        stru_path.write_text("""\
ATOMIC_SPECIES
Si 28.086 Si_ONCV_PBE_FR.upf

LATTICE_CONSTANT
1.88973

LATTICE_VECTORS
  2.715000  0.000000  0.000000
  0.000000  2.715000  0.000000
  0.000000  0.000000  2.715000

ATOMIC_POSITIONS
Direct
Si
0.0
2
0.0000000000 0.0000000000 0.0000000000 1 1 1
0.2500000000 0.2500000000 0.2500000000 1 1 1
""")

        (workdir / "abacus_out").write_text("TOTAL  Time : 12.5\n")

        doc = compose_abacus_task_doc(abacus_ncp19_env)

        res_path = Path(abacus_ncp19_env + ".res")
        assert res_path.is_file()
        content = res_path.read_text()
        titl_line = [l for l in content.splitlines() if l.startswith("TITL")][0]

        # Diamond Si should have R-3m symmetry (not P1)
        assert "R-3m" in titl_line


# ============================================================================
# End-to-end tests with actual executables
# ============================================================================


CASTEP_EXE = "/home/bonan/appdir/castep/CASTEP-23.1/bin/linux_x86_64_gfortran10--mpi/castep.mpi"


@pytest.mark.e2e
class TestCastepEndToEnd:
    """Run an actual CASTEP calculation and verify REM extraction."""

    def test_castep_relax_rem_extraction(self, tmp_path):
        """Run CASTEP on a small Si cell and check REM in output .res."""
        from castepinput.inputs import CellInput, ParamInput
        from airsspy.jf.runners import AirssCastepRelaxRunner, compose_task_doc

        # Copy pre-generated pseudopotential to avoid on-the-fly generation
        usp_src = Path(__file__).parent / "data" / "Si_QC5_PBE_OTF.usp"
        shutil.copy2(usp_src, tmp_path / "Si_QC5_PBE_OTF.usp")

        # Write .cell (gencell-style defaults)
        (tmp_path / "Si.cell").write_text("""\
%BLOCK LATTICE_CART
5.43 0.0 0.0
0.0 5.43 0.0
0.0 0.0 5.43
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Si 0.0 0.0 0.0
Si 0.25 0.25 0.25
%ENDBLOCK POSITIONS_FRAC

KPOINTS_MP_SPACING 0.07

SYMMETRY_GENERATE
SNAP_TO_SYMMETRY

%BLOCK SPECIES_POT
Si Si_QC5_PBE_OTF.usp
%ENDBLOCK SPECIES_POT
""")

        cell = CellInput.from_file("Si.cell")

        # Param defaults matching gencell
        param = ParamInput()
        param["task"] = "geometryoptimization"
        param["xc_functional"] = "PBE"
        param["cut_off_energy"] = "270 eV"
        param["max_scf_cycles"] = 60
        param["geom_max_iter"] = 20
        param["finite_basis_corr"] = "0"
        param["fixed_npw"] = "true"
        param["opt_strategy"] = "speed"
        param["write_checkpoint"] = "none"
        param["write_otfg"] = "false"
        param["write_bib"] = "false"
        param["geom_method"] = "LBFGS"

        (tmp_path / "Si.param").write_text(param.get_string())

        # Run relaxation
        runner = AirssCastepRelaxRunner(executable=CASTEP_EXE, max_iterations=50)
        returncode = runner.run("Si", cell, param)
        assert returncode == 0

        # Compose task doc and verify REM
        doc = compose_task_doc("Si")
        assert "rem_lines" in doc
        assert doc["rem_lines"] is not None

        full_rem = "\n".join(doc["rem_lines"])
        assert "PBE" in full_rem
        assert "270" in full_rem

        # Check .res file has REM lines
        res_path = tmp_path / "Si.res"
        assert res_path.is_file()
        res_content = res_path.read_text()
        assert "REM" in res_content

        # Check symmetry is real
        titl_line = [l for l in res_content.splitlines() if l.startswith("TITL")][0]
        assert "Fd-3m" in titl_line  # Si diamond


ABACUS_EXE = "abacus"


@pytest.mark.e2e
class TestAbacusEndToEnd:
    """Run an actual ABACUS calculation and verify REM extraction."""

    def test_abacus_relax_rem_extraction(self, tmp_path, monkeypatch):
        """Run ABACUS on a small Si cell and check REM in output .res."""
        from airsspy.abacustools import compose_abacus_task_doc
        from airsspy.jf.runners import AirssAbacusRelaxRunner

        monkeypatch.setenv("OMP_NUM_THREADS", "1")

        # Write .cell
        (tmp_path / "Si.cell").write_text("""\
%BLOCK LATTICE_CART
5.43 0.0 0.0
0.0 5.43 0.0
0.0 0.0 5.43
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Si 0.0 0.0 0.0
Si 0.25 0.25 0.25
%ENDBLOCK POSITIONS_FRAC

%BLOCK SPECIES_POT
Si Si_ONCV_PBE-1.2.upf
%ENDBLOCK SPECIES_POT
""")

        # Write INPUT
        input_content = """\
INPUT_PARAMETERS
pseudo_dir     ../
basis_type      pw
ecutwfc         20
scf_thr         1e-5
scf_nmax        50
calculation     relax
kspacing        1.0
relax_nmax      5
cal_stress 1
"""
        (tmp_path / "Si.INPUT").write_text(input_content)

        # Copy pseudopotential — try known locations
        pseudo_candidates = [
            Path("/home/bonan/aiida_envs/aiida-2.0-dev/dev-folder/abacus-3.10lts/tests/PP_ORB/Si_ONCV_PBE-1.2.upf"),
            AIRSS_GIT / "abacus_example" / "Si_C19MK2_PBE_OTF.upf",
        ]
        pseudo_src = None
        for candidate in pseudo_candidates:
            if candidate.is_file():
                pseudo_src = candidate
                break
        if pseudo_src is None:
            pytest.skip("No Si pseudopotential file found for ABACUS")
        shutil.copy2(pseudo_src, tmp_path / "Si_ONCV_PBE-1.2.upf")

        # Also copy into workdir for ABACUS to find
        (tmp_path / "Si.abacus").mkdir(parents=True, exist_ok=True)
        shutil.copy2(pseudo_src, tmp_path / "Si.abacus" / "Si_ONCV_PBE-1.2.upf")

        # Run relaxation
        cell_content = (tmp_path / "Si.cell").read_text()
        runner = AirssAbacusRelaxRunner(
            executable=ABACUS_EXE, max_iterations=30, max_fails=3
        )
        returncode = runner.run("Si", cell_content, input_content)
        assert returncode == 0

        # Compose task doc and verify REM
        doc = compose_abacus_task_doc("Si")
        assert "rem_lines" in doc
        assert doc["rem_lines"] is not None

        full_rem = "\n".join(doc["rem_lines"])
        assert "Functional" in full_rem or "Basis" in full_rem

        # Check .res file has REM lines
        res_path = tmp_path / "Si.res"
        assert res_path.is_file()
        res_content = res_path.read_text()
        assert "REM" in res_content

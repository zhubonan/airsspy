"""Tests for RES ↔ extxyz conversion."""

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

from airsspy.convert import (
    _parse_res_forces,
    extract_structure,
    extxyz_to_res,
    res_to_extxyz,
)
from airsspy.restools import _read_res, parse_titl

# ---------------------------------------------------------------------------
# Sample data
# ---------------------------------------------------------------------------

RES_BASIC = """\
TITL Si-001 -0.05 40.0 -42.500000 0 0 4 (Fd-3m) n - 1
REM Test structure
REM AIRSS Version 0.9.1
CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0
LATT -1
SFAC Si
Si     1  0.0000000000  0.0000000000  0.0000000000 1.0
Si     1  0.2500000000  0.2500000000  0.2500000000 1.0
Si     1  0.5000000000  0.5000000000  0.5000000000 1.0
Si     1  0.7500000000  0.7500000000  0.7500000000 1.0
END
"""

RES_WITH_SPIN = """\
TITL Fe-001 0.0 30.0 -25.000000 4.0 2.0 2 (Im-3m) n - 3
CELL 1.0  2.87 2.87 2.87 90.0 90.0 90.0
LATT -1
SFAC Fe
Fe     1  0.0  0.0  0.0  1.0  2.5
Fe     1  0.5  0.5  0.5  1.0 -2.5
END
"""

RES_WITH_FORCES = """\
TITL Si-002 -0.05 40.0 -43.000000 0 0 4 (Pm-3m) n - 1
CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0
LATT -1
SFAC Si
Si     1  0.0000000000  0.0000000000  0.0000000000 1.0    0.010000   -0.020000    0.030000
Si     1  0.2500000000  0.2500000000  0.2500000000 1.0   -0.010000    0.020000   -0.030000
Si     1  0.5000000000  0.5000000000  0.5000000000 1.0    0.005000   -0.005000    0.005000
Si     1  0.7500000000  0.7500000000  0.7500000000 1.0   -0.005000    0.005000   -0.005000
END
"""

RES_WITH_SPIN_FORCES = """\
TITL Fe-002 0.0 30.0 -26.000000 4.0 2.0 2 (Im-3m) n - 1
CELL 1.0  2.87 2.87 2.87 90.0 90.0 90.0
LATT -1
SFAC Fe
Fe     1  0.0  0.0  0.0  1.0  2.5    0.100000   -0.200000    0.300000
Fe     1  0.5  0.5  0.5  1.0 -2.5   -0.100000    0.200000   -0.300000
END
"""

PACKED_RES = RES_BASIC.strip() + "\n" + RES_WITH_FORCES.strip() + "\n"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _read_single_res_from_dir(output_dir, label):
    """Read the .res file for a given label from the output directory."""
    res_file = output_dir / f"{label}.res"
    assert res_file.exists(), f"Expected {res_file} to exist"
    return res_file.read_text()


# ---------------------------------------------------------------------------
# Force parsing tests
# ---------------------------------------------------------------------------


class TestParseForces:
    def test_no_forces(self):
        lines = RES_BASIC.strip().splitlines()
        assert _parse_res_forces(lines) is None

    def test_with_forces_no_spin(self):
        lines = RES_WITH_FORCES.strip().splitlines()
        forces = _parse_res_forces(lines)
        assert forces is not None
        assert len(forces) == 4
        assert forces[0] == pytest.approx([0.01, -0.02, 0.03])
        assert forces[1] == pytest.approx([-0.01, 0.02, -0.03])
        assert _read_res(lines)["spins"] == []

    def test_with_spin_no_forces(self):
        lines = RES_WITH_SPIN.strip().splitlines()
        assert _parse_res_forces(lines) is None

    def test_with_spin_and_forces(self):
        lines = RES_WITH_SPIN_FORCES.strip().splitlines()
        forces = _parse_res_forces(lines)
        assert forces is not None
        assert len(forces) == 2
        assert forces[0] == pytest.approx([0.1, -0.2, 0.3])


# ---------------------------------------------------------------------------
# RES → extxyz tests
# ---------------------------------------------------------------------------


class TestResToExtxyz:
    def test_basic(self, tmp_path):
        res_file = tmp_path / "test.res"
        xyz_file = tmp_path / "test.xyz"
        res_file.write_text(RES_BASIC)

        n = res_to_extxyz(res_file, xyz_file)
        assert n == 1

        atoms = read(str(xyz_file))
        assert len(atoms) == 4
        assert atoms.get_potential_energy() == pytest.approx(-42.5)
        assert atoms.info["pressure"] == pytest.approx(-0.05)
        assert atoms.info["label"] == "Si-001"
        assert atoms.info["copies"] == 1
        assert atoms.info["symm"] == "(Fd-3m)"

    def test_with_forces(self, tmp_path):
        res_file = tmp_path / "test.res"
        xyz_file = tmp_path / "test.xyz"
        res_file.write_text(RES_WITH_FORCES)

        res_to_extxyz(res_file, xyz_file)

        atoms = read(str(xyz_file))
        forces = atoms.get_forces()
        assert forces.shape == (4, 3)
        assert forces[0] == pytest.approx([0.01, -0.02, 0.03], abs=1e-4)
        assert atoms.get_initial_magnetic_moments().tolist() == pytest.approx(
            [0.0, 0.0, 0.0, 0.0]
        )

    def test_with_spin(self, tmp_path):
        res_file = tmp_path / "test.res"
        xyz_file = tmp_path / "test.xyz"
        res_file.write_text(RES_WITH_SPIN)

        res_to_extxyz(res_file, xyz_file)

        atoms = read(str(xyz_file))
        assert atoms.info["spin"] == pytest.approx(4.0)
        assert atoms.info["spin_abs"] == pytest.approx(2.0)
        assert atoms.info["copies"] == 3
        magmoms = atoms.get_initial_magnetic_moments()
        assert magmoms[0] == pytest.approx(2.5)
        assert magmoms[1] == pytest.approx(-2.5)

    def test_with_spin_and_forces(self, tmp_path):
        res_file = tmp_path / "test.res"
        xyz_file = tmp_path / "test.xyz"
        res_file.write_text(RES_WITH_SPIN_FORCES)

        res_to_extxyz(res_file, xyz_file)

        atoms = read(str(xyz_file))
        forces = atoms.get_forces()
        assert forces[0] == pytest.approx([0.1, -0.2, 0.3], abs=1e-4)
        magmoms = atoms.get_initial_magnetic_moments()
        assert magmoms[0] == pytest.approx(2.5)

    def test_packed(self, tmp_path):
        res_file = tmp_path / "packed.res"
        xyz_file = tmp_path / "packed.xyz"
        res_file.write_text(PACKED_RES)

        n = res_to_extxyz(res_file, xyz_file)
        assert n == 2

        atoms_list = read(str(xyz_file), index=":")
        assert len(atoms_list) == 2
        assert atoms_list[0].info["label"] == "Si-001"
        assert atoms_list[1].info["label"] == "Si-002"

    def test_packed_without_final_end(self, tmp_path):
        res_file = tmp_path / "packed.res"
        xyz_file = tmp_path / "packed.xyz"
        res_file.write_text(RES_BASIC.strip() + "\n" + RES_WITH_FORCES.replace("END\n", ""))

        n = res_to_extxyz(res_file, xyz_file)

        assert n == 2
        atoms_list = read(str(xyz_file), index=":")
        assert [atoms.info["label"] for atoms in atoms_list] == ["Si-001", "Si-002"]

    def test_rem_preserved(self, tmp_path):
        res_file = tmp_path / "test.res"
        xyz_file = tmp_path / "test.xyz"
        res_file.write_text(RES_BASIC)

        res_to_extxyz(res_file, xyz_file)

        atoms = read(str(xyz_file))
        assert "rem" in atoms.info
        rem = atoms.info["rem"]
        assert any("Test structure" in r for r in rem)
        assert any("AIRSS" in r for r in rem)


# ---------------------------------------------------------------------------
# extxyz → RES tests (unpacked: one .res file per structure)
# ---------------------------------------------------------------------------


class TestExtxyzToRes:
    def _make_si_atoms(self, with_forces=False, with_spin=False):
        atoms = Atoms(
            "Si4",
            scaled_positions=[
                [0.0, 0.0, 0.0],
                [0.25, 0.25, 0.25],
                [0.5, 0.5, 0.5],
                [0.75, 0.75, 0.75],
            ],
            cell=[5.43, 5.43, 5.43],
            pbc=True,
        )
        atoms.info["label"] = "Si-001"
        atoms.info["pressure"] = -0.05
        atoms.info["spin"] = 0.0
        atoms.info["spin_abs"] = 0.0
        atoms.info["symm"] = "Fd-3m"
        atoms.info["copies"] = 1

        energy = -42.5

        if with_forces:
            forces = np.array(
                [
                    [0.01, -0.02, 0.03],
                    [-0.01, 0.02, -0.03],
                    [0.005, -0.005, 0.005],
                    [-0.005, 0.005, -0.005],
                ]
            )
            calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
        else:
            calc = SinglePointCalculator(atoms, energy=energy)
        atoms.calc = calc

        if with_spin:
            atoms.set_initial_magnetic_moments([2.5, -2.5, 1.0, -1.0])

        return atoms

    def test_basic(self, tmp_path):
        xyz_file = tmp_path / "test.xyz"
        out_dir = tmp_path / "res_output"

        atoms = self._make_si_atoms()
        write(str(xyz_file), atoms, format="extxyz")

        n = extxyz_to_res(xyz_file, out_dir)
        assert n == 1

        content = _read_single_res_from_dir(out_dir, "Si-001")
        assert "TITL" in content
        assert "Si-001" in content
        assert "CELL" in content
        assert "END" in content

        # Parse with _read_res to verify cryan compatibility
        lines = content.strip().splitlines()
        data = _read_res(lines)
        assert data["titl"] is not None
        assert data["titl"].label == "Si-001"
        assert data["titl"].pressure == pytest.approx(-0.05, abs=0.01)
        assert len(data["species"]) == 4

    def test_with_forces(self, tmp_path):
        xyz_file = tmp_path / "test.xyz"
        out_dir = tmp_path / "res_output"

        atoms = self._make_si_atoms(with_forces=True)
        write(str(xyz_file), atoms, format="extxyz")

        extxyz_to_res(xyz_file, out_dir)

        content = _read_single_res_from_dir(out_dir, "Si-001")
        lines = content.strip().splitlines()
        atom_lines = [line.split() for line in lines if line.startswith("Si")]
        assert all(len(tokens) == 10 for tokens in atom_lines)
        assert all(float(tokens[6]) == pytest.approx(0.0) for tokens in atom_lines)
        forces = _parse_res_forces(lines)
        assert forces is not None
        assert len(forces) == 4
        assert forces[0] == pytest.approx([0.01, -0.02, 0.03], abs=1e-4)

    def test_with_spin(self, tmp_path):
        xyz_file = tmp_path / "test.xyz"
        out_dir = tmp_path / "res_output"

        atoms = self._make_si_atoms(with_spin=True)
        write(str(xyz_file), atoms, format="extxyz")

        extxyz_to_res(xyz_file, out_dir)

        content = _read_single_res_from_dir(out_dir, "Si-001")
        lines = content.strip().splitlines()
        data = _read_res(lines)
        assert len(data["spins"]) == 4
        assert data["spins"][0] == pytest.approx(2.5)

    def test_cryan_compatible(self, tmp_path):
        """Verify output .res is parseable by cryan-style parsers."""
        xyz_file = tmp_path / "test.xyz"
        out_dir = tmp_path / "res_output"

        atoms = self._make_si_atoms()
        write(str(xyz_file), atoms, format="extxyz")

        extxyz_to_res(xyz_file, out_dir)

        content = _read_single_res_from_dir(out_dir, "Si-001")
        lines = content.strip().splitlines()

        # Parse TITL line — cryan expects standard format
        titl_line = [line for line in lines if "TITL" in line][0]
        ti = parse_titl(titl_line)
        assert ti.label == "Si-001"
        assert ti.natoms == 4
        assert "n - 1" in titl_line

    def test_multiple_structures_unpacked(self, tmp_path):
        """Multiple structures become individual .res files."""
        out_dir = tmp_path / "res_output"
        xyz_file = tmp_path / "multi.xyz"

        atoms_list = []
        for i in range(3):
            atoms = Atoms("Si2", positions=[[0, 0, 0], [1.3, 1.3, 1.3]], cell=[5.43] * 3, pbc=True)
            atoms.info["label"] = f"Si-{i:03d}"
            atoms.info["pressure"] = 0.0
            atoms.info["spin"] = 0.0
            atoms.info["spin_abs"] = 0.0
            atoms.info["symm"] = "P1"
            atoms.info["copies"] = 1
            calc = SinglePointCalculator(atoms, energy=-40.0 - i)
            atoms.calc = calc
            atoms_list.append(atoms)

        write(str(xyz_file), atoms_list, format="extxyz")
        n = extxyz_to_res(xyz_file, out_dir)
        assert n == 3

        # Each structure should have its own file
        for i in range(3):
            content = _read_single_res_from_dir(out_dir, f"Si-{i:03d}")
            assert f"Si-{i:03d}" in content


# ---------------------------------------------------------------------------
# Round-trip tests
# ---------------------------------------------------------------------------


class TestRoundTrip:
    def test_res_to_xyz_to_res(self, tmp_path):
        """res → extxyz → res: verify data preserved."""
        res_in = tmp_path / "in.res"
        xyz_mid = tmp_path / "mid.xyz"
        res_out_dir = tmp_path / "res_output"

        res_in.write_text(RES_WITH_FORCES)

        res_to_extxyz(res_in, xyz_mid)
        extxyz_to_res(xyz_mid, res_out_dir)

        # Parse input
        in_data = _read_res(RES_WITH_FORCES.strip().splitlines())
        in_forces = _parse_res_forces(RES_WITH_FORCES.strip().splitlines())

        # Parse output
        content = _read_single_res_from_dir(res_out_dir, "Si-002")
        out_data = _read_res(content.strip().splitlines())
        out_forces = _parse_res_forces(content.strip().splitlines())

        # TITL metadata
        assert out_data["titl"].label == in_data["titl"].label
        assert out_data["titl"].pressure == pytest.approx(in_data["titl"].pressure, abs=0.01)
        assert out_data["titl"].enthalpy == pytest.approx(in_data["titl"].enthalpy, abs=0.001)
        assert out_data["titl"].natoms == in_data["titl"].natoms

        # Positions
        for i in range(4):
            for j in range(3):
                assert out_data["scaled_positions"][i][j] == pytest.approx(
                    in_data["scaled_positions"][i][j], abs=1e-6
                )

        # Forces
        assert out_forces is not None
        for i in range(4):
            assert out_forces[i] == pytest.approx(in_forces[i], abs=1e-4)

    def test_xyz_to_res_to_xyz(self, tmp_path):
        """extxyz → res → extxyz: verify energy and forces preserved."""
        xyz_in = tmp_path / "in.xyz"
        res_mid_dir = tmp_path / "res_output"
        # Collect the single .res back into a packed file for res_to_extxyz
        packed_res = tmp_path / "packed.res"
        xyz_out = tmp_path / "out.xyz"

        atoms = Atoms(
            "Si2",
            positions=[[0, 0, 0], [1.3575, 1.3575, 1.3575]],
            cell=[5.43, 5.43, 5.43],
            pbc=True,
        )
        forces = np.array([[0.01, -0.02, 0.03], [-0.01, 0.02, -0.03]])
        calc = SinglePointCalculator(atoms, energy=-42.5, forces=forces)
        atoms.calc = calc
        atoms.info["label"] = "Si-test"
        atoms.info["pressure"] = 1.5
        atoms.info["spin"] = 2.0
        atoms.info["spin_abs"] = 1.0
        atoms.info["symm"] = "P1"
        atoms.info["copies"] = 5

        write(str(xyz_in), atoms, format="extxyz")

        # extxyz → unpacked res
        extxyz_to_res(xyz_in, res_mid_dir)

        # Collect unpacked res into a packed file for round-trip
        res_content = _read_single_res_from_dir(res_mid_dir, "Si-test")
        packed_res.write_text(res_content)

        # packed res → extxyz
        res_to_extxyz(packed_res, xyz_out)

        atoms2 = read(str(xyz_out))
        assert atoms2.get_potential_energy() == pytest.approx(-42.5, abs=0.001)
        assert atoms2.get_forces()[0] == pytest.approx([0.01, -0.02, 0.03], abs=1e-4)
        assert atoms2.info["pressure"] == pytest.approx(1.5, abs=0.01)
        assert atoms2.info["label"] == "Si-test"
        assert atoms2.info["copies"] == 5

    def test_packed_roundtrip(self, tmp_path):
        """Multiple structures round-trip correctly."""
        res_in = tmp_path / "in.res"
        xyz_mid = tmp_path / "mid.xyz"
        res_out_dir = tmp_path / "res_output"

        res_in.write_text(PACKED_RES)

        res_to_extxyz(res_in, xyz_mid)
        extxyz_to_res(xyz_mid, res_out_dir)

        # Verify extxyz has both structures
        atoms_list = read(str(xyz_mid), index=":")
        assert len(atoms_list) == 2
        assert atoms_list[0].info["label"] == "Si-001"
        assert atoms_list[1].info["label"] == "Si-002"

        # Verify unpacked output: two individual .res files
        content1 = _read_single_res_from_dir(res_out_dir, "Si-001")
        content2 = _read_single_res_from_dir(res_out_dir, "Si-002")
        assert "TITL" in content1
        assert "TITL" in content2
        assert "Si-001" in content1
        assert "Si-002" in content2


# ---------------------------------------------------------------------------
# Extract single structure by label
# ---------------------------------------------------------------------------


class TestExtractStructure:
    def test_extract_from_packed_res_to_res(self, tmp_path):
        """Extract one structure from packed .res → single .res."""
        packed = tmp_path / "packed.res"
        out = tmp_path / "Si-002.res"
        packed.write_text(PACKED_RES)

        found = extract_structure(packed, "Si-002", out)
        assert found is True
        assert out.exists()

        content = out.read_text()
        assert "Si-002" in content
        assert "Si-001" not in content
        assert "END" in content

    def test_extract_from_packed_res_to_extxyz(self, tmp_path):
        """Extract one structure from packed .res → extxyz."""
        packed = tmp_path / "packed.res"
        out = tmp_path / "Si-001.xyz"
        packed.write_text(PACKED_RES)

        found = extract_structure(packed, "Si-001", out)
        assert found is True

        atoms = read(str(out))
        assert atoms.info["label"] == "Si-001"
        assert len(atoms) == 4

    def test_extract_from_extxyz_to_res(self, tmp_path):
        """Extract one structure from extxyz → .res."""
        xyz_file = tmp_path / "multi.xyz"
        out = tmp_path / "Si-002.res"

        atoms_list = []
        for name in ["Si-001", "Si-002", "Si-003"]:
            atoms = Atoms("Si2", positions=[[0, 0, 0], [1.3, 1.3, 1.3]], cell=[5.43] * 3, pbc=True)
            atoms.info["label"] = name
            atoms.info["pressure"] = 0.0
            atoms.info["spin"] = 0.0
            atoms.info["spin_abs"] = 0.0
            atoms.info["symm"] = "P1"
            atoms.info["copies"] = 1
            calc = SinglePointCalculator(atoms, energy=-40.0)
            atoms.calc = calc
            atoms_list.append(atoms)
        write(str(xyz_file), atoms_list, format="extxyz")

        found = extract_structure(xyz_file, "Si-002", out)
        assert found is True

        content = out.read_text()
        assert "Si-002" in content
        data = _read_res(content.strip().splitlines())
        assert data["titl"].label == "Si-002"

    def test_extract_from_extxyz_to_extxyz(self, tmp_path):
        """Extract one structure from extxyz → extxyz."""
        xyz_file = tmp_path / "multi.xyz"
        out = tmp_path / "extracted.xyz"

        atoms_list = []
        for name in ["Fe-001", "Fe-002"]:
            atoms = Atoms("Fe2", positions=[[0, 0, 0], [1.4, 1.4, 1.4]], cell=[2.87] * 3, pbc=True)
            atoms.info["label"] = name
            atoms.info["pressure"] = 0.0
            atoms.info["spin"] = 2.0
            atoms.info["spin_abs"] = 1.0
            atoms.info["symm"] = "Im-3m"
            atoms.info["copies"] = 1
            forces = np.array([[0.1, -0.2, 0.3], [-0.1, 0.2, -0.3]])
            calc = SinglePointCalculator(atoms, energy=-25.0, forces=forces)
            atoms.calc = calc
            atoms_list.append(atoms)
        write(str(xyz_file), atoms_list, format="extxyz")

        found = extract_structure(xyz_file, "Fe-002", out)
        assert found is True

        atoms = read(str(out))
        assert atoms.info["label"] == "Fe-002"
        assert atoms.get_potential_energy() == pytest.approx(-25.0)
        assert atoms.get_forces()[0] == pytest.approx([0.1, -0.2, 0.3], abs=1e-4)

    def test_extract_not_found(self, tmp_path):
        """Label not present returns False."""
        packed = tmp_path / "packed.res"
        out = tmp_path / "out.res"
        packed.write_text(PACKED_RES)

        found = extract_structure(packed, "nonexistent", out)
        assert found is False
        assert not out.exists()

    def test_extract_with_forces_preserved(self, tmp_path):
        """Forces survive extraction from packed .res."""
        packed = tmp_path / "packed.res"
        out = tmp_path / "Si-002.res"
        packed.write_text(PACKED_RES)

        extract_structure(packed, "Si-002", out)

        content = out.read_text()
        forces = _parse_res_forces(content.strip().splitlines())
        assert forces is not None
        assert forces[0] == pytest.approx([0.01, -0.02, 0.03], abs=1e-4)

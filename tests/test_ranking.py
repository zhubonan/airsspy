"""Tests for the ranking module."""

import io
import os
import tempfile

import pytest

from airsspy.ranking import (
    StructureRecord,
    _compute_distance_fingerprint,
    _parse_res_fast,
    _reduce_formula,
    eliminate_similar,
    format_header,
    format_rank_line,
    rank_structures,
    read_res_file,
    read_res_stream,
    summary_structures,
)


# ---------------------------------------------------------------------------
# Sample RES data
# ---------------------------------------------------------------------------

RES_SI_1 = """\
TITL Si-001 -0.05 40.0 -42.500000 0 0 4 (Fd-3m) n - 1
CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0
LATT -1
SFAC Si
Si     1  0.0000000000  0.0000000000  0.0000000000 1.0
Si     1  0.2500000000  0.2500000000  0.2500000000 1.0
Si     1  0.5000000000  0.5000000000  0.5000000000 1.0
Si     1  0.7500000000  0.7500000000  0.7500000000 1.0
END
"""

RES_SI_2 = """\
TITL Si-002 -0.05 42.0 -43.200000 0 0 4 (Pm-3m) n - 1
CELL 1.0  5.50 5.50 5.50 90.0 90.0 90.0
LATT -1
SFAC Si
Si     1  0.0000000000  0.0000000000  0.0000000000 1.0
Si     1  0.5000000000  0.5000000000  0.5000000000 1.0
Si     1  0.2500000000  0.2500000000  0.2500000000 1.0
Si     1  0.7500000000  0.7500000000  0.7500000000 1.0
END
"""

RES_SIO2 = """\
TITL SiO2-001 0.0 45.0 -80.500000 0 0 6 (P1) n - 1
CELL 1.0  4.0 5.0 6.0 90.0 90.0 90.0
LATT -1
SFAC Si O
Si     1  0.0  0.0  0.0  1.0
O      2  0.5  0.5  0.5  1.0
O      2  0.3  0.3  0.3  1.0
Si     1  0.6  0.6  0.6  1.0
O      2  0.1  0.1  0.1  1.0
O      2  0.8  0.8  0.8  1.0
END
"""

RES_SI_COPIES = """\
TITL Si-003 0.0 41.0 -42.800000 0 0 4 (P1) n - 5
CELL 1.0  5.45 5.45 5.45 90.0 90.0 90.0
LATT -1
SFAC Si
Si     1  0.0  0.0  0.0  1.0
Si     1  0.25 0.25 0.25 1.0
Si     1  0.5  0.5  0.5  1.0
Si     1  0.75 0.75 0.75 1.0
END
"""

RES_WITH_SPIN = """\
TITL Fe-001 0.0 30.0 -25.000000 4.0 2.0 2 (Im-3m) n - 1
CELL 1.0  2.87 2.87 2.87 90.0 90.0 90.0
LATT -1
SFAC Fe
Fe     1  0.0  0.0  0.0  1.0
Fe     1  0.5  0.5  0.5  1.0
END
"""


# ---------------------------------------------------------------------------
# StructureRecord tests
# ---------------------------------------------------------------------------


class TestStructureRecord:
    def test_n_formula_units_single_element(self):
        rec = StructureRecord(
            label="Si-001",
            pressure=0.0,
            volume=40.0,
            enthalpy=-42.5,
            spin=0.0,
            spin_abs=0.0,
            natoms=4,
            symm="(Fd-3m)",
            species_counts={"Si": 4},
        )
        assert rec.n_formula_units == 4
        assert rec.enthalpy_per_fu == pytest.approx(-42.5 / 4)

    def test_n_formula_units_compound(self):
        rec = StructureRecord(
            label="SiO2-001",
            pressure=0.0,
            volume=45.0,
            enthalpy=-80.5,
            spin=0.0,
            spin_abs=0.0,
            natoms=6,
            symm="(P1)",
            species_counts={"Si": 2, "O": 4},
        )
        assert rec.n_formula_units == 2
        assert rec.reduced_formula == "SiO2"
        assert rec.enthalpy_per_fu == pytest.approx(-80.5 / 2)

    def test_hill_system_with_carbon(self):
        rec = StructureRecord(
            label="test",
            pressure=0.0,
            volume=1.0,
            enthalpy=0.0,
            species_counts={"O": 2, "C": 1, "H": 4},
        )
        # C first, H second, O last (cryan convention)
        assert rec.reduced_formula == "CH4O2"

    def test_hill_system_no_carbon(self):
        rec = StructureRecord(
            label="test",
            pressure=0.0,
            volume=1.0,
            enthalpy=0.0,
            species_counts={"Cl": 2, "Na": 2},
        )
        # By atomic number: Na(11) before Cl(17)
        assert rec.reduced_formula == "NaCl"

    def test_volume_per_fu(self):
        rec = StructureRecord(
            label="test",
            pressure=0.0,
            volume=45.0,
            enthalpy=-80.5,
            species_counts={"Si": 2, "O": 4},
        )
        assert rec.volume_per_fu == pytest.approx(45.0 / 2)


# ---------------------------------------------------------------------------
# Formula tests
# ---------------------------------------------------------------------------


class TestReduceFormula:
    def test_single_element(self):
        assert _reduce_formula({"Si": 4}) == "Si4"

    def test_binary(self):
        # Si(14) before O(last) → SiO2
        assert _reduce_formula({"Si": 2, "O": 4}) == "SiO2"

    def test_ternary(self):
        # Li(3), Fe(26), O(last) → LiFeO2
        assert _reduce_formula({"Li": 2, "Fe": 2, "O": 4}) == "LiFeO2"

    def test_carbon_first(self):
        # C first, H second, O last
        assert _reduce_formula({"H": 4, "O": 1, "C": 1}) == "CH4O"

    def test_empty(self):
        assert _reduce_formula({}) == ""


# ---------------------------------------------------------------------------
# Parsing tests
# ---------------------------------------------------------------------------


class TestParsing:
    def test_parse_res_fast_si(self):
        lines = RES_SI_1.strip().splitlines()
        rec = _parse_res_fast(lines)
        assert rec is not None
        assert rec.label == "Si-001"
        assert rec.pressure == -0.05
        assert rec.volume == 40.0
        assert rec.enthalpy == -42.5
        assert rec.natoms == 4
        assert rec.species_counts == {"Si": 4}
        assert rec.copies == 1
        assert rec.symm == "(Fd-3m)"

    def test_parse_res_fast_sio2(self):
        lines = RES_SIO2.strip().splitlines()
        rec = _parse_res_fast(lines)
        assert rec is not None
        assert rec.label == "SiO2-001"
        assert rec.species_counts == {"Si": 2, "O": 4}
        assert rec.reduced_formula == "SiO2"

    def test_parse_res_fast_copies(self):
        lines = RES_SI_COPIES.strip().splitlines()
        rec = _parse_res_fast(lines)
        assert rec is not None
        assert rec.copies == 5

    def test_parse_res_fast_spin(self):
        lines = RES_WITH_SPIN.strip().splitlines()
        rec = _parse_res_fast(lines)
        assert rec is not None
        assert rec.spin == 4.0
        assert rec.spin_abs == 2.0

    def test_parse_empty_returns_none(self):
        rec = _parse_res_fast([])
        assert rec is None

    def test_read_res_stream(self):
        packed = RES_SI_1 + RES_SI_2
        stream = io.StringIO(packed)
        records = read_res_stream(stream)
        assert len(records) == 2
        assert records[0].label == "Si-001"
        assert records[1].label == "Si-002"

    def test_read_res_stream_multiple_compositions(self):
        packed = RES_SI_1 + RES_SIO2
        stream = io.StringIO(packed)
        records = read_res_stream(stream)
        assert len(records) == 2
        assert records[0].reduced_formula == "Si4"
        assert records[1].reduced_formula == "SiO2"


# ---------------------------------------------------------------------------
# Ranking tests
# ---------------------------------------------------------------------------


class TestRanking:
    def _make_si_record(self, label, enthalpy, copies=1):
        return StructureRecord(
            label=label,
            pressure=-0.05,
            volume=40.0,
            enthalpy=enthalpy,
            natoms=4,
            species_counts={"Si": 4},
            copies=copies,
        )

    def test_rank_single_composition(self):
        rec1 = self._make_si_record("Si-001", -42.5)
        rec2 = self._make_si_record("Si-002", -43.2)
        ranked = rank_structures([rec1, rec2])
        assert len(ranked) == 2
        # Most stable first
        assert ranked[0]["label"] == "Si-002"
        # First entry shows absolute
        assert ranked[0]["display_enthalpy"] == pytest.approx(-43.2 / 4)
        # Second entry shows relative
        assert ranked[1]["display_enthalpy"] > 0

    def test_rank_absolute_mode(self):
        rec1 = self._make_si_record("Si-001", -42.5)
        rec2 = self._make_si_record("Si-002", -43.2)
        ranked = rank_structures([rec1, rec2], absolute=True)
        assert ranked[0]["display_enthalpy"] == pytest.approx(-43.2 / 4)
        assert ranked[1]["display_enthalpy"] == pytest.approx(-42.5 / 4)

    def test_rank_delta_e_filter(self):
        rec1 = self._make_si_record("Si-001", -42.5)
        rec2 = self._make_si_record("Si-002", -30.0)
        # per_fu difference is 12.5/4 = 3.125 eV/atom per formula unit
        # relative_enthalpy_per_atom = (3.125 * 4) / 4 = 3.125 eV/atom
        ranked = rank_structures([rec1, rec2], delta_e=0.5)
        assert len(ranked) == 1
        assert ranked[0]["label"] == "Si-001"

    def test_rank_formula_filter(self):
        rec_si = self._make_si_record("Si-001", -42.5)
        rec_sio2 = StructureRecord(
            label="SiO2-001",
            pressure=0.0,
            volume=45.0,
            enthalpy=-80.5,
            natoms=6,
            species_counts={"Si": 2, "O": 4},
        )
        ranked = rank_structures([rec_si, rec_sio2], formula_filter="SiO2")
        assert len(ranked) == 1
        assert ranked[0]["formula"] == "SiO2"

    def test_rank_top_n(self):
        records = [self._make_si_record(f"Si-{i}", -40.0 - i) for i in range(10)]
        ranked = rank_structures(records, top_n=3)
        assert len(ranked) == 3

    def test_rank_multi_composition(self):
        rec_si = self._make_si_record("Si-001", -42.5)
        rec_sio2 = StructureRecord(
            label="SiO2-001",
            pressure=0.0,
            volume=45.0,
            enthalpy=-80.5,
            natoms=6,
            species_counts={"Si": 2, "O": 4},
        )
        ranked = rank_structures([rec_si, rec_sio2])
        # Sorted globally by enthalpy_per_fu
        assert len(ranked) == 2

    def test_rank_preserves_copies(self):
        rec = self._make_si_record("Si-003", -42.8, copies=5)
        ranked = rank_structures([rec])
        assert ranked[0]["copies"] == 5


# ---------------------------------------------------------------------------
# Summary tests
# ---------------------------------------------------------------------------


class TestSummary:
    def test_summary_returns_one_per_composition(self):
        rec1 = StructureRecord(
            label="Si-001", pressure=0.0, volume=40.0, enthalpy=-42.5,
            natoms=4, species_counts={"Si": 4}, copies=1,
        )
        rec2 = StructureRecord(
            label="Si-002", pressure=0.0, volume=42.0, enthalpy=-43.2,
            natoms=4, species_counts={"Si": 4}, copies=1,
        )
        rec3 = StructureRecord(
            label="SiO2-001", pressure=0.0, volume=45.0, enthalpy=-80.5,
            natoms=6, species_counts={"Si": 2, "O": 4}, copies=1,
        )
        ranked, total = summary_structures([rec1, rec2, rec3])
        # One per composition
        assert len(ranked) == 2
        # SiO2 (most negative per fu) and Si should both appear
        formulas = {r["formula"] for r in ranked}
        assert "SiO2" in formulas
        assert "Si4" in formulas
        # Total copies
        assert total == 3

    def test_summary_accumulates_copies(self):
        rec1 = StructureRecord(
            label="Si-001", pressure=0.0, volume=40.0, enthalpy=-42.5,
            natoms=4, species_counts={"Si": 4}, copies=3,
        )
        rec2 = StructureRecord(
            label="Si-002", pressure=0.0, volume=42.0, enthalpy=-43.2,
            natoms=4, species_counts={"Si": 4}, copies=2,
        )
        ranked, total = summary_structures([rec1, rec2])
        assert len(ranked) == 1
        assert ranked[0]["label"] == "Si-002"  # most stable
        assert ranked[0]["group_copies"] == 5
        assert total == 5


# ---------------------------------------------------------------------------
# Formatting tests
# ---------------------------------------------------------------------------


class TestFormatting:
    def test_format_rank_line_basic(self):
        rec = {
            "label": "Si-001",
            "pressure": -0.05,
            "volume_per_fu": 10.0,
            "display_enthalpy": -10.625,
            "nfu": 4,
            "formula": "Si",
            "symm": "(Fd-3m)",
            "copies": 1,
        }
        line = format_rank_line(rec)
        assert "Si-001" in line
        assert "Fd-3m" in line
        assert "(" not in line

    def test_format_rank_line_with_spin(self):
        rec = {
            "label": "Fe-001",
            "pressure": 0.0,
            "volume_per_fu": 15.0,
            "display_enthalpy": -12.5,
            "spin_per_fu": 2.0,
            "spin_abs_per_fu": 1.0,
            "nfu": 2,
            "formula": "Fe",
            "symm": "(Im-3m)",
            "copies": 1,
        }
        line = format_rank_line(rec, show_spin=True)
        assert "Fe-001" in line
        assert "Im-3m" in line
        assert "(" not in line

    def test_format_rank_line_long_label(self):
        rec = {
            "label": "very-long-structure-name-that-exceeds-20-chars",
            "pressure": 0.0,
            "volume_per_fu": 10.0,
            "display_enthalpy": -10.0,
            "nfu": 1,
            "formula": "Si",
            "symm": "(P1)",
            "copies": 1,
        }
        line_short = format_rank_line(rec, long_labels=False)
        line_long = format_rank_line(rec, long_labels=True)
        assert len(line_long) > len(line_short)


# ---------------------------------------------------------------------------
# File I/O tests
# ---------------------------------------------------------------------------


class TestFileIO:
    def test_read_res_file_single(self, tmp_path):
        """Test reading a single .res file."""
        res_file = tmp_path / "test.res"
        res_file.write_text(RES_SI_1)
        records = read_res_file(str(res_file))
        assert len(records) == 1
        assert records[0].label == "Si-001"
        assert records[0].source == str(res_file)

    def test_read_res_file_packed(self, tmp_path):
        """Test reading a packed .res file with multiple structures."""
        res_file = tmp_path / "packed.res"
        res_file.write_text(RES_SI_1 + RES_SI_2)
        records = read_res_file(str(res_file))
        assert len(records) == 2
        assert records[0].label == "Si-001"
        assert records[1].label == "Si-002"

    def test_read_res_file_no_end(self, tmp_path):
        """Test reading a .res file without trailing END."""
        res_file = tmp_path / "noend.res"
        res_file.write_text(RES_SI_1.replace("END\n", ""))
        records = read_res_file(str(res_file))
        assert len(records) == 1
        assert records[0].label == "Si-001"


# ---------------------------------------------------------------------------
# Header formatting tests
# ---------------------------------------------------------------------------


class TestFormatHeader:
    def test_header_default(self):
        header = format_header()
        assert "structure" in header
        assert "P/GPa" in header
        assert "V/A^3" in header
        assert "H/eV" in header
        assert "nfu" in header
        assert "formula" in header

    def test_header_with_spin(self):
        header = format_header(show_spin=True)
        assert "S" in header
        assert "|S|" in header

    def test_header_summary_mode(self):
        header = format_header(summary_mode=True)
        assert "tot#" in header

    def test_header_summary_with_spin(self):
        header = format_header(show_spin=True, summary_mode=True)
        assert "S" in header
        assert "tot#" in header


# ---------------------------------------------------------------------------
# Summary formatting tests
# ---------------------------------------------------------------------------


class TestFormatSummaryLine:
    def test_format_rank_line_summary(self):
        rec = {
            "label": "Si-001",
            "pressure": -0.05,
            "volume_per_fu": 10.0,
            "display_enthalpy": -10.625,
            "spin_per_fu": 0.0,
            "spin_abs_per_fu": 0.0,
            "nfu": 4,
            "formula": "Si",
            "symm": "(Fd-3m)",
            "copies": 3,
            "group_copies": 5,
        }
        line = format_rank_line(rec, summary_mode=True)
        assert "Si-001" in line
        # Should have two number columns (copies and group_copies)
        assert "3" in line
        assert "5" in line


# ---------------------------------------------------------------------------
# Eliminate similar tests
# ---------------------------------------------------------------------------


class TestEliminateSimilar:
    def test_eliminate_identical_structures(self):
        """Identical structures should be merged."""
        lines1 = [
            "TITL Si-001 -0.05 40.0 -42.500000 0 0 4 (Fd-3m) n - 1",
            "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Si",
            "Si     1  0.0000000000  0.0000000000  0.0000000000 1.0",
            "Si     1  0.2500000000  0.2500000000  0.2500000000 1.0",
            "Si     1  0.5000000000  0.5000000000  0.5000000000 1.0",
            "Si     1  0.7500000000  0.7500000000  0.7500000000 1.0",
        ]
        lines2 = [
            "TITL Si-002 -0.05 40.0 -43.000000 0 0 4 (Fd-3m) n - 1",
            "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Si",
            "Si     1  0.0000000000  0.0000000000  0.0000000000 1.0",
            "Si     1  0.2500000000  0.2500000000  0.2500000000 1.0",
            "Si     1  0.5000000000  0.5000000000  0.5000000000 1.0",
            "Si     1  0.7500000000  0.7500000000  0.7500000000 1.0",
        ]
        rec1 = _parse_res_fast(lines1)
        rec2 = _parse_res_fast(lines2)
        assert rec1 is not None
        assert rec2 is not None

        result = eliminate_similar([rec1, rec2], threshold=0.1)
        # Identical coordinates should merge
        assert len(result) == 1
        assert result[0].copies == 2

    def test_eliminate_different_structures(self):
        """Very different structures should not merge."""
        lines1 = [
            "TITL Si-001 -0.05 40.0 -42.500000 0 0 4 (P1) n - 1",
            "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Si",
            "Si     1  0.0  0.0  0.0  1.0",
            "Si     1  0.1  0.1  0.1  1.0",
            "Si     1  0.2  0.2  0.2  1.0",
            "Si     1  0.3  0.3  0.3  1.0",
        ]
        lines2 = [
            "TITL Si-002 -0.05 40.0 -43.000000 0 0 4 (P1) n - 1",
            "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Si",
            "Si     1  0.0  0.0  0.0  1.0",
            "Si     1  0.5  0.5  0.5  1.0",
            "Si     1  0.7  0.7  0.7  1.0",
            "Si     1  0.9  0.9  0.9  1.0",
        ]
        rec1 = _parse_res_fast(lines1)
        rec2 = _parse_res_fast(lines2)
        assert rec1 is not None
        assert rec2 is not None

        result = eliminate_similar([rec1, rec2], threshold=0.01)
        # Very different distances, small threshold → no merge
        assert len(result) == 2

    def test_eliminate_no_raw_lines(self):
        """Records without raw_lines should not merge (fingerprint=None)."""
        rec1 = StructureRecord(
            label="Si-001", pressure=0.0, volume=40.0, enthalpy=-42.5,
            natoms=4, species_counts={"Si": 4}, copies=1,
        )
        rec2 = StructureRecord(
            label="Si-002", pressure=0.0, volume=40.0, enthalpy=-43.0,
            natoms=4, species_counts={"Si": 4}, copies=1,
        )
        # No _raw_lines, so fingerprint will be None
        result = eliminate_similar([rec1, rec2], threshold=0.1)
        assert len(result) == 2


# ---------------------------------------------------------------------------
# Fingerprint tests
# ---------------------------------------------------------------------------


class TestFingerprint:
    def test_fingerprint_includes_periodic_images(self):
        """Fingerprint should include periodic image distances for small cells."""
        # BCC Fe with a=2.87 Å — many periodic images within 4 Å
        lines = [
            "TITL Fe-001 0.0 30.0 -25.000000 0 0 2 (Im-3m) n - 1",
            "CELL 1.0  2.87 2.87 2.87 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Fe",
            "Fe     1  0.0  0.0  0.0  1.0",
            "Fe     1  0.5  0.5  0.5  1.0",
        ]
        rec = _parse_res_fast(lines)
        assert rec is not None
        fp = _compute_distance_fingerprint(rec, cutoff=4.0, zweight=False)
        assert fp is not None
        # With a 2.87 Å cell, get_all_neighbors(4.0) should find multiple
        # periodic images — more than the single minimum-image pair
        assert len(fp) > 1

    def test_fingerprint_zweighting(self):
        """Z-weighting should scale distances by atomic number."""
        lines = [
            "TITL SiO2-001 0.0 120.0 -80.5 0 0 6 (P1) n - 1",
            "CELL 1.0  5.0 5.0 5.0 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Si O",
            "Si     1  0.0  0.0  0.0  1.0",
            "O      2  0.5  0.5  0.5  1.0",
            "O      2  0.3  0.3  0.3  1.0",
            "Si     1  0.6  0.6  0.6  1.0",
            "O      2  0.1  0.1  0.1  1.0",
            "O      2  0.8  0.8  0.8  1.0",
        ]
        rec = _parse_res_fast(lines)
        assert rec is not None
        fp_unweighted = _compute_distance_fingerprint(rec, cutoff=4.0, zweight=False)
        fp_weighted = _compute_distance_fingerprint(rec, cutoff=4.0, zweight=True)
        assert fp_unweighted is not None
        assert fp_weighted is not None
        assert len(fp_unweighted) == len(fp_weighted)
        # Z-weighted distances should differ from unweighted
        assert fp_weighted != fp_unweighted
        # For SiO2: zmax=14 (Si). Si-Si weight = 14²/(14·14) = 1.0,
        # Si-O weight = 14²/(14·8) = 1.75, O-O weight = 14²/(8·8) = 3.0625.
        # All weighted distances should be >= corresponding unweighted distances.
        for w, u in zip(fp_weighted, fp_unweighted):
            assert w >= u - 1e-10

    def test_fingerprint_no_raw_lines_returns_none(self):
        """Records without raw_lines should return None."""
        rec = StructureRecord(
            label="Si-001", pressure=0.0, volume=40.0, enthalpy=-42.5,
            natoms=4, species_counts={"Si": 4},
        )
        assert _compute_distance_fingerprint(rec) is None


# ---------------------------------------------------------------------------
# Edge case parsing tests
# ---------------------------------------------------------------------------


class TestParsingEdgeCases:
    def test_parse_titl_bad_pressure(self):
        """TITL with non-numeric pressure should default to 0.0."""
        lines = [
            "TITL Si-001 bad 40.0 -42.5 0 0 4 (P1) n - 1",
            "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Si",
            "Si     1  0.0  0.0  0.0  1.0",
        ]
        rec = _parse_res_fast(lines)
        assert rec is not None
        assert rec.pressure == 0.0

    def test_parse_titl_bad_volume(self):
        """TITL with non-numeric volume should default to 0.0."""
        lines = [
            "TITL Si-001 -0.05 bad -42.5 0 0 4 (P1) n - 1",
            "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Si",
            "Si     1  0.0  0.0  0.0  1.0",
        ]
        rec = _parse_res_fast(lines)
        assert rec is not None
        assert rec.volume == 0.0

    def test_parse_titl_bad_enthalpy(self):
        """TITL with non-numeric enthalpy should default to 0.0."""
        lines = [
            "TITL Si-001 -0.05 40.0 bad 0 0 4 (P1) n - 1",
            "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Si",
            "Si     1  0.0  0.0  0.0  1.0",
        ]
        rec = _parse_res_fast(lines)
        assert rec is not None
        assert rec.enthalpy == 0.0

    def test_parse_titl_bad_spin(self):
        """TITL with bad spin values should default to 0.0."""
        lines = [
            "TITL Si-001 -0.05 40.0 -42.5 bad bad 4 (P1) n - 1",
            "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Si",
            "Si     1  0.0  0.0  0.0  1.0",
        ]
        rec = _parse_res_fast(lines)
        assert rec is not None
        assert rec.spin == 0.0
        assert rec.spin_abs == 0.0

    def test_parse_titl_bad_natoms(self):
        """TITL with bad natoms should default to species count."""
        lines = [
            "TITL Si-001 -0.05 40.0 -42.5 0 0 bad (P1) n - 1",
            "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Si",
            "Si     1  0.0  0.0  0.0  1.0",
            "Si     1  0.5  0.5  0.5  1.0",
        ]
        rec = _parse_res_fast(lines)
        assert rec is not None
        assert rec.natoms == 2  # From species count

    def test_parse_titl_no_copies_suffix(self):
        """TITL without 'n - <N>' should default to copies=1."""
        lines = [
            "TITL Si-001 -0.05 40.0 -42.5 0 0 4 (P1)",
            "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Si",
            "Si     1  0.0  0.0  0.0  1.0",
        ]
        rec = _parse_res_fast(lines)
        assert rec is not None
        assert rec.copies == 1

    def test_parse_titl_short(self):
        """TITL with fewer than 5 fields after TITL keyword."""
        lines = [
            "TITL Si-001 -0.05",
            "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Si",
            "Si     1  0.0  0.0  0.0  1.0",
        ]
        rec = _parse_res_fast(lines)
        # Too few tokens → no label parsed, returns None
        assert rec is None

    def test_parse_titl_bad_copies(self):
        """TITL with 'n - bad' copies suffix should default to 1."""
        lines = [
            "TITL Si-001 -0.05 40.0 -42.5 0 0 4 (P1) n - bad",
            "CELL 1.0  5.43 5.43 5.43 90.0 90.0 90.0",
            "LATT -1",
            "SFAC Si",
            "Si     1  0.0  0.0  0.0  1.0",
        ]
        rec = _parse_res_fast(lines)
        assert rec is not None
        assert rec.copies == 1

    def test_n_formula_units_empty_counts(self):
        """Empty species_counts should return nfu=1."""
        rec = StructureRecord(
            label="test", pressure=0.0, volume=1.0, enthalpy=0.0,
            species_counts={},
        )
        assert rec.n_formula_units == 1

    def test_n_formula_units_single_species_count_one(self):
        """Single species with count 1 should return nfu=1."""
        rec = StructureRecord(
            label="test", pressure=0.0, volume=1.0, enthalpy=0.0,
            species_counts={"Si": 1},
        )
        assert rec.n_formula_units == 1

    def test_enthalpy_per_fu_zero_nfu(self):
        """Zero nfu edge case — should return raw enthalpy."""
        rec = StructureRecord(
            label="test", pressure=0.0, volume=1.0, enthalpy=-5.0,
            species_counts={},
        )
        # nfu=1 for empty counts, so this tests the normal path
        assert rec.enthalpy_per_fu == -5.0

"""Tests for analysis modules (collect and query)."""

import pytest

from airsspy.analysis.collect import (
    collect_res_in_df,
    get_minsep_range,
    read_ca,
    read_stream,
)
from airsspy.restools import RESFile

RES_EXAMPLE_1 = """\
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

RES_EXAMPLE_2 = """\
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

PACKED_RES = RES_EXAMPLE_1 + RES_EXAMPLE_2


def test_read_stream():
    """Test reading multiple RES from a stream."""
    lines = PACKED_RES.strip().splitlines()
    titl_list, atoms_list = read_stream(lines)
    assert len(titl_list) == 2
    assert len(atoms_list) == 2
    assert titl_list[0].label == "Si-001"
    assert titl_list[1].label == "Si-002"
    assert len(atoms_list[0]) == 4


def test_read_ca_no_spin():
    """Test reading ca command output without spin."""
    ca_lines = [
        "Si-001 -0.05 40.0 -42.5 1 Si4 (Fd-3m) 3",
        "Si-002 -0.05 42.0 -43.2 1 Si4 (Pm-3m) 1",
    ]
    df = read_ca(ca_lines)
    assert len(df) == 2
    assert df.iloc[0]["label"] == "Si-001"
    # Second row H should be relative (added to first row's H)
    assert df.iloc[1]["H"] == pytest.approx(-43.2 + -42.5)


def test_read_ca_with_spin():
    """Test reading ca command output with spin columns."""
    ca_lines = [
        "Fe-001 0.0 30.0 -50.0 2.0 1.5 2 Fe2 (Immm) 5",
        "Fe-002 0.0 31.0 -49.0 1.0 0.5 2 Fe2 (Pm-3m) 3",
    ]
    df = read_ca(ca_lines)
    assert len(df) == 2
    assert "spin" in df.columns
    assert "aspin" in df.columns
    assert df.iloc[0]["spin"] == 2.0


def test_read_ca_empty():
    """Test reading empty ca output."""
    df = read_ca([])
    assert df.empty


def test_collect_res_in_df():
    """Test collecting RESFile objects into a DataFrame."""
    r1 = RESFile.from_string(RES_EXAMPLE_1)
    r2 = RESFile.from_string(RES_EXAMPLE_2)

    df = collect_res_in_df([r1, r2], norm_mode="per_atom")
    assert len(df) == 2
    assert "H" in df.columns
    assert "V" in df.columns
    # Should be sorted by H (ascending)
    assert df.iloc[0]["H"] < df.iloc[1]["H"]


def test_collect_res_in_df_per_fu():
    """Test per-formula-unit normalisation."""
    r1 = RESFile.from_string(RES_EXAMPLE_1)
    df = collect_res_in_df([r1], norm_mode="per_formula_unit")
    assert len(df) == 1


def test_collect_res_empty():
    """Test collecting empty list."""
    df = collect_res_in_df([])
    assert df.empty


def test_get_minsep_range():
    """Test minsep range computation from ensemble."""
    minseps = [
        {"Si-Si": 2.3, "Si-O": 1.6},
        {"Si-Si": 2.5, "Si-O": 1.7},
        {"Si-Si": 2.1, "Si-O": 1.8},
    ]
    result = get_minsep_range(minseps)
    assert result["Si-Si"] == [2.1, 2.5]
    assert result["Si-O"] == [1.6, 1.8]


def test_get_minsep_range_with_cap():
    """Test minsep range with capping."""
    minseps = [
        {"Si-Si": 2.3},
        {"Si-Si": 0.5},  # Below cap
        {"Si-Si": 5.0},  # Above cap
    ]
    result = get_minsep_range(minseps, cap=(1.0, 4.0))
    assert result["Si-Si"][0] == 1.0  # Capped from below
    assert result["Si-Si"][1] == 4.0  # Capped from above


def test_resfile_from_packed(tmp_path):
    """Test reading packed RES file with multiple structures."""
    packed_file = tmp_path / "packed.res"
    packed_file.write_text(PACKED_RES)

    results = RESFile.from_packed(str(packed_file), include_structure=True)
    assert len(results) == 2
    assert results[0].label == "Si-001"
    assert results[1].label == "Si-002"


def test_resfile_from_packed_titl_only(tmp_path):
    """Test reading packed RES with only_titl=True."""
    packed_file = tmp_path / "packed.res"
    packed_file.write_text(PACKED_RES)

    results = RESFile.from_packed(str(packed_file), only_titl=True)
    assert len(results) == 2
    assert results[0].label == "Si-001"

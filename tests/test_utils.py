# -*- coding: utf-8 -*-
###########################################################################
# airss-ase                                                               #
# Copyright (C) 2019  Bonan Zhu                                           #
#                                                                         #
# This program is free software; you can redistribute it and/or modify    #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation; either version 2 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# This program is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License along #
# with this program; if not, write to the Free Software Foundation, Inc., #
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             #
###########################################################################
"""
Tests for the utils module
"""

import io
import pytest
from airsspy.utils import (
    calc_kpt_tuple_recip,
    count_pattern_in_file,
    extract_number_from_string,
    filter_out_stream,
    find_pattern_in_file,
    format_time_elapsed,
    safe_cast_float,
    safe_cast_int,
    stream_to_list,
    trim_stream,
    unique,
)
from ase import Atoms
import numpy as np


def test_trim_stream():
    """Test stream trimming functionality"""
    test_content = """Start of content
This should be kept
%BLOCK LATTICE
This should be removed
More removed content
%ENDBLOCK LATTICE
This should be kept again
End of content"""

    stream = io.StringIO(test_content)
    trimmed = trim_stream(stream, r"^%BLOCK [Ll][Aa][Tt]", r"^%ENDBLOCK [Ll][Aa][Tt]")

    result = trimmed.read()
    # The trim function extracts content between start and end patterns
    assert "This should be removed" in result
    assert "More removed content" in result
    assert "This should be kept" not in result
    assert "This should be kept again" not in result


def test_filter_out_stream():
    """Test stream filtering functionality"""
    test_content = """Start of content
This should be kept
%BLOCK LATTICE
This should be removed
More removed content
%ENDBLOCK LATTICE
This should be kept again
End of content"""

    stream = io.StringIO(test_content)
    filtered = filter_out_stream(stream, r"^%BLOCK [Ll][Aa][Tt]", r"^%ENDBLOCK [Ll][Aa][Tt]")

    result = filtered.read()
    assert "This should be kept" in result
    assert "This should be kept again" in result
    assert "This should be removed" not in result


def test_unique():
    """Test unique function"""
    items = [1, 2, 3, 2, 1, 4, 3, 5]
    result = unique(items)
    assert result == [1, 2, 3, 4, 5]

    # Test with strings
    str_items = ["Al", "Si", "Al", "C", "Si", "Al"]
    str_result = unique(str_items)
    assert str_result == ["Al", "Si", "C"]


def test_calc_kpt_tuple_recip():
    """Test k-point calculation"""
    # Create a simple cubic cell
    atoms = Atoms(cell=[4.0, 4.0, 4.0], pbc=True)

    # Test different spacings
    kpts = calc_kpt_tuple_recip(atoms, mp_spacing=0.25)
    assert len(kpts) == 3
    assert all(isinstance(k, int) for k in kpts)
    assert all(k > 0 for k in kpts)

    # Test rounding options
    kpts_up = calc_kpt_tuple_recip(atoms, mp_spacing=0.3, rounding="up")
    kpts_down = calc_kpt_tuple_recip(atoms, mp_spacing=0.3, rounding="down")

    # Should be roughly the same but potentially different due to rounding
    assert len(kpts_up) == 3
    assert len(kpts_down) == 3


def test_safe_cast_float():
    """Test safe float casting"""
    assert safe_cast_float("3.14") == 3.14
    assert safe_cast_float("-2.5") == -2.5
    assert safe_cast_float("invalid", 999.0) == 999.0
    assert safe_cast_float("", 0.0) == 0.0


def test_safe_cast_int():
    """Test safe integer casting"""
    assert safe_cast_int("42") == 42
    assert safe_cast_int("-7") == -7
    assert safe_cast_int("invalid", 999) == 999
    assert safe_cast_int("", 0) == 0


def test_extract_number_from_string():
    """Test number extraction from strings"""
    assert extract_number_from_string("Energy: -15.5 eV") == -15.5
    assert extract_number_from_string("Pressure: 2.5 GPa") == 2.5
    assert extract_number_from_string("No numbers here", 0.0) == 0.0
    assert extract_number_from_string("1.23e-4") == 1.23e-4


def test_stream_to_list():
    """Test stream to list conversion"""
    test_content = "Line 1\nLine 2\nLine 3\n"
    stream = io.StringIO(test_content)

    result = stream_to_list(stream)
    assert result == ["Line 1", "Line 2", "Line 3"]


def test_format_time_elapsed():
    """Test time formatting"""
    assert format_time_elapsed(30.5) == "30.5s"
    assert format_time_elapsed(150.0) == "2m 30s"
    assert format_time_elapsed(3661.0) == "1h 1m 1s"


def test_find_pattern_in_file(tmpfile):
    """Test finding patterns in files"""
    test_content = """This is a test file
With multiple lines
Some lines contain ENERGY: -10.5
Others have FORCE: 0.25
Final ENERGY: -8.2 eV
"""

    with open(tmpfile, 'w') as f:
        f.write(test_content)

    # Find energy lines
    energy_lines = find_pattern_in_file(tmpfile, r"ENERGY:")
    assert len(energy_lines) == 2
    assert "ENERGY: -10.5" in energy_lines[0]
    assert "ENERGY: -8.2" in energy_lines[1]

    # Test max_matches
    energy_lines_limited = find_pattern_in_file(tmpfile, r"ENERGY:", max_matches=1)
    assert len(energy_lines_limited) == 1


def test_count_pattern_in_file(tmpfile):
    """Test counting patterns in files"""
    test_content = """This is a test file
With multiple lines
Some lines contain ENERGY
Others have different content
Final ENERGY line
Another ENERGY mention
"""

    with open(tmpfile, 'w') as f:
        f.write(test_content)

    count = count_pattern_in_file(tmpfile, r"ENERGY")
    assert count == 3

    # Test with non-existent pattern
    count_none = count_pattern_in_file(tmpfile, r"NONEXISTENT")
    assert count_none == 0


def test_trim_stream_with_extra_remove():
    """Test trim_stream with extra patterns to remove"""
    test_content = """Start of content
Keep this line
FIX_VOL line should be removed
ANG line should be removed
Keep this line too
%ENDBLOCK LATTICE"""

    stream = io.StringIO(test_content)
    trimmed = trim_stream(
        stream,
        r"^Keep",
        r"%ENDBLOCK",
        extra_remove=["FIX_VOL", "ANG"]
    )

    result = trimmed.read()
    assert "Keep this line" in result
    assert "Keep this line too" in result
    assert "FIX_VOL line should be removed" not in result
    assert "ANG line should be removed" not in result

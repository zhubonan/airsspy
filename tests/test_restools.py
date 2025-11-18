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
Tests for the restool module
"""

import pytest
from airsspy.restools import (
    RESFile,
    TitlInfo,
    extract_res,
    format_minsep,
    get_minsep,
    parse_titl,
    read_res_atoms,
    save_airss_res,
    unique,
)
from tempfile import mkstemp
import numpy as np


res_example = """
TITL Al-574-3531-1 -0.05 60.5091828526 -53.2326560919 0 0 8 (Fm-3m) n - 1
REM
REM Run started: at Thu May 16 21:47:56 BST 2019 in /home/bonan/appdir/airss-git/examples/1.1
REM AIRSS Version 0.9.1
REM Al.pp (a228ffcbdc33ac12a8acef1e350d2b41)
REM
CELL 1.54180    3.81630    3.81630    4.92683  104.96324  104.96324  109.47112
LATT -1
SFAC Al
Al     1  0.7051152169851  0.6124565323292  0.0296981852527 1.0
Al     1  0.8301152168366  0.2374565298261  0.2796981816833 1.0
Al     1  0.2051152186055  0.1124565295978  0.0296981856128 1.0
Al     1  0.3301152106585  0.7374565288889  0.2796981811681 1.0
Al     1  0.9551152048079  0.8624565202265  0.5296981643404 1.0
Al     1  0.5801152067938  0.9874565210782  0.7796981682082 1.0
Al     1  0.4551152132493  0.3624565200782  0.5296981654291 1.0
Al     1  0.0801152120634  0.4874565179751  0.7796981683054 1.0
END
"""


def test_extract_res(tmpfile):
    """Test extract infomation from res file"""
    tmpfile = mkstemp()[1]
    with open(tmpfile, "w") as fh:
        fh.write(res_example)
    data = extract_res(tmpfile)
    assert "rem" in data
    assert data["uid"] == "Al-574-3531-1"
    assert data["P"] == -0.05
    assert data["V"] == 60.5091828526
    assert data["H"] == -53.2326560919
    assert data["nat"] == 8
    assert data["sym"] == "(Fm-3m)"


def test_save_airss_res(al_atoms, tmpfile):
    """Test saving AIRSS style SHELX file and extract out"""

    infodict = {
        "uid": "bla",
        "P": 1.2,
        "V": 20.0,
        "H": 1010.0,
        "nat": 10,
        "sym": "(Pm-3m)",
    }

    save_airss_res(al_atoms, infodict, fname=tmpfile, force_write=True)
    extracted = extract_res(tmpfile)
    for k, v in infodict.items():
        assert v == extracted[k]


def test_parse_titl():
    """Test parsing TITL line"""
    line = "TITL Al-574-3531-1 -0.05 60.5091828526 -53.2326560919 0 0 8 (Fm-3m) n - 1"
    titl = parse_titl(line)

    assert isinstance(titl, TitlInfo)
    assert titl.label == "Al-574-3531-1"
    assert titl.pressure == -0.05
    assert titl.volume == 60.5091828526
    assert titl.enthalpy == -53.2326560919
    assert titl.spin == 0.0
    assert titl.spin_abs == 0.0
    assert titl.natoms == 8
    assert titl.symm == "(Fm-3m)"
    assert titl.flag1 == "n"
    assert titl.flag2 == "-"
    assert titl.flag3 == "1"


def test_read_res_atoms():
    """Test reading RES file as Atoms"""
    lines = res_example.strip().split('\n')
    titl, atoms = read_res_atoms(lines)

    assert isinstance(titl, TitlInfo)
    assert titl.natoms == 8
    assert len(atoms) == 8
    assert atoms.get_chemical_formula() == "Al8"


def test_unique():
    """Test unique function"""
    items = ["Al", "Si", "Al", "C", "Si", "Al"]
    result = unique(items)
    assert result == ["Al", "Si", "C"]


def test_get_minsep():
    """Test minimum separation calculation"""
    species = ["Al", "Si", "Al"]
    # Simple distance matrix for testing
    dist_matrix = np.array([
        [0.0, 2.5, 3.0],
        [2.5, 0.0, 2.8],
        [3.0, 2.8, 0.0]
    ])

    minsep = get_minsep(species, dist_matrix)
    expected = {
        "Al-Al": 3.0,  # distance between the two Al atoms
        "Al-Si": 2.5,  # minimum Al-Si distance
    }

    assert minsep == expected


def test_format_minsep():
    """Test formatting minimum separations"""
    minsep = {"Al-Al": 2.5, "Al-Si": 3.1}
    result = format_minsep(minsep)
    assert "Al-Al=2.50" in result
    assert "Al-Si=3.10" in result


def test_resfile_from_string():
    """Test RESFile creation from string"""
    res_obj = RESFile.from_string(res_example)

    assert res_obj.label == "Al-574-3531-1"
    assert res_obj.natoms == 8
    assert res_obj.pressure == -0.05
    assert res_obj.volume == 60.5091828526
    assert res_obj.enthalpy == -53.2326560919


def test_resfile_from_file(tmpfile):
    """Test RESFile creation from file"""
    # Write test RES file
    with open(tmpfile, "w") as f:
        f.write(res_example)

    res_obj = RESFile.from_file(tmpfile)

    assert res_obj.label == "Al-574-3531-1"
    assert res_obj.natoms == 8


def test_resfile_properties():
    """Test RESFile properties"""
    res_obj = RESFile.from_string(res_example)

    assert res_obj.formula == "Al8"
    assert res_obj.reduced_formula == "Al"
    assert res_obj.spin == 0.0
    assert res_obj.spin_abs == 0.0


def test_resfile_to_res_lines():
    """Test RESFile to RES lines conversion"""
    res_obj = RESFile.from_string(res_example)
    lines = res_obj.to_res_lines()

    assert isinstance(lines, list)
    assert any("TITL" in line for line in lines)
    assert any("CELL" in line for line in lines)
    assert any("END" in line for line in lines)

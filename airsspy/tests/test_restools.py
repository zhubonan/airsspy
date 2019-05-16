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

from ..restools import extract_res, save_airss_res
from tempfile import mkstemp


res_example = \
"""
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
    with open(tmpfile, 'w') as fh:
        fh.write(res_example)
    data = extract_res(tmpfile)
    assert 'rem' in data
    assert data['uid'] == 'Al-574-3531-1'
    assert data['P'] == -0.05
    assert data['V'] == 60.5091828526
    assert data['H'] == -53.2326560919
    assert data['nat'] == 8
    assert data['sym'] == '(Fm-3m)'


def test_save_airss_res(al_atoms, tmpfile):
    """Test saving AIRSS style SHELX file and extract out"""

    infodict = {
        'uid': 'bla',
        'P': 1.2,
        'V': 20.,
        'H': 1010.,
        'nat': 10,
        'sym': '(Pm-3m)'
    }

    save_airss_res(al_atoms, infodict, fname=tmpfile, force_write=True)
    extracted = extract_res(tmpfile)
    for k, v in infodict.items():
        assert v == extracted[k]

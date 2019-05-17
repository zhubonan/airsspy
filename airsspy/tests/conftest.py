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
Test configuration
"""
from ase import Atoms
from tempfile import mkstemp
import os
import pytest


@pytest.fixture
def al_atoms():
    return Atoms(
        'Al2',
        cell=[2, 2, 2],
        positions=[[0., 0., 0.], [1., 1., 1.]],
        pbc=True)


@pytest.fixture
def tmpfile():
    fname = mkstemp()[1]
    yield fname
    os.remove(fname)

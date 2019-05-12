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
Test the Buildcell class
"""

from __future__ import absolute_import
import pytest
from distutils import spawn
from ..seed import TemplateAtoms
from ..build import Buildcell


@pytest.fixture
def template_c2():
    c2 = TemplateAtoms('C2')
    c2.build.slack = 1
    c2.build.overlap = 1
    c2.build.minsep = 1.5
    c2.build.symmops = (2, 4)
    c2.build.varvol = 20
    return c2


@pytest.mark.skipif(
    spawn.find_executable('buildcell') is None,
    reason='No buildcell executable in PATH')
def test_generate(template_c2):

    bc = Buildcell(template_c2)
    atoms = bc.generate()
    assert atoms
    assert bc.bc_err
    assert bc.bc_out

    # The method of the template atoms should also work
    atoms = template_c2.get_random_atoms()
    assert atoms

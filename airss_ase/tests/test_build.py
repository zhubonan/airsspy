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


@pytest.mark.skipif(spawn.find_executable('buildcell') is None,
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

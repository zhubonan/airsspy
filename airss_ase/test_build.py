"""
Test the Buildcell class
"""

import pytest
from .seed import TemplateAtoms


@pytest.fixture
def template_c2():
    c2 = TemplateAtoms('C2')
    c2.build_param.slack = 1
    c2.build_param.overlap = 1
    c2.build_param.minsep = 1.5
    c2.build_param.symmops = (2, 4)
    c2.build_param.varvol = 20
    return c2


def test_generate(template_c2):
    from .build import Buildcell

    bc = Buildcell(template_c2)
    atoms = bc.generate()
    assert atoms
    assert bc.bc_err
    assert bc.bc_out

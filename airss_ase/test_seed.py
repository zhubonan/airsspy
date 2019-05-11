import pytest


def test_bc_param():
    from .seed import BuildcellParam

    bcp = BuildcellParam()
    bcp.fix = True
    bcp.nforms = 3
    assert 'NFORMS' in bcp.to_string()
    bcp.minsep = [2, {'Ce-O': (2, 3)}]
    assert 'Ce-O=2-3' in bcp.to_string()


def test_nested_range():

    from .seed import BuildcellParam
    bcp = BuildcellParam()

    # Test different possible configurations
    bcp.minsep = 2
    assert 'MINSEP=2' in bcp.to_string()

    bcp.minsep = (2, 3)
    assert 'MINSEP=2-3' in bcp.to_string()

    bcp.minsep = (2, {'Ce-O': 1})
    assert 'MINSEP=2 Ce-O=1' in bcp.to_string()

    bcp.minsep = ((2, 3), {'Ce-O': 1})
    assert 'MINSEP=2-3 Ce-O=1' in bcp.to_string()

    bcp.minsep = (2, {'Ce-O': (1, 2)})
    assert 'MINSEP=2 Ce-O=1-2' in bcp.to_string()


def test_atom_param():
    """
    Test the single atom specification
    """
    from .seed import SingeAtomParam

    param = SingeAtomParam()
    param.num = 3
    param.posamp = (1, 2)
    param.tagname = 'O1'
    string = param.get_append_string()
    assert string.startswith('# O1 %')
    assert 'POSAMP=1-2' in string
    assert 'NUM=3' in string


def test_template_atom():

    from .seed import TemplateAtom
    ta = TemplateAtom(symbol='C')
    ta.xamp = 0
    ta.tagname = 'C1'
    assert ta.get_append_string() == '# C1 % XAMP=0'


def test_template_atom_from_tmp():

    from .seed import TemplateAtoms
    bc = TemplateAtoms(symbols='C2')
    c1 = bc[0]
    c1.posamp = 1
    assert c1.get_append_string() == '# C0 % POSAMP=1'
    assert bc.arrays['buildtag'][0].get_append_string() == '# C0 % POSAMP=1'


def test_tuple2range():
    """
    Test tuple/number to string conversion
    """
    from .seed import tuple2range
    assert tuple2range(2) == '2'
    assert tuple2range((2, 3)) == '2-3'


def test_template_atom_write():
    """
    Test writing output from a template atom
    """
    from .seed import TemplateAtoms

    atoms = TemplateAtoms('C4')
    atoms.build_param.symmops = (2, 3)
    atoms.build_param.sgrank = 2
    atoms.build_param.minsep = [2, {'C-C': (2, 3)}]
    atoms[0].posamp = 3
    atoms[0].xamp = (2, 3)
    atoms[1].tagname = 'CX0'
    buff = '\n'.join(atoms.get_seed_lines())

    assert '#SGRANK=2' in buff
    assert '#SYMMOPS=2-3' in buff
    assert '#MINSEP=2 C-C=2-3' in buff
    assert 'XAMP=2-3' in buff
    assert '# C0' in buff
    assert '# CX0 ' in buff

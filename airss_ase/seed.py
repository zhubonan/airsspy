"""
Classes for preparing AIRSS seed
"""
import numpy as np
from ase import Atoms, Atom
import numbers
import collections
from ase.constraints import FixConstraint, FixBondLengths


class TemplateAtoms(Atoms):
    """Subclass of ase.atoms.Atoms object. Template for generating random cells
    """

    def __init__(self, *args, **kwargs):
        """Initialise an TemplateAtoms for buildcell.
        Same arguments signature as ase.Atoms object

        Special attribute:
          build_param : BuildcellParam instance for storing parameter of buildcell

        An array of SingeAtomParam objects are added automatically.
        Each one can be retrieved/replaced with
          set_buildcell_tag
          get_buildcell_tag
        """
        super(TemplateAtoms, self).__init__(*args, **kwargs)
        self.build_param = BuildcellParam()
        # Construct tags for each Atom
        tags = []
        symbols = self.get_chemical_symbols()
        for i in range(len(self)):
            tag = SingeAtomParam()
            # Make independent tags initially
            tag.tagname = symbols[i] + str(i)
            tags.append(tag)

        self.new_array('buildtag', tags, dtype=object, shape=None)

    def set_buildcell_tag(self, tag, index):
        """Set buildcell tags for individual atom
        if the SingeAtomParam object has no tagname, set automatically"""
        if tag.tagname is None:
            tag.tagname = self.get_chemical_symbols()[index]
        self.arrays['buildtag'][index] = tag

    def get_buidcell_tag(self, index):
        """
        Return the buildcell tag for the atom of the index.
        Can be used for in-place change
        """
        return self.arrays['buildtag'][index]

    def __getitem__(self, i):
        """Return a subset of the atoms.

        i -- scalar integer, list of integers, or slice object
        describing which atoms to return.

        If i is a scalar, return an Atom object. If i is a list or a
        slice, return an Atoms object with the same cell, pbc, and
        other associated info as the original Atoms object. The
        indices of the constraints will be shuffled so that they match
        the indexing in the subset returned.

        """

        if isinstance(i, numbers.Integral):
            natoms = len(self)
            if i < -natoms or i >= natoms:
                raise IndexError('Index out of range.')

            return TemplateAtom(atoms=self, index=i)
        elif isinstance(i, list) and len(i) > 0:
            # Make sure a list of booleans will work correctly and not be
            # interpreted at 0 and 1 indices.
            i = np.array(i)

        import copy

        conadd = []
        # Constraints need to be deepcopied, but only the relevant ones.
        for con in copy.deepcopy(self.constraints):
            if isinstance(con, (FixConstraint, FixBondLengths)):
                try:
                    con.index_shuffle(self, i)
                    conadd.append(con)
                except IndexError:
                    pass

        atoms = self.__class__(
            cell=self._cell,
            pbc=self._pbc,
            info=self.info,
            # should be communicated to the slice as well
            celldisp=self._celldisp)
        # TODO: Do we need to shuffle indices in adsorbate_info too?

        atoms.arrays = {}
        for name, a in self.arrays.items():
            atoms.arrays[name] = a[i].copy()

        atoms.constraints = conadd
        return atoms


# Singular, plural, default value:
bc_param_names = {
    "FIX": 'tag',
    "CFIX": 'tag',
    "NFORM": 'range',
    "SUPERCELL": 'range',
    "SLAB": 'tag',
    "CLUSTER": 'tag'
}


def tagproperty(name, doc):
    """Set a tag-like property"""

    def getter(self):
        return self.get(name)

    def setter(self, value):
        self.type_registry.update({name: 'tag'})
        if value is True:
            self.set_tag(name)
        elif value is False:
            self.delete(name)

    def deleter(self):
        self.delete(name)

    return property(getter, setter, deleter, doc)


def genericproperty(name, doc):
    """Set a range-like property"""

    def getter(self):
        return self.get_prop(name)

    def setter(self, value):
        self.type_registry.update({name: 'generic'})
        self.set_prop(name, value)

    def deleter(self):
        self.delete_prop(name)

    return property(getter, setter, deleter, doc)


def rangeproperty(name, doc):
    def getter(self):
        return self.get_prop(name)

    def setter(self, value):
        if isinstance(value, (tuple, list)):
            if len(value) != 2:
                raise ValueError("A tuple/list of two element must be used.")
            if any([not isinstance(x, numbers.Number) for x in value]):
                raise ValueError("Both elements need to be a number")
        self.type_registry.update({name: 'range'})
        self.set_prop(name, value)

    def deleter(self):
        self.delete_prop(name)

    return property(getter, setter, deleter, doc)


def nestedrangeproperty(name, doc):
    def getter(self):
        return self.get_prop(name)

    def setter(self, value):
        if isinstance(value, (tuple, list)):
            if len(value) != 2:
                raise RuntimeError("A tuple/list of two element must be used.")
        self.type_registry.update({name: 'nested_range'})
        self.set_prop(name, value)

    def deleter(self):
        self.delete_prop(name)

    return property(getter, setter, deleter, doc)


class BuildcellParam(object):
    """
    A class for storing parameters for the Buldcell program
    """

    def __init__(self):
        self.data = dict()
        self.type_registry = dict()

    def clear_all(self):
        """Set all property to be None"""
        self.data.clear()

    def get_prop(self, value):
        return self.data.__getitem__(value)

    def set_prop(self, name, value):
        return self.data.__setitem__(name, value)

    def set_tag(self, tag):
        self.data.__setitem__(tag, '')

    def delete_prop(self, name):
        self.data.pop(name)

    def to_string(self):
        """Return the string that should go into the .cell file"""
        lines = []
        for key, value in self.data.items():
            name = key.upper()
            type_string = self.type_registry[key]
            if value is False or value is None:
                continue
            if name == 'FIX':
                continue
            if type_string == 'tag':
                lines.append("#{}".format(name))
            # Allow direct passing of string
            elif type_string == 'generic':
                lines.append("#{}={}".format(name, value))
                continue
            elif type_string in ('range', 'nested_range'):
                # Check if there is a dictionary to unpack
                # The value can be
                if not isinstance(value, (list, tuple)):
                    line = '#{}={}'.format(name, value)
                else:
                    # The value is a list/tuple
                    # is a simple range?
                    if all([isinstance(x, numbers.Number) for x in value]):
                        line = '#{}={}'.format(name, tuple2range(value))
                    # Deal with a nested range
                    else:
                        line = '#{}={}'.format(name, tuple2range(value[0]))
                        tokens = [line]
                        if value[1]:
                            for k_, v_ in value[1].items():
                                tokens.append(k_ + '=' + tuple2range(v_))
                        line = ' '.join(tokens)
                lines.append(line)

        return '\n'.join(lines) + '\n'

    fix = tagproperty('FIX', 'Fix the cell')
    cfix = tagproperty('CFIX', 'Fix the caxis')
    cluster = tagproperty('CLUSTER', 'We are predicting CLUSTER')
    nforms = genericproperty('NFORMS', 'Number of formula units')
    minsep = nestedrangeproperty('MINSEP', 'Minimum separation constraints')
    posamp = nestedrangeproperty('POSAMP', 'Position amplitudes')
    symmops = rangeproperty('SYMMOPS',
                            'Number of symmetry operation requested cell')
    minamp = rangeproperty('MINAMP', 'Minimum aplitude of randomisation')
    zamp = rangeproperty('ZAMP', 'Randomisation amplitude in Z')
    xamp = rangeproperty('XAMP', 'Randomisation amplitude in X')
    yamp = rangeproperty('YAMP', 'Randomisation amplitude in Y')
    angamp = rangeproperty('ANGAMP',
                           'Angular randomisation amplitude from fragments')
    sgrank = rangeproperty('SGRANK', 'Minimum rank of the spacegroup')
    species = nestedrangeproperty('SPECIES', 'Species to be put into the cell')
    varvol = genericproperty(
        'VARVOL', 'Target volume of cell with the original configuration')
    slack = rangeproperty(
        'SLACK', 'Slack the hard sphere potentials enforcing the MINSEP')
    overlap = rangeproperty(
        'OVERLAP', 'Threhold of the overlap for the hard sphere potentials')
    #focus = genericproperty('FOCUS', 'Focus on a specific compositions')
    compact = tagproperty('COMPACT', 'Compact the cell using Niggli reduction')
    cons = genericproperty('CONS', 'Parameter for cell shape constraint')


def test_bc_param():

    bcp = BuildcellParam()
    bcp.fix = True
    bcp.nforms = 3
    assert 'NFORMS' in bcp.to_string()
    bcp.minsep = [2, {'Ce-O': (2, 3)}]
    assert 'Ce-O=2-3' in bcp.to_string()


def test_nested_range():

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


class SingeAtomParam(object):
    """Paramter for a single auto"""

    def __init__(self, *args, **kwargs):
        self.prop_data = dict()
        self.type_registry = dict()
        self.disabled = False

    def clear_all(self):
        """Set all property to be None"""
        self.prop_data.clear()

    def get_prop(self, value):
        return self.prop_data.__getitem__(value)

    def set_prop(self, name, value):
        return self.prop_data.__setitem__(name, value)

    def set_tag(self, tag):
        self.prop_data.__setitem__(tag, '')

    def delete_prop(self, name):
        self.prop_data.pop(name)

    tagname = genericproperty('tagname', 'Name of the tag')
    posamp = rangeproperty('POSAMP', 'Position amplitude')
    minamp = rangeproperty('POSAMP', 'Minimum positional amplitude')
    zamp = rangeproperty('ZAMP', 'Amplitude in Z')
    xamp = rangeproperty('XAMP', 'Amplitude in X')
    yamp = rangeproperty('YAMP', 'Amplitude in Y')
    num = rangeproperty('NUM', 'Number of atoms/fragments')
    atatom = rangeproperty('ADATOM', 'Add atoms after making supercell')
    fix = tagproperty('FIX', 'FIX this atom')

    def get_append_string(self):
        """
        Return the per entry string for this atom to be appended after its
        line in POSITIONS_FRAC / POSITIONS_ABS block
        """
        if self.disabled is True:
            return ''

        tokens = []
        # Set the tag
        tagname = self.prop_data.get('tagname')
        if not tagname:
            raise ValueError('The tagname property must be set')

        tokens.append('# {} %'.format(self.prop_data['tagname']))
        for key, value in self.prop_data.items():
            if key == "tagname":
                continue

            type_string = self.type_registry[key]
            name = key.upper()

            # Process the value based on type string
            if type_string == 'tag':
                tokens.append(name)
            elif type_string == 'range':
                if isinstance(value, (list, tuple)):
                    tokens.append('{}={}-{}'.format(name, *value))
                else:
                    tokens.append('{}={}'.format(name, value))
            else:
                tokens.append("{}={}".format(name, value))

        string = ' '.join(tokens)
        return string


def test_atom_param():
    """
    Test the single atom specification
    """
    param = SingeAtomParam()
    param.num = 3
    param.posamp = (1, 2)
    param.tagname = 'O1'
    string = param.get_append_string()
    assert string.startswith('# O1 %')
    assert 'POSAMP=1-2' in string
    assert 'NUM=3' in string


class TemplateAtom(Atom, SingeAtomParam):
    """
    Element atoms in a AIRSS seed
    """

    def __init__(self, *args, **kwargs):
        super(TemplateAtom, self).__init__(*args, **kwargs)
        SingeAtomParam.__init__(self, *args, **kwargs)
        if self.atoms is not None:
            self.prop_data = self.atoms.arrays['buildtag'][
                self.index].prop_data
            self.type_registry = self.atoms.arrays['buildtag'][
                self.index].type_registry


def test_template_atom():

    ta = TemplateAtom(symbol='C')
    ta.xamp = 0
    ta.tagname = 'C1'
    assert ta.get_append_string() == '# C1 % XAMP=0'


def test_template_atom_from_tmp():
    bc = TemplateAtoms(symbols='C2')
    c1 = bc[0]
    c1.posamp = 1
    assert c1.get_append_string() == '# C0 % POSAMP=1'
    assert bc.arrays['buildtag'][0].get_append_string() == '# C0 % POSAMP=1'


def tuple2range(value):
    """
    Return the string for a given value. If the value is a tuple
    make it a range.
    """
    if isinstance(value, (list, tuple)):
        return "{}-{}".format(value[0], value[1])
    else:
        return str(value)


# The following function is modified base on ase.io.castep.write_castep_cell
def write_buildcell_seed(fd,
                         atoms,
                         positions_frac=False,
                         castep_cell=None,
                         force_write=False):
    """
    This CASTEP export function write minimal information to
    a .cell file. If the atoms object is a trajectory, it will
    take the last image.

    Note that function has been altered in order to require a filedescriptor
    rather than a filename. This allows to use the more generic write()
    function from formats.py

    Note that the "force_write" keywords has no effect currently.
    """
    import numpy as np
    from ase.constraints import FixedLine, FixAtoms, FixCartesian
    if atoms is None:
        print('Atoms object not initialized')
        return False
    if isinstance(atoms, list):
        if len(atoms) > 1:
            atoms = atoms[-1]

# deprecated; should be handled on the more generic write() level
#    if os.path.isfile(filename) and not force_write:
#        print('ase.io.castep.write_param: Set optional argument')
#        print('force_write=True to overwrite %s.' % filename)
#        return False

#    fd = open(filename, 'w')

# I have to disable extra # comments
# fd.write('#######################################################\n')
# fd.write('#CASTEP cell file: %s\n' % fd.name)
# fd.write('#Created using the Atomic Simulation Environment (ASE)#\n')
# fd.write('#######################################################\n\n')
    fd.write('%BLOCK LATTICE_CART\n')
    cell = np.matrix(atoms.get_cell())
    for line in atoms.get_cell():
        fd.write('    %.10f %.10f %.10f\n' % tuple(line))
    # NOTE Addition for buildcell tag
    if atoms.build_param.fix is True:
        fd.write('#FIX\n')
    if atoms.build_param.cfix is True:
        fd.write('#CFIX\n')
    fd.write('%ENDBLOCK LATTICE_CART\n\n\n')

    if positions_frac:
        keyword = 'POSITIONS_FRAC'
        positions = np.array(atoms.get_positions() * cell.I)

    else:
        keyword = 'POSITIONS_ABS'
        positions = atoms.get_positions()

    if (hasattr(atoms, 'calc') and hasattr(atoms.calc, 'param')
            and hasattr(atoms.calc.param, 'task')):
        _spin_pol = any([
            getattr(atoms.calc.param, i).value
            for i in ['spin_polarized', 'spin_polarised']
        ])
    else:
        _spin_pol = True

    # Gather the data that will be used to generate the block
    pos_block_data = []
    pos_block_format = '%s %8.6f %8.6f %8.6f'
    if atoms.has('castep_custom_species'):
        pos_block_data.append(atoms.get_array('castep_custom_species'))
    else:
        pos_block_data.append(atoms.get_chemical_symbols())
    pos_block_data += [xlist for xlist in zip(*positions)]

    if atoms.get_initial_magnetic_moments().any() and _spin_pol:
        pos_block_data.append(atoms.get_initial_magnetic_moments())
        pos_block_format += ' SPIN=%4.2f'

    pos_block = [(pos_block_format % line_data)
                 for line_data in zip(*pos_block_data)]

    # Adding the CASTEP labels output
    if atoms.has('castep_labels'):
        labels = atoms.get_array('castep_labels')
        for l_i, label in enumerate(labels):
            # avoid empty labels that crash CASTEP runs
            if label and label != 'NULL':
                pos_block[l_i] += ' LABEL=%s' % label

    # NOTE Adding buildcell tags
    build_tags = [tag.to_string() for tag in atoms.get_array('buildtag')]
    for l_i, tag in enumerate(build_tags):
        pos_block[l_i] += tag

    fd.write('%%BLOCK %s\n' % keyword)
    for line in pos_block:
        fd.write('    %s\n' % line)
    fd.write('%%ENDBLOCK %s\n\n' % keyword)

    # NOTE Adding BuildcellParam
    fd.write(atoms.build_param.to_string())

    # if atoms, has a CASTEP calculator attached, then only
    # write constraints if really necessary
    if (hasattr(atoms, 'calc') and hasattr(atoms.calc, 'param')
            and hasattr(atoms.calc.param, 'task')):
        task = atoms.calc.param.task
        if atoms.calc.param.task.value is None:
            suppress_constraints = True
        elif task.value.lower() not in [
                'geometryoptimization',
                # well, CASTEP understands US and UK english...
                'geometryoptimisation',
                'moleculardynamics',
                'transitionstatesearch',
                'phonon'
        ]:
            suppress_constraints = True
        else:
            suppress_constraints = False
    else:
        suppress_constraints = True

    constraints = atoms.constraints
    if len(constraints) and not suppress_constraints:
        fd.write('%BLOCK IONIC_CONSTRAINTS \n')
        count = 0
        for constr in constraints:
            if (not isinstance(constr, FixAtoms)
                    and not isinstance(constr, FixCartesian)
                    and not isinstance(constr, FixedLine)
                    and not suppress_constraints):
                print('Warning: you have constraints in your atoms, that are')
                print('         not supported by the CASTEP ase interface')
                break
            if isinstance(constr, FixAtoms):
                # sorry, for this complicated block
                # reason is that constraint.index can either
                # hold booleans or integers and in both cases
                # it is an numpy array, so no simple comparison works
                for n, val in enumerate(constr.index):
                    if val.dtype.name.startswith('bool'):
                        if not val:
                            continue
                        symbol = atoms.get_chemical_symbols()[n]
                        nis = atoms.calc._get_number_in_species(n)
                    elif val.dtype.name.startswith('int'):
                        symbol = atoms.get_chemical_symbols()[val]
                        nis = atoms.calc._get_number_in_species(val)
                    else:
                        raise UserWarning('Unrecognized index in' +
                                          ' constraint %s' % constr)
                    fd.write('%6d %3s %3d   1 0 0 \n' %
                             (count + 1, symbol, nis))
                    fd.write('%6d %3s %3d   0 1 0 \n' %
                             (count + 2, symbol, nis))
                    fd.write('%6d %3s %3d   0 0 1 \n' %
                             (count + 3, symbol, nis))
                    count += 3
            elif isinstance(constr, FixCartesian):
                n = constr.a
                symbol = atoms.get_chemical_symbols()[n]
                nis = atoms.calc._get_number_in_species(n)
                # fix_cart = - constr.mask + 1
                # just use the logical opposite
                fix_cart = np.logical_not(constr.mask)
                if fix_cart[0]:
                    count += 1
                    fd.write('%6d %3s %3d   1 0 0 \n' % (count, symbol, nis))
                if fix_cart[1]:
                    count += 1
                    fd.write('%6d %3s %3d   0 1 0 \n' % (count, symbol, nis))
                if fix_cart[2]:
                    count += 1
                    fd.write('%6d %3s %3d   0 0 1 \n' % (count, symbol, nis))
            elif isinstance(constr, FixedLine):
                n = constr.a
                symbol = atoms.get_chemical_symbols()[n]
                nis = atoms.calc._get_number_in_species(n)
                direction = constr.dir
                # print(direction)
                ((i1, v1), (i2, v2)) = sorted(enumerate(direction),
                                              key=lambda x: abs(x[1]),
                                              reverse=True)[:2]
                # print(sorted(enumerate(direction), key = lambda x:x[1])[:2])
                # print(sorted(enumerate(direction), key = lambda x:x[1]))

                # print(v1)
                # print(v2)
                n1 = np.array([v2, v1, 0])
                n1 = n1 / np.linalg.norm(n1)

                n2 = np.cross(direction, n1)
                count += 1
                fd.write('%6d %3s %3d   %f %f %f \n' %
                         (count, symbol, nis, n1[0], n1[1], n1[2]))

                count += 1
                fd.write('%6d %3s %3d   %f %f %f \n' %
                         (count, symbol, nis, n2[0], n2[1], n2[2]))
        fd.write('%ENDBLOCK IONIC_CONSTRAINTS \n')

    if castep_cell is None:
        if hasattr(atoms, 'calc') and hasattr(atoms.calc, 'cell'):
            castep_cell = atoms.calc.cell
        else:
            # fd.close()
            return True

    for option in castep_cell._options.values():
        if option.value is not None:
            #            print(option.value)
            if option.type == 'Block':
                fd.write('%%BLOCK %s\n' % option.keyword.upper())
                fd.write(option.value)
                fd.write('\n%%ENDBLOCK %s\n\n' % option.keyword.upper())
            else:
                fd.write('%s : %s\n\n' %
                         (option.keyword.upper(), option.value))


#    fd.close()
    return True


class Buildcell:
    def __init__(self, atoms):
        """Initialise an Buildcell object"""
        self.atoms = atoms
        self.proc = None
        self.res = None
        self.bc_out = None
        self.bc_err = None

    def generate(self, timeout=10, write_cell=None):
        """Generate a random atom based on a template
        timeout: time to wait for buildcell binary
        write_seed : Name of the output cell to be written"""
        import io
        import ase.io.castep
        import subprocess as sbp
        tmp = io.StringIO()
        write_buildcell_seed(tmp, self.atoms)
        bc = sbp.Popen('buildcell',
                       universal_newlines=True,
                       stdin=sbp.PIPE,
                       stdout=sbp.PIPE,
                       stderr=sbp.PIPE)
        self.proc = bc
        tmp.seek(0)
        self.seed = tmp.read()
        try:
            self.bc_out, self.bc_err = bc.communicate(input=self.seed,
                                                      timeout=timeout)
        except sbp.TimeoutExpired:
            bc.kill()
            self.bc_out, self.bc_err = bc.communicate()
            print('Generation Failed to finished. Output captured')
            return

        tmp.close()
        tmp = io.StringIO()
        tmp.write(self.bc_out)
        tmp.seek(0)
        res = ase.io.castep.read_castep_cell(tmp)

        # Write the output from buildcell
        if write_cell:
            with open(write_cell + '.cell', 'w') as output:
                tmp.seek(0)
                for line in tmp:
                    output.write(line)
        tmp.close()

        # Detach unnecessary calculator
        res.calc = None
        # Store as attribute
        self.res = res
        return res

    def write_seed(self, seedname):
        """Write the seed for buildcell to the disk"""
        with open(seedname + '.cell', 'w') as tmp:
            write_buildcell_seed(tmp, self.atoms)

    def gen_and_view(self, viewer=None, wrap=False, timeout=20):
        from ase.visualize import view
        res = self.generate(timeout=timeout)
        if not res:
            return
        if wrap:
            res.wrap()
        if viewer:
            view(res, viewer=viewer)
        else:
            view(res)

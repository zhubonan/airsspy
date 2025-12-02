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
Classes for preparing AIRSS seed
"""
import numbers
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from ase import Atom, Atoms
from ase.atoms import Atoms as ASEAtoms
from ase.constraints import FixBondLengths, FixConstraint
from castepinput import CellInput

from .build import BuildcellError


class SeedAtoms(Atoms):
    """Subclass of ase.atoms.Atoms object. Template for generating random cells"""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        """Initialise an SeedAtoms for buildcell.
        Same arguments signature as ase.Atoms object

        Special attribute:
          gentags : BuildcellParam instance for storing parameter of buildcell

        You can set the global buildcell tags by setting attributes of the gentags
        instance.

        Per-atom generation tags are saved in an array called 'atom_gentags'.
        They can be retrieved using ``get_atom_tag`` or set with ``set_atom_tag``.
        The ``atom_tags`` attribute also provide a convient way of accessing.
        """
        super().__init__(*args, **kwargs)
        self.gentags: BuildcellParam = BuildcellParam()
        # Construct tags for each Atom
        tags: List[SeedAtomTag] = []
        symbols = self.get_chemical_symbols()
        for i in range(len(self)):
            tag = SeedAtomTag()
            # Make independent tags initially
            tag.tagname = symbols[i] + str(i)
            tags.append(tag)

        self.new_array("atom_gentags", tags, dtype=object, shape=None)

    def set_atom_tag(self, tag: "SeedAtomTag", index: int) -> None:
        """Set buildcell tags for individual atom
        if the SeedAtomTag object has no tagname, set automatically"""
        if tag.tagname is None:
            tag.tagname = self.get_chemical_symbols()[index]
        self.arrays["atom_gentags"][index] = tag

    def get_atom_tag(self, index: int) -> Any:
        """
        Return the buildcell tag for the atom of the index.
        Can be used for in-place change
        """
        return self.arrays["atom_gentags"][index]  # type: ignore[return-value]

    @property
    def atom_tags(self) -> np.ndarray:
        """Array of tags, each for one Atom"""
        return self.arrays["atom_gentags"]

    def write_seed(self, fpath: str) -> None:
        """Write the seed to file"""
        with open(fpath, "w") as fhandle:
            fhandle.write("\n".join(self.get_cell_inp_lines()))

    def get_cell_inp(self) -> CellInput:
        """Return the python object represent the cell"""
        return get_cell_inp(self)

    def get_cell_inp_lines(self) -> List[str]:
        """
        Return a list of strings of the seed file
        """
        return get_cell_inp_lines(self)

    def build_random_atoms(
        self, timeout: int = 10, also_buildcell: bool = False, fail_ok: bool = True
    ) -> Optional[Union[ASEAtoms, Tuple[ASEAtoms, Any]]]:
        """
        Returns the randomize Atoms built using ``buildcell`` program
        """
        from .build import Buildcell

        buildcell = Buildcell(self)
        try:
            rand_atoms = buildcell.generate(timeout)
        except BuildcellError as e:
            if fail_ok:
                return None
            raise e
        if also_buildcell:
            return rand_atoms, buildcell
        return rand_atoms

    def __getitem__(
        self, i: Union[int, Sequence[int], slice]
    ) -> Union["SeedAtom", "SeedAtoms"]:  # type: ignore[override]
        """Return a subset of the atoms.

        i -- scalar integer, list of integers, or slice object
        describing which atoms to return.

        If i is a scalar, return an Atom object. If i is a list or a
        slice, return an Atoms object with the same cell, pbc, and
        other associated info as the original Atoms object. The
        indices of the constraints will be shuffled so that they match
        the indexing in the subset returned.

        """

        if isinstance(i, int):
            natoms = len(self)
            if i < -natoms or i >= natoms:
                raise IndexError("Index out of range.")

            return SeedAtom(atoms=self, index=i)
        elif isinstance(i, list) and i:
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
            cell=self.get_cell(),
            pbc=self.get_pbc(),
            info=self.info,
            # should be communicated to the slice as well
            celldisp=self.get_celldisp(),
        )
        # TODO: Do we need to shuffle indices in adsorbate_info too?

        atoms.arrays = {}
        for name, array in self.arrays.items():
            atoms.arrays[name] = array[i].copy()

        atoms.constraints = conadd
        return atoms


class BoolTag:
    """Descriptor for tag-like properties (boolean flags)"""

    def __init__(self, doc: str = "", storage_name: Optional[str] = None) -> None:
        self.doc = doc
        self.name: str = ""
        self._storage_name = storage_name

    def __set_name__(self, owner: type, name: str) -> None:
        self.name = name
        if self._storage_name is None:
            self._storage_name = name.upper()

    @property
    def storage_name(self) -> str:
        """Storage name for the property, guaranteed to be set after __set_name__"""
        assert self._storage_name is not None, "Descriptor not initialized"
        return self._storage_name

    def __get__(self, instance: Optional["TagHolder"], owner: type) -> Optional[bool]:
        if instance is None:
            return None
        return instance.get_tag(self.storage_name)

    def __set__(self, instance: "TagHolder", value: bool) -> None:
        instance.type_registry.update({self.storage_name: "tag"})
        if value is True:
            instance.set_tag(self.storage_name)
        elif value is False:
            instance.delete(self.storage_name)

    def __delete__(self, instance: "TagHolder") -> None:
        instance.delete(self.storage_name)


class GenericTag:
    """Descriptor for generic properties (any value type)"""

    def __init__(self, doc: str = "", storage_name: Optional[str] = None) -> None:
        self.doc = doc
        self.name: str = ""
        self._storage_name = storage_name

    def __set_name__(self, owner: type, name: str) -> None:
        self.name = name
        if self._storage_name is None:
            self._storage_name = name.upper()

    @property
    def storage_name(self) -> str:
        """Storage name for the property, guaranteed to be set after __set_name__"""
        assert self._storage_name is not None, "Descriptor not initialized"
        return self._storage_name

    def __get__(self, instance: Optional["TagHolder"], owner: type) -> Any:
        if instance is None:
            return None
        return instance.get_prop(self.storage_name)

    def __set__(self, instance: "TagHolder", value: Any) -> None:
        instance.type_registry.update({self.storage_name: "generic"})
        instance.set_prop(self.storage_name, value)

    def __delete__(self, instance: "TagHolder") -> None:
        instance.delete_prop(self.storage_name)


class RangeTag:
    """Descriptor for range properties (number or tuple of two numbers)"""

    def __init__(self, doc: str = "", storage_name: Optional[str] = None) -> None:
        self.doc = doc
        self.name: str = ""
        self._storage_name = storage_name

    def __set_name__(self, owner: type, name: str) -> None:
        self.name = name
        if self._storage_name is None:
            self._storage_name = name.upper()

    @property
    def storage_name(self) -> str:
        """Storage name for the property, guaranteed to be set after __set_name__"""
        assert self._storage_name is not None, "Descriptor not initialized"
        return self._storage_name

    def __get__(self, instance: Optional["TagHolder"], owner: type) -> Any:
        if instance is None:
            return None
        return instance.get_prop(self.storage_name)

    def __set__(
        self,
        instance: "TagHolder",
        value: Union[
            numbers.Number, Tuple[numbers.Number, numbers.Number], List[numbers.Number]
        ],
    ) -> None:
        if isinstance(value, (tuple, list)):
            if len(value) != 2:
                raise ValueError("A tuple/list of two element must be used.")
            if any(not isinstance(x, numbers.Number) for x in value):
                raise ValueError("Both elements need to be a number")
        instance.type_registry.update({self.storage_name: "range"})
        instance.set_prop(self.storage_name, value)

    def __delete__(self, instance: "TagHolder") -> None:
        instance.delete_prop(self.storage_name)


class NestedRangeTag:
    """Descriptor for nested range properties (tuple/list of two elements)"""

    def __init__(self, doc: str = "", storage_name: Optional[str] = None) -> None:
        self.doc = doc
        self.name: str = ""
        self._storage_name = storage_name

    def __set_name__(self, owner: type, name: str) -> None:
        self.name = name
        if self._storage_name is None:
            self._storage_name = name.upper()

    @property
    def storage_name(self) -> str:
        """Storage name for the property, guaranteed to be set after __set_name__"""
        assert self._storage_name is not None, "Descriptor not initialized"
        return self._storage_name

    def __get__(self, instance: Optional["TagHolder"], owner: type) -> Any:
        if instance is None:
            return None
        return instance.get_prop(self.storage_name)

    def __set__(
        self, instance: "TagHolder", value: Union[Tuple[Any, Any], List[Any]]
    ) -> None:
        if isinstance(value, (tuple, list)):
            if len(value) != 2:
                raise RuntimeError("A tuple/list of two element must be used.")
        instance.type_registry.update({self.storage_name: "nested_range"})
        instance.set_prop(self.storage_name, value)

    def __delete__(self, instance: "TagHolder") -> None:
        instance.delete_prop(self.storage_name)


class TagHolder:
    """Container for the tags"""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        """A container for tags of a single SeedAtom"""
        self.prop_data: Dict[str, Any] = {}
        self.type_registry: Dict[str, str] = {}
        self.disabled: bool = False

    def get_prop_dict(self) -> Dict[str, Any]:
        return self.prop_data

    def clear_all(self) -> None:
        """Set all property to be None"""
        self.prop_data.clear()

    def get_prop(self, value: str) -> Any:
        """Get property"""
        return self.prop_data.get(value)

    def set_prop(self, name: str, value: Any) -> None:
        """Set property"""
        return self.prop_data.__setitem__(name, value)

    def set_tag(self, tag: str) -> None:
        """Set a tag-like property"""
        self.prop_data.__setitem__(tag, "")

    def get_tag(self, tag: str) -> Optional[bool]:
        """Set a tag-like property"""
        value = self.prop_data.get(tag)
        if value == "":
            return True
        return None

    def delete_prop(self, name: str) -> None:
        """Deleta a property"""
        self.prop_data.pop(name)

    def delete(self, name: str) -> None:
        """Delete a property by name"""
        self.delete_prop(name)

    def to_string(self) -> str:
        raise NotImplementedError

    def __repr__(self) -> str:
        string = self.to_string()
        string = string.replace("\n", " ")
        string = string.replace("#", "")
        string = string.strip()
        if len(string) > 60:
            string = string[:60] + "..."
        return f"{type(self).__name__}<{string}>"


class BuildcellParam(TagHolder):
    """
    A class for storing parameters for the Buldcell program
    """

    def to_string(self) -> str:
        """Return the string that should go into the .cell file"""
        lines = []
        for key, value in self.prop_data.items():
            name = key.upper()
            type_string = self.type_registry[key]
            if value is False or value is None:
                continue
            if name == "FIX":
                continue
            if type_string == "tag":
                lines.append(f"#{name}")
            # Allow direct passing of string
            elif type_string == "generic":
                lines.append(f"#{name}={value}")
                continue
            elif type_string in ("range", "nested_range"):
                # Check if there is a dictionary to unpack
                # The value can be
                if not isinstance(value, (list, tuple)):
                    line = f"#{name}={value}"
                else:
                    # The value is a list/tuple
                    # is a simple range?
                    if all(isinstance(x, numbers.Number) for x in value):
                        line = f"#{name}={tuple2range(value)}"
                    # Deal with a nested range
                    else:
                        line = f"#{name}={tuple2range(value[0])}"
                        tokens = [line]
                        if value[1]:
                            for k_tmp, v_tmp in value[1].items():
                                tokens.append(k_tmp + "=" + tuple2range(v_tmp))
                        line = " ".join(tokens)
                lines.append(line)

        return "\n".join(lines) + "\n"

    fix = BoolTag("Fix the cell")
    abfix = BoolTag("Fix ab axes")
    adjgen = GenericTag("Adjust the general positions")
    autoslack = BoolTag("")
    breakamp = GenericTag("Amplitude for breaking symmetry")
    celladapt = GenericTag("")
    cellamp = GenericTag("Amplitude for cell")
    cellcon = GenericTag("")
    coord = RangeTag("")
    cylinder = GenericTag(
        "Confining cylinder(positive) or attractive line potential (neagtive)"
    )
    flip = BoolTag("Enable mirror reflection of fragments", storage_name="flip")
    maxbangle = GenericTag("")
    maxtime = GenericTag("")
    minbangle = GenericTag("")
    focus = GenericTag("Focus on composition?")
    molecules = GenericTag("")
    nocompact = GenericTag("No compact cell")
    nopush = GenericTag("No pushing")
    octet = GenericTag("")
    permfrac = GenericTag("")
    permute = GenericTag("Enable permutation or specify species")
    rad = GenericTag("")
    rash = GenericTag("")
    rash_angamp = GenericTag("")
    rash_posamp = GenericTag("")
    remove = GenericTag("")
    slab = GenericTag("")
    species = GenericTag("")
    sphere = GenericTag("")
    spin = GenericTag("")
    supercell = GenericTag("")
    surface = BoolTag("")
    symm = GenericTag("")
    symmno = GenericTag("")
    symmorphic = BoolTag("")
    system = GenericTag("Crystal system")
    targvol = RangeTag("Target volume")
    three = GenericTag("User three body hard sphere potential")
    tight = BoolTag("Tigh packing?")
    vacancies = GenericTag("Introduct vacancies")
    vacuum = GenericTag("Add vacuum")
    width = GenericTag("Width of a confining slab spacer")

    cfix = BoolTag("Fix the caxis")
    cluster = BoolTag("We are predicting CLUSTER")
    nform = RangeTag(
        "Number of formula units. "
        "This must be set otherwise the number of atoms is n times "
        "the number of symmetries"
    )
    minsep = NestedRangeTag("Minimum separation constraints")
    posamp = NestedRangeTag("Position amplitudes")
    symmops = RangeTag("Number of symmetry operation requested cell")
    minamp = RangeTag("Minimum aplitude of randomisation")
    zamp = RangeTag("Randomisation amplitude in Z")
    xamp = RangeTag("Randomisation amplitude in X")
    yamp = RangeTag("Randomisation amplitude in Y")
    angamp = RangeTag("Angular randomisation amplitude from fragments")
    sgrank = RangeTag("Minimum rank of the spacegroup")
    varvol = GenericTag("Target volume of cell with the original configuration")
    slack = RangeTag("Slack the hard sphere potentials enforcing the MINSEP")
    overlap = RangeTag("Threhold of the overlap for the hard sphere potentials")
    compact = BoolTag("Compact the cell using Niggli reduction")
    cons = GenericTag("Parameter for cell shape constraint")
    natom = RangeTag("Number of atoms in cell, if not explicitly defined")
    formula = GenericTag("Chemical formula (e.g., Si2O4)")
    seed = GenericTag("Random seed for reproducible structures")
    vol = RangeTag("Target cell volume")
    nfails = GenericTag("Number of failures before giving up")
    hole = GenericTag("Create a hole in the structure")
    holepos = GenericTag("Position for hole constraint")
    shift = GenericTag("Coordinate shift")


class SeedAtomTag(TagHolder):
    """Tags for a single auto"""

    tagname = GenericTag("Name of the tag", storage_name="tagname")
    posamp = RangeTag("Position amplitude")
    minamp = RangeTag("Minimum positional amplitude")
    zamp = RangeTag("Amplitude in Z")
    xamp = RangeTag("Amplitude in X")
    yamp = RangeTag("Amplitude in Y")
    num = RangeTag("Number of atoms/fragments")
    adatom = BoolTag("Add atoms after making supercell")
    fix = BoolTag("FIX this atom")
    nomove = BoolTag("Do not move this atom (even in push)")
    rad = GenericTag("Radius of ion")
    occ = GenericTag("Occupation, can be fractional (e.g 1/3)")
    perm = BoolTag("")
    athole = BoolTag("Place at hole position")
    coord = RangeTag("Coordination of the ion")
    angamp = RangeTag("Angular randomisation magnitude (for fragments)")
    vol = GenericTag("Atomic volume")
    mult = GenericTag("Multiplicity")
    spin = GenericTag("Magnetic spin moment")

    def to_string(self) -> str:
        """
        Return the per entry string for this atom to be appended after its
        line in POSITIONS_FRAC / POSITIONS_ABS block
        """
        if self.disabled is True:
            return ""

        tokens = []
        # Set the tag
        tagname = self.prop_data.get("tagname")
        if not tagname:
            raise ValueError("The tagname property must be set")

        tokens.append("# {} %".format(self.prop_data["tagname"]))
        for key, value in self.prop_data.items():
            if key == "tagname":
                continue

            type_string = self.type_registry[key]
            name = key.upper()

            # Process the value based on type string
            if type_string == "tag":
                tokens.append(name)
            elif type_string == "range":
                if isinstance(value, (list, tuple)):
                    tokens.append("{}={}-{}".format(name, *value))
                else:
                    tokens.append(f"{name}={value}")
            else:
                tokens.append(f"{name}={value}")

        string = " ".join(tokens)
        return string


class SeedAtom(Atom, SeedAtomTag):
    """
    Element atoms in a AIRSS seed
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        SeedAtomTag.__init__(self, *args, **kwargs)
        if self.atoms is not None:
            self.prop_data = self.atoms.arrays["atom_gentags"][self.index].prop_data
            self.type_registry = self.atoms.arrays["atom_gentags"][
                self.index
            ].type_registry


def tuple2range(
    value: Union[numbers.Number, List[numbers.Number], Tuple[numbers.Number, ...]],
) -> str:
    """
    Return the string for a given value. If the value is a tuple
    make it a range.
    """
    if isinstance(value, (list, tuple)):
        return f"{value[0]}-{value[1]}"
    return str(value)


def get_cell_inp(atoms: "SeedAtoms") -> CellInput:
    """Get the CellInput holder for a given seed"""
    cell = CellInput()

    # Prepare the cell out
    cell.set_cell(atoms.get_cell())

    # Prepare the positions
    if atoms.positions.size > 0:
        species = atoms.get_chemical_symbols()
        pos_line_tags = list(atoms.atom_tags)
        tags_lines = [tag.to_string() for tag in pos_line_tags]
        cell.set_positions(species, atoms.get_positions(), tags_lines)

    return cell


def get_cell_inp_lines(atoms: "SeedAtoms") -> List[str]:
    """
    Write the seed to a file handle
    """
    cell = get_cell_inp(atoms)
    # Insert tags in the cell block
    tags = atoms.gentags.get_prop_dict()
    for tag in tags:
        if tag in ["FIX", "CFIX", "ABFIX"]:
            cell["lattice_cart"].append("#" + tag)

    lines = []
    lines.extend(cell.get_file_lines())
    lines.extend(atoms.gentags.to_string().split("\n"))
    return lines

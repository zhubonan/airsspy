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
Tools for handling res files
"""
import os
import re
from collections import namedtuple
from typing import Any, Dict, List, Optional, Union

import numpy as np
from ase import Atoms
from ase.geometry import cellpar_to_cell
from ase.io import write

# pymatgen dependencies for advanced features
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

try:
    from spglib import get_spacegroup

    SPGLIB_AVAILABLE = True
except ImportError:
    SPGLIB_AVAILABLE = False
    get_spacegroup = None


def extract_res(fname: str) -> Dict[str, Union[str, float, int, List[str]]]:
    """
    Extract information from res file.
    Returns a dictionary.
    The structure of he res file is not extracted
    """
    rems: List[str] = []
    title: str = ""
    with open(fname) as fh:
        for line in fh:
            if "TITL" in line:
                title = line.strip()
            if "REM" in line:
                rems.append(line.replace("REM", "").strip())
            # Break when data starts
            if "cell" in line:
                break
    entries = title.split()
    res: Dict[str, Union[str, float, int, List[str]]] = {}
    res["rem"] = rems
    res["uid"] = entries[1]
    res["P"] = float(entries[2])
    res["V"] = float(entries[3])
    res["H"] = float(entries[4])
    res["nat"] = int(entries[7])
    res["sym"] = entries[8]
    res["fname"] = fname
    return res


def save_airss_res(
    atoms: Atoms,
    info_dict: Dict[str, Any],
    fname: Optional[str] = None,
    force_write: bool = False,
) -> None:
    """
    Save the relaxed structure in res format which is compatible with the
    ``cryan`` program.
    """

    # Prepare output file
    if fname is None:
        fname = info_dict["uid"] + ".res"
    if os.path.isfile(fname) and not force_write:
        raise FileExistsError("Switch on force_write to overwrite existing files")

    P, V, H = info_dict["P"], info_dict["V"], info_dict["H"]
    # Get number of atoms, spin
    nat, sg = info_dict["nat"], info_dict["sym"]
    # Construct title line
    pvh = f" {P:.3f} {V:.3f} {H:.6f} "
    title = (
        "TITL "
        + info_dict["uid"]
        + " "
        + pvh
        + " 0 "
        + " 0 "
        + " "
        + str(nat)
        + " "
        + sg
        + " n - 1\n"
    )
    rems = info_dict.get("rem", [])

    # Write to the top of res file
    restmp = info_dict["uid"] + ".rtmp"
    write(restmp, atoms, format="res")
    with open(restmp) as resin:
        with open(fname, "w") as fout:
            # Write the title
            fout.write(title + "\n")
            for line in rems:
                fout.write("REM " + line + "\n")
            # Write the data lines
            for n, line in enumerate(resin):
                if n > 0:
                    fout.write(line)

    os.remove(restmp)
    return


# RES file parsing structures
TITLE_KEYS = [
    "label",
    "pressure",
    "volume",
    "enthalpy",
    "spin",
    "spin_abs",
    "natoms",
    "symm",
    "flag1",
    "flag2",
    "flag3",
]
TitlInfo = namedtuple("TitlInfo", TITLE_KEYS)

# Regular expressions for RES file parsing
RES_COORD_PATT = re.compile(
    r"""(\w+)\s+
                            ([0-9]+)\s+
                            ([0-9\-\.]+)\s+
                            ([0-9\-\.]+)\s+
                            ([0-9\-\.]+)\s+
                            ([0-9\-\.]+)""",
    re.VERBOSE,
)

RES_COORD_PATT_WITH_SPIN = re.compile(
    r"""(\w+)\s+
                            ([0-9]+)\s+
                            ([0-9\-\.]+)\s+
                            ([0-9\-\.]+)\s+
                            ([0-9\-\.]+)\s+
                            ([0-9\-\.]+)\s+
                            ([0-9\-\.]+)""",
    re.VERBOSE,
)


def parse_titl(line: str) -> TitlInfo:
    """Parse TITL line and return a TitlInfo object"""
    tokens = line.split()[1:]
    return TitlInfo(
        label=tokens[0],
        pressure=float(tokens[1]),
        volume=float(tokens[2]),
        enthalpy=float(tokens[3]),
        spin=float(tokens[4]),
        spin_abs=float(tokens[5]),
        natoms=int(tokens[6]),
        symm=tokens[7],
        flag1=tokens[8],
        flag2=tokens[9],
        flag3=tokens[10],
    )


def _read_res(lines: List[str]) -> Dict[str, Any]:
    """
    Read a res file from a list of lines

    Args:
        lines: A list of lines containing RES data

    Returns:
        Dictionary containing parsed RES file data
    """
    abc = []
    ang = []
    species = []
    coords = []

    line_no = 0
    title_items = None
    rem_lines = []
    spins = []

    while line_no < len(lines):
        line = lines[line_no]
        tokens = line.split()
        if not tokens:
            line_no += 1
            continue

        if tokens[0] == "TITL":
            title_items = parse_titl(line)

        elif tokens[0] == "CELL" and len(tokens) == 8:
            abc = [float(tok) for tok in tokens[2:5]]
            ang = [float(tok) for tok in tokens[5:8]]

        elif tokens[0] == "SFAC":
            for atom_line in lines[line_no:]:
                if line.strip() == "END":
                    break

                match = RES_COORD_PATT_WITH_SPIN.search(atom_line)
                if match:
                    has_spin = True
                else:
                    has_spin = False
                    match = RES_COORD_PATT.search(atom_line)

                if match:
                    species.append(match.group(1))
                    xyz = match.groups()[2:5]
                    coords.append([float(c) for c in xyz])
                    if has_spin:
                        spins.append(float(match.group(7)))
                line_no += 1

        elif tokens[0] == "REM":
            rem_lines.append(line[4:].strip())

        line_no += 1

    return {
        "titl": title_items,
        "species": species,
        "scaled_positions": coords,
        "cellpar": list(abc) + list(ang),
        "rem_lines": rem_lines,
        "spins": spins,
    }


def _get_res_lines(
    titl,
    species: List[str],
    scaled_positions: List[List[float]],
    cellpar: List[float],
    rem_lines: Optional[List[str]] = None,
    spins: Optional[List[float]] = None,
) -> List[str]:
    """
    Write RES format lines using given data

    Args:
        titl: A TitlInfo object or list of title items
        species: A list of species for each site
        scaled_positions: Scaled (fractional) atomic positions
        cellpar: Cell parameters in a, b, c, alpha, beta, gamma
        rem_lines: Lines for the REM information
        spins: A list of spins to be added to each site if given

    Returns:
        A list of lines for the RES file
    """
    lines = []

    if isinstance(titl, TitlInfo):
        titl_list = [
            titl.label,
            titl.pressure,
            titl.volume,
            titl.enthalpy,
            titl.spin,
            titl.spin_abs,
            titl.natoms,
            titl.symm,
            titl.flag1,
            titl.flag2,
            titl.flag3,
        ]
    else:
        titl_list = titl

    if len(titl_list) != 11:
        raise ValueError("TITL must be a list of length 11.")

    titl_str = "{} {:.3f} {:.3f} {:.4f} {:.2f} {:.2f} {} ({}) {} {} {}".format(
        *titl_list
    )
    lines.append("TITL " + titl_str)

    if rem_lines:
        for line in rem_lines:
            lines.append("REM " + line)

    lines.append(
        "CELL 1.0 {:<12.6f} {:<12.6f} {:<12.6f} {:<12.6f} {:<12.6f} {:<12.6f}".format(
            *cellpar
        )
    )
    lines.append("LATT -1")

    unique_species = unique(species)
    lines.append("SFAC " + " ".join(unique_species))
    lookup = {s: i + 1 for i, s in enumerate(unique_species)}

    for i, (symbol, pos) in enumerate(zip(species, scaled_positions)):
        line = "{:<4} {:<2} {:>20.13f} {:>20.13f} {:>20.13f} 1.0".format(
            symbol, lookup[symbol], *pos
        )
        if spins:
            line = line + f" {spins[i]:>8.3f}"
        lines.append(line)

    lines.append("END")
    return lines


def unique(items: List[Any]) -> List[Any]:
    """Get a list of ordered unique items"""
    out = []
    for item in items:
        if item not in out:
            out.append(item)
    return out


def read_res_atoms(lines: List[str]) -> tuple[TitlInfo, Atoms]:
    """Read a RES file, return as (TitlInfo, ase.Atoms)"""
    out = _read_res(lines)
    return out["titl"], Atoms(
        symbols=out["species"],
        scaled_positions=out["scaled_positions"],
        cell=cellpar_to_cell(out["cellpar"]),
        pbc=True,
    )


def read_res_pmg(
    lines: List[str],
) -> tuple[TitlInfo, List[str], Optional[Structure], List[float]]:
    """Read a RES file, return as (TitlInfo, rem_lines, pymatgen.Structure, spins)"""
    out = _read_res(lines)
    cell = cellpar_to_cell(out["cellpar"])
    structure = Structure(
        cell, out["species"], out["scaled_positions"], coords_are_cartesian=False
    )
    return out["titl"], out["rem_lines"], structure, out["spins"]


def get_spacegroup_atoms(
    atoms: Atoms, symprec: float = 1e-5, angle_tolerance: float = -1.0
) -> str:
    """Get spacegroup of atoms using spglib"""
    if not SPGLIB_AVAILABLE:
        raise ImportError("spglib is required for spacegroup detection")

    return get_spacegroup(
        (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers()),
        symprec=symprec,
        angle_tolerance=angle_tolerance,
    )


def get_minsep(species: List[str], distance_matrix: np.ndarray) -> Dict[str, float]:
    """
    Calculate minimum separations given species and distance matrix

    Args:
        species: A list of species symbols
        distance_matrix: The distance matrix

    Returns:
        Dictionary of {set(s1, s2): minsep}
    """
    species = np.asarray(species)
    nspec = distance_matrix.shape[0]
    all_minseps = {}

    for i in range(nspec):
        for j in range(i + 1, nspec):
            dist = distance_matrix[i, j]
            spair = sorted([str(species[i]), str(species[j])])
            pair = f"{spair[0]}-{spair[1]}"
            if pair in all_minseps:
                if all_minseps[pair] > dist:
                    all_minseps[pair] = dist
            else:
                all_minseps[pair] = dist
    return all_minseps


def format_minsep(minsep: Dict[str, float]) -> str:
    """Return string representation of the minimum separations"""
    string = ""
    for key, value in minsep.items():
        if isinstance(value, (list, tuple)):
            string += f"{key}={value[0]:.2f}-{value[1]:.2f} "
        else:
            string += f"{key}={value:.2f} "
    return string.strip()


class RESFile:
    """
    Class representing a RES file.

    The SHELX file contains both the structure and computed properties
    as well as metadata.
    """

    def __init__(
        self,
        structure: Union[Atoms, Structure, None],
        data: Dict[str, Any],
        lines: Optional[List[str]] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ):
        """
        Initialize a RESFile object.

        Args:
            structure: ASE Atoms or pymatgen Structure instance
            data: Dictionary containing underlying data
            lines: Raw lines of the RES file
            metadata: Additional metadata
        """
        if isinstance(structure, Atoms):
            structure = AseAtomsAdaptor.get_structure(structure)

        self.structure = structure
        self.lines = lines

        # Set default values if not provided
        if "volume" not in data:
            data["volume"] = structure.volume if structure else None
        if "natoms" not in data:
            data["natoms"] = len(structure) if structure else 0
        if "symm" not in data:
            data["symm"] = structure.get_space_group_info()[0] if structure else None

        self._data = data
        self.metadata = metadata if metadata else {}

    @property
    def rem(self) -> Optional[List[str]]:
        """REM lines"""
        return self._data.get("rem")

    @property
    def atoms(self) -> Optional[Atoms]:
        """Returns an ASE Atoms object"""
        if not self.structure:
            return None
        return AseAtomsAdaptor.get_atoms(self.structure)

    @property
    def data(self) -> Dict[str, Any]:
        """Underlying data of the object"""
        return self._data

    @property
    def label(self) -> Optional[str]:
        """Label of the structure"""
        return self._data.get("label")

    @property
    def name(self) -> Optional[str]:
        """Alias for label"""
        return self.label

    @property
    def enthalpy(self) -> Optional[float]:
        """Enthalpy as reported"""
        return self._data.get("enthalpy")

    @property
    def volume(self) -> Optional[float]:
        """Volume as reported"""
        return self._data.get("volume")

    @property
    def pressure(self) -> float:
        """External pressure as reported"""
        return self._data.get("pressure", 0.0)

    @property
    def natoms(self) -> Optional[int]:
        """Number of atoms"""
        return self._data.get("natoms")

    @property
    def symm(self) -> Optional[str]:
        """Symmetry as reported"""
        return self._data.get("symm")

    @property
    def spin(self) -> float:
        """Spin as reported"""
        return self._data.get("spin", 0.0)

    @property
    def spins(self) -> List[float]:
        """Spin values for each atom"""
        return self._data.get("spins", [])

    @property
    def spin_abs(self) -> float:
        """Absolute integrated spin"""
        return self._data.get("spin_abs", 0.0)

    @property
    def composition(self) -> Optional[Any]:
        """Composition of the structure"""
        if not self.structure:
            return None
        return self.structure.composition

    @property
    def formula(self) -> str:
        """Formula of the structure"""
        if not self.structure:
            return "Unknown"
        return self.composition.formula.replace(" ", "")

    @property
    def reduced_formula(self) -> str:
        """Reduced formula of the structure"""
        if not self.structure:
            return "Unknown"
        return self.composition.reduced_formula

    @classmethod
    def from_string(cls, string: str) -> "RESFile":
        """Construct from a string"""
        return cls.from_lines(string.split("\n"))

    @classmethod
    def from_lines(
        cls, lines: List[str], include_structure: bool = True, only_titl: bool = False
    ) -> "RESFile":
        """Construct from lines"""
        if include_structure:
            titls, rem_lines, structure, spins = read_res_pmg(lines)
            data = {
                "rem": rem_lines,
                "spins": spins,
                **titls._asdict(),
            }
        elif only_titl:
            for line in lines:
                if "TITL" in line:
                    titls = parse_titl(line)
                    break
            else:
                titls = None
            structure = None
            data = titls._asdict() if titls else {}
        else:
            output = _read_res(lines)
            data = {
                "rem_lines": output.get("rem_lines", []),
                "spins": output.get("spins", []),
            }
            if output["titl"]:
                data.update(output["titl"]._asdict())
            structure = None

        obj = cls(structure, data, lines=lines)
        return obj

    def load_structure(self) -> None:
        """Load structure from the lines"""
        if not self.lines:
            raise ValueError("No lines available to load structure from")
        new_obj = self.from_lines(self.lines, include_structure=True)
        self.structure = new_obj.structure
        self._data = new_obj.data

    @classmethod
    def from_file(
        cls, fname: str, include_structure: bool = True, only_titl: bool = False
    ) -> "RESFile":
        """Construct from a file"""
        with open(fname) as fhandle:
            return cls.from_lines(
                fhandle.readlines(),
                include_structure=include_structure,
                only_titl=only_titl,
            )

    def to_res_lines(self) -> List[str]:
        """Get the raw RES representation of this object"""

        species = [site.symbol for site in self.structure.species]
        frac_pos = [row.tolist() for row in self.structure.frac_coords]
        cellpar = self.structure.lattice.parameters

        titl = [
            self.label,
            self.pressure,
            self.volume,
            self.enthalpy,
            self.spin,
            self.spin_abs,
            self.natoms,
            self.symm,
            "n",
            "-",
            "1",
        ]

        lines = _get_res_lines(titl, species, frac_pos, cellpar, self.rem, self.spins)
        lines.append("")  # Add trailing newline
        return lines

    def get_minsep(self, string: bool = False) -> Union[Dict[str, float], str]:
        """Return species-wise minimum separations"""

        minsep = get_minsep(self.structure.species, self.structure.distance_matrix)
        if string:
            return format_minsep(minsep)
        return minsep

    def __repr__(self) -> str:
        return f"<RESFile with label={self.label}, formula={self.formula}, enthalpy={self.enthalpy}...>"

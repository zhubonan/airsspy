"""ABACUS output parsing and result composition.

Utilities for parsing ABACUS log files, STRU structure files, and
composing result documents compatible with the AIRSS jobflow pipeline.

Supports both ABACUS develop (v3.9.x) and LTS (v3.10.x) versions which
differ in log output format for energy and SCF convergence messages.
"""

import logging
import re
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

# Conversion constants
BOHR_TO_ANG = 0.529177249
BOHR3_TO_ANG3 = 0.148184743
GPA_TO_EV_PER_ANG3 = 0.006241510219780177
GPA_KBAR_TO_EV_PER_ANG3 = GPA_TO_EV_PER_ANG3


def parse_abacus_log(logfile: str) -> dict:
    """Parse an ABACUS running log file for key results.

    Handles format differences between ABACUS develop (v3.9.x) and
    LTS (v3.10.x).

    Args:
        logfile: Path to the ABACUS log file (e.g. OUT.ABACUS/running_cell-relax.log).

    Returns:
        Dictionary with keys: energy, pressure, volume, converged,
        scf_converged, n_ionic_steps.
    """
    result = {
        "energy": None,
        "pressure": None,
        "volume": None,
        "converged": False,
        "scf_converged": False,
        "n_ionic_steps": 0,
    }

    if not Path(logfile).is_file():
        logger.warning("ABACUS log file not found: %s", logfile)
        return result

    with open(logfile) as fh:
        content = fh.read()

    # Energy — try in order of specificity
    # 1. !FINAL_ETOT_IS <val> eV (both versions)
    m = re.search(r"!FINAL_ETOT_IS\s+([-eE+0-9.]+)\s*eV", content)
    if m:
        result["energy"] = float(m.group(1))
    else:
        # 2. #TOTAL ENERGY# <val> eV (develop)
        m = re.search(r"#TOTAL ENERGY#\s+([-eE+0-9.]+)\s*eV", content)
        if m:
            result["energy"] = float(m.group(1))
        else:
            # 3. final etot is <val> eV (LTS)
            m = re.search(r"final etot is\s+([-eE+0-9.]+)\s*eV", content)
            if m:
                result["energy"] = float(m.group(1))

    # Pressure — use the last occurrence (final ionic step value)
    # LTS: #TOTAL-PRESSURE# (EXCLUDE KINETIC PART OF IONS): <val> GPa
    matches = re.findall(
        r"#TOTAL-PRESSURE#.*?([-eE+0-9.]+)\s*GPa", content, re.IGNORECASE
    )
    if matches:
        result["pressure"] = float(matches[-1])
    else:
        # develop: TOTAL-PRESSURE: <val> KBAR
        matches = re.findall(r"TOTAL-PRESSURE:\s*([-eE+0-9.]+)\s*KBAR", content)
        if matches:
            result["pressure"] = float(matches[-1]) / 10.0

    # Volume
    m = re.search(r"Cell volume \(A\^3\)\s*=\s*([-eE+0-9.]+)", content)
    if m:
        result["volume"] = float(m.group(1))
    else:
        m = re.search(r"Cell volume \(Bohr\^3\)\s*=\s*([-eE+0-9.]+)", content)
        if m:
            result["volume"] = float(m.group(1)) * BOHR3_TO_ANG3

    # Ionic relaxation convergence
    if "Relaxation is converged!" in content:
        result["converged"] = True

    # SCF convergence
    if "#SCF IS CONVERGED#" in content:
        result["scf_converged"] = True
    elif "charge density convergence is achieved" in content:
        result["scf_converged"] = True

    # Ionic step count
    result["n_ionic_steps"] = content.count("STEP OF RELAXATION")

    return result


_AXES = {"x": 0, "y": 1, "z": 2}
_VECTOR_INPUT_KEYS = {"kspacing", "press", "force", "stress"}


def _axis_permutation(axis_map: Optional[str]) -> list[int]:
    """Return old-axis indices ordered by new x,y,z axes."""
    if axis_map is None or not str(axis_map).strip():
        return [0, 1, 2]

    text = str(axis_map).strip().lower()
    if text in ("identity", "none"):
        return [0, 1, 2]

    mapping: dict[int, int] = {}
    for item in text.replace(",", " ").split():
        if ":" not in item:
            raise ValueError(
                "cell_axis_map entries must use old:new syntax, e.g. 'z:x'"
            )
        old, new = [part.strip() for part in item.split(":", 1)]
        if old not in _AXES or new not in _AXES:
            raise ValueError(f"Invalid axis in cell_axis_map entry: {item!r}")
        new_axis = _AXES[new]
        if new_axis in mapping:
            raise ValueError(f"Duplicate target axis in cell_axis_map: {new!r}")
        mapping[new_axis] = _AXES[old]

    remaining_old = [axis for axis in range(3) if axis not in mapping.values()]
    for new_axis in range(3):
        if new_axis not in mapping:
            mapping[new_axis] = remaining_old.pop(0)
    perm = [mapping[i] for i in range(3)]
    if sorted(perm) != [0, 1, 2]:
        raise ValueError(f"cell_axis_map is not one-to-one: {axis_map!r}")
    return perm


def apply_cell_axis_map_to_cell_text(
    cell_content: str,
    cell_axis_map: Optional[str] = None,
) -> str:
    """Permute lattice and coordinates in CASTEP cell text.

    The map uses ``old:new`` axis names.  For example ``z:x`` makes the old
    AIRSS slab-normal/vacuum axis become the new ABACUS x axis, while the
    unspecified old x/y axes fill new y/z in order.
    """
    perm = _axis_permutation(cell_axis_map)
    if perm == [0, 1, 2]:
        return cell_content
    if any(
        line.strip().upper().startswith("%BLOCK LATTICE_ABC")
        for line in cell_content.splitlines()
    ):
        raise ValueError("cell_axis_map requires LATTICE_CART, not LATTICE_ABC")

    lines = cell_content.splitlines()
    out: list[str] = []
    in_block: Optional[str] = None
    lattice_rows: list[list[float]] = []

    def _flush_lattice() -> None:
        if not lattice_rows:
            return
        lattice = np.asarray(lattice_rows, dtype=float)
        mapped = lattice[perm, :][:, perm]
        for vec in mapped:
            out.append("  " + "  ".join(f"{v:.10f}" for v in vec))
        lattice_rows.clear()

    for line in lines:
        stripped = line.strip()
        upper = stripped.upper()
        if upper.startswith("%BLOCK LATTICE_CART"):
            in_block = "lattice_cart"
            out.append(line)
            continue
        if upper.startswith("%BLOCK POSITIONS_FRAC"):
            in_block = "positions_frac"
            out.append(line)
            continue
        if upper.startswith("%BLOCK POSITIONS_ABS"):
            in_block = "positions_abs"
            out.append(line)
            continue
        if upper.startswith("%ENDBLOCK"):
            if in_block == "lattice_cart":
                _flush_lattice()
            in_block = None
            out.append(line)
            continue

        if in_block == "lattice_cart" and stripped and not stripped.startswith("#"):
            lattice_rows.append([float(x) for x in line.split()[:3]])
            continue
        if in_block in ("positions_frac", "positions_abs") and stripped:
            parts = line.split()
            if len(parts) >= 4 and parts[0][0].isalpha():
                coords = np.asarray([float(x) for x in parts[1:4]], dtype=float)[perm]
                rest = " " + " ".join(parts[4:]) if len(parts) > 4 else ""
                out.append(
                    f"{parts[0]} {coords[0]:.10f} {coords[1]:.10f} {coords[2]:.10f}{rest}"
                )
                continue
        out.append(line)
    return "\n".join(out) + ("\n" if cell_content.endswith("\n") else "")


def apply_cell_axis_map_to_abacus_input(
    input_content: str,
    cell_axis_map: Optional[str] = None,
) -> str:
    """Permute direction-vector ABACUS INPUT values to match a cell axis map."""
    perm = _axis_permutation(cell_axis_map)
    if perm == [0, 1, 2]:
        return input_content

    out: list[str] = []
    for line in input_content.splitlines():
        data, sep, comment = line.partition("#")
        parts = data.split()
        if len(parts) >= 4 and parts[0].lower() in _VECTOR_INPUT_KEYS:
            try:
                for value in parts[1:4]:
                    float(value)
            except ValueError:
                out.append(line)
                continue
            mapped = [parts[1:4][old] for old in perm]
            trailing = parts[4:]
            new_data = " ".join([parts[0], *mapped, *trailing])
            if sep:
                new_data += " #" + comment
            out.append(new_data)
        else:
            out.append(line)
    return "\n".join(out) + ("\n" if input_content.endswith("\n") else "")


def cell_to_stru(cell_content: str, cell_axis_map: Optional[str] = None) -> str:
    """Convert CASTEP .cell content to ABACUS STRU format.

    Parses the LATTICE_CART, POSITIONS_FRAC, and SPECIES_POT blocks
    from a .cell file and produces an ABACUS STRU file string.

    Args:
        cell_content: Content of the .cell file.
        cell_axis_map: Optional ``old:new`` axis map applied before conversion.

    Returns:
        STRU file content as a string.
    """
    cell_content = apply_cell_axis_map_to_cell_text(cell_content, cell_axis_map)
    lines = cell_content.splitlines()

    # Parse lattice (LATTICE_CART or LATTICE_ABC)
    lattice = []
    in_block = None  # 'cart' or 'abc'
    abc_vals = []
    for line in lines:
        stripped = line.strip()
        upper = stripped.upper()
        if upper.startswith("%BLOCK LATTICE_CART"):
            in_block = "cart"
            continue
        elif upper.startswith("%BLOCK LATTICE_ABC"):
            in_block = "abc"
            continue
        if upper.startswith("%ENDBLOCK") and in_block:
            in_block = None
            continue
        if in_block == "cart" and stripped:
            lattice.append([float(x) for x in stripped.split()[:3]])
        elif in_block == "abc" and stripped:
            abc_vals.append([float(x) for x in stripped.split()[:3]])

    if abc_vals and not lattice:
        # Convert ABC (lengths + angles) to Cartesian vectors
        from ase.cell import Cell

        lattice = Cell.new(abc_vals[0] + abc_vals[1]).array.tolist()

    if len(lattice) != 3:
        raise ValueError(f"Expected 3 lattice vectors, got {len(lattice)}")

    def _parse_positions_block(block_name):
        parsed_elements = []
        parsed_positions = []
        in_positions = False
        for line in lines:
            stripped = line.strip()
            upper = stripped.upper()
            if upper.startswith(f"%BLOCK {block_name}"):
                in_positions = True
                continue
            if upper.startswith("%ENDBLOCK"):
                if in_positions:
                    break
                continue
            if in_positions and stripped:
                parts = stripped.split()
                parsed_elements.append(parts[0])
                parsed_positions.append([float(x) for x in parts[1:4]])
        return parsed_elements, parsed_positions

    # Parse positions. CRUD rebuilds claimed structures as POSITIONS_ABS, while
    # hand-written ABACUS seeds commonly use POSITIONS_FRAC.
    elements, positions = _parse_positions_block("POSITIONS_FRAC")
    position_mode = "Direct"
    if not elements:
        elements, positions = _parse_positions_block("POSITIONS_ABS")
        position_mode = "Cartesian"

    # Parse SPECIES_POT — pairs of (orbital, pseudopotential) per element
    # Format (LCAO): Element orbital_file / Element pseudopotential_file
    # Format (PW):   Element pseudopotential_file
    species_pot = {}
    in_block = False
    pot_lines = []
    for line in lines:
        stripped = line.strip()
        if stripped.upper().startswith("%BLOCK SPECIES_POT"):
            in_block = True
            continue
        if stripped.upper().startswith("%ENDBLOCK"):
            if in_block:
                break
            continue
        if in_block and stripped:
            # Strip the leading element symbol
            parts = stripped.split(None, 1)
            pot_lines.append(parts[1] if len(parts) > 1 else stripped)

    # Group by element:
    # If entries come in pairs (orbital + upf per element), use (orb, upf)
    # If one entry per element, use (None, upf) — PW mode, no orbital file
    unique_elements = list(dict.fromkeys(elements))  # preserve order
    n_elem = len(unique_elements)
    if len(pot_lines) >= n_elem * 2:
        # Pairs: (orbital, upf) for each element
        for idx, elem in enumerate(unique_elements):
            orb_idx = idx * 2
            upf_idx = idx * 2 + 1
            if orb_idx < len(pot_lines) and upf_idx < len(pot_lines):
                species_pot[elem] = (pot_lines[orb_idx], pot_lines[upf_idx])
    elif len(pot_lines) >= n_elem:
        # Single: just upf for each element (PW mode)
        for idx, elem in enumerate(unique_elements):
            if idx < len(pot_lines):
                species_pot[elem] = (None, pot_lines[idx])

    # Build STRU content
    out_lines = []
    out_lines.append("ATOMIC_SPECIES")
    for elem in unique_elements:
        if elem in species_pot:
            orb, upf = species_pot[elem]
            mass = _ATOMIC_MASSES.get(elem, 1.0)
            out_lines.append(f"{elem} {mass} {upf}")
    out_lines.append("")

    out_lines.append("LATTICE_CONSTANT")
    out_lines.append(f"{1.0 / BOHR_TO_ANG:.10f}")
    out_lines.append("")

    out_lines.append("LATTICE_VECTORS")
    for vec in lattice:
        out_lines.append("  ".join(f"{v:.10f}" for v in vec))
    out_lines.append("")

    out_lines.append("ATOMIC_POSITIONS")
    out_lines.append(position_mode)
    for elem in unique_elements:
        out_lines.append(elem)
        out_lines.append("0.0")  # magnetization
        n = elements.count(elem)
        out_lines.append(str(n))
        for e, pos in zip(elements, positions):
            if e == elem:
                out_lines.append(f"{pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f} 1 1 1")

    # Only add NUMERICAL_ORBITAL for LCAO basis (where orbital files exist)
    orbital_entries = [
        species_pot[elem][0]
        for elem in unique_elements
        if elem in species_pot and species_pot[elem][0] is not None
    ]
    if orbital_entries:
        out_lines.append("")
        out_lines.append("NUMERICAL_ORBITAL")
        for orb in orbital_entries:
            out_lines.append(orb)

    return "\n".join(out_lines) + "\n"


# Minimal atomic masses for ABACUS ATOMIC_SPECIES block
_ATOMIC_MASSES = {
    "H": 1.008,
    "He": 4.003,
    "Li": 6.941,
    "Be": 9.012,
    "B": 10.811,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "Ne": 20.180,
    "Na": 22.990,
    "Mg": 24.305,
    "Al": 26.982,
    "Si": 28.086,
    "P": 30.974,
    "S": 32.065,
    "Cl": 35.453,
    "Ar": 39.948,
    "K": 39.098,
    "Ca": 40.078,
    "Sc": 44.956,
    "Ti": 47.867,
    "V": 50.942,
    "Cr": 51.996,
    "Mn": 54.938,
    "Fe": 55.845,
    "Co": 58.933,
    "Ni": 58.693,
    "Cu": 63.546,
    "Zn": 65.380,
    "Ga": 69.723,
    "Ge": 72.630,
    "As": 74.922,
    "Se": 78.971,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.468,
    "Sr": 87.620,
    "Y": 88.906,
    "Zr": 91.224,
    "Nb": 92.906,
    "Mo": 95.950,
    "Tc": 98.000,
    "Ru": 101.070,
    "Rh": 102.906,
    "Pd": 106.420,
    "Ag": 107.868,
    "Cd": 112.414,
    "In": 114.818,
    "Sn": 118.711,
    "Sb": 121.760,
    "Te": 127.600,
    "I": 126.904,
    "Xe": 131.293,
    "Cs": 132.905,
    "Ba": 137.327,
    "La": 138.905,
    "Ce": 140.116,
    "Pr": 140.908,
    "Nd": 144.242,
    "Pm": 145.000,
    "Sm": 150.360,
    "Eu": 151.964,
    "Gd": 157.250,
    "Tb": 158.925,
    "Dy": 162.500,
    "Ho": 164.930,
    "Er": 167.259,
    "Tm": 168.934,
    "Yb": 173.045,
    "Lu": 174.967,
    "Hf": 178.490,
    "Ta": 180.948,
    "W": 183.840,
    "Re": 186.207,
    "Os": 190.230,
    "Ir": 192.217,
    "Pt": 195.084,
    "Au": 196.967,
    "Hg": 200.590,
    "Tl": 204.383,
    "Pb": 207.200,
    "Bi": 208.980,
    "Po": 209.000,
    "At": 210.000,
    "Rn": 222.000,
    "Fr": 223.000,
    "Ra": 226.000,
    "Ac": 227.000,
    "Th": 232.038,
    "Pa": 231.036,
    "U": 238.029,
    "Np": 237.000,
    "Pu": 244.000,
}


def parse_abacus_stru(stru_path: str):
    """Parse an ABACUS STRU file to extract structure information.

    Args:
        stru_path: Path to the STRU file.

    Returns:
        Tuple of (elements, positions, cell) where:
        - elements: list of element symbols (str)
        - positions: numpy array of fractional positions (N x 3)
        - cell: numpy array of cell vectors in Angstrom (3 x 3)
    """
    with open(stru_path) as fh:
        lines = fh.readlines()

    lattice_constant = 1.0
    lattice_vectors = None
    coord_type = "Direct"
    elements = []
    positions = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("LATTICE_CONSTANT"):
            # Next non-empty line is the value
            i += 1
            while i < len(lines) and not lines[i].strip():
                i += 1
            if i < len(lines):
                lattice_constant = float(lines[i].strip())
        elif line.startswith("LATTICE_VECTORS"):
            vecs = []
            i += 1
            while i < len(lines) and len(vecs) < 3:
                line = lines[i].strip()
                if line:
                    vecs.append([float(x) for x in line.split()[:3]])
                i += 1
            lattice_vectors = np.array(vecs)
            # Convert from Bohr to Angstrom if LATTICE_CONSTANT ~ 1.8897
            if lattice_constant > 1.5:
                lattice_vectors *= lattice_constant * BOHR_TO_ANG
            else:
                lattice_vectors *= lattice_constant
            continue
        elif line.startswith("ATOMIC_POSITIONS"):
            i += 1
            # Next non-empty line is the coordinate type
            while i < len(lines) and not lines[i].strip():
                i += 1
            if i < len(lines):
                coord_type = lines[i].strip()
            i += 1
            # Parse element blocks
            while i < len(lines):
                line = lines[i].strip()
                if not line or line.startswith("#"):
                    i += 1
                    continue
                # Check if this is a new section header
                if line in (
                    "ATOMIC_SPECIES",
                    "LATTICE_CONSTANT",
                    "LATTICE_VECTORS",
                    "ATOMIC_POSITIONS",
                    "NUMERICAL_ORBITAL",
                ):
                    break
                # Element name (strip ABACUS #label suffix)
                current_element = line.split()[0]
                i += 1
                # Magnetization (skip)
                while i < len(lines) and not lines[i].strip():
                    i += 1
                i += 1
                # Number of atoms
                while i < len(lines) and not lines[i].strip():
                    i += 1
                n_atoms = int(lines[i].strip().split()[0])
                i += 1
                # Coordinate lines
                for _ in range(n_atoms):
                    while i < len(lines) and not lines[i].strip():
                        i += 1
                    if i < len(lines):
                        parts = lines[i].strip().split()
                        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                        elements.append(current_element)
                        positions.append([x, y, z])
                        i += 1

        i += 1

    positions = np.array(positions)

    # Convert Cartesian to fractional if needed
    if coord_type.lower() == "cartesian" and lattice_vectors is not None:
        inv_cell = np.linalg.inv(lattice_vectors)
        positions = positions @ inv_cell

    return elements, positions, lattice_vectors


def detect_logfile(workdir: str, input_path: str) -> Optional[str]:
    """Detect the ABACUS log file path based on calculation type.

    Args:
        workdir: Path to the ABACUS working directory.
        input_path: Path to the INPUT file.

    Returns:
        Path to the log file, or None if not found.
    """
    calc_type = None
    if Path(input_path).is_file():
        with open(input_path) as fh:
            for line in fh:
                m = re.match(r"^\s*calculation\s+(\S+)", line)
                if m:
                    calc_type = m.group(1)
                    break

    if calc_type is None:
        calc_type = "cell-relax"

    logbase = f"running_{calc_type}.log"
    logfile = Path(workdir) / "OUT.ABACUS" / logbase
    if logfile.is_file():
        return str(logfile)

    # Fallback: try any running_*.log
    out_dir = Path(workdir) / "OUT.ABACUS"
    if out_dir.is_dir():
        logs = sorted(out_dir.glob("running_*.log"))
        if logs:
            return str(logs[-1])

    return None


def _parse_external_pressure_gpa(input_path: str) -> Optional[float]:
    """Parse hydrostatic external pressure from ABACUS INPUT press1/2/3.

    ABACUS takes ``press1``, ``press2`` and ``press3`` in kbar.  The AIRSS
    result line stores scalar pressure in GPa.
    """
    path = Path(input_path)
    if not path.is_file():
        return None

    values: list[float] = []
    with open(path) as fh:
        for line in fh:
            stripped = line.split("#", 1)[0].strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            if parts[0].lower() in ("press1", "press2", "press3"):
                try:
                    values.append(float(parts[1]))
                except ValueError:
                    pass

    if not values:
        return None
    return sum(values) / len(values) / 10.0


def extract_abacus_rem(struct_name: str) -> dict:
    """Extract quality-affecting computational parameters from ABACUS output.

    Reads the ABACUS INPUT file and log file to extract metadata for
    the REM block of a .res file.

    Args:
        struct_name: Structure name (without extension).

    Returns:
        Dictionary with keys: functional, cutoff, basis_type, kspacing,
        nkpts, psps, orbital_info.
    """
    rem: dict = {}
    workdir = f"{struct_name}.abacus"

    # Parse INPUT file
    input_path = struct_name + ".INPUT"
    if Path(input_path).is_file():
        with open(input_path) as fh:
            for line in fh:
                m = re.match(r"^\s*ecutwfc\s+(\S+)", line)
                if m:
                    rem["cutoff"] = float(m.group(1))
                m = re.match(r"^\s*basis_type\s+(\S+)", line)
                if m:
                    rem["basis_type"] = m.group(1).strip()
                m = re.match(r"^\s*kspacing\s+(\S+)", line)
                if m:
                    rem["kspacing"] = float(m.group(1))
                m = re.match(r"^\s*dft_functional\s+(\S+)", line)
                if m:
                    rem["functional"] = "Functional " + m.group(1).strip()

    # Parse log file
    logfile = detect_logfile(workdir, input_path)
    if logfile and Path(logfile).is_file():
        with open(logfile) as fh:
            in_atom_type = False
            atom_type_orbital: list[str] = []
            for line in fh:
                # Functional from pseudopotential section
                if (
                    "exchange-correlation functional" in line
                    and "functional" not in rem
                ):
                    m = re.search(r"=\s*(\S+)", line)
                    if m:
                        rem["functional"] = "Functional " + m.group(1)
                # Number of k-points (last occurrence)
                if "nkstot" in line:
                    m = re.search(r"nkstot\s*=\s*(\d+)", line)
                    if m:
                        rem["nkpts"] = int(m.group(1))
                # Pseudopotential file names
                if "Read in pseudopotential file is" in line:
                    m = re.search(r"is\s+(\S+)", line)
                    if m:
                        if "psps" not in rem:
                            rem["psps"] = []
                        rem["psps"].append(m.group(1))
                # Orbital zeta info for LCAO basis
                if re.match(r"\s*READING ATOM TYPE", line):
                    in_atom_type = True
                    atom_type_orbital = []
                if in_atom_type:
                    if re.match(r"\s*atom label", line):
                        label = line.split("=")[-1].strip()
                        atom_type_orbital.append(label)
                    if "number of zeta" in line:
                        m = re.search(r"L=(\d+),\s*number of zeta\s*=\s*(\d+)", line)
                        if m and atom_type_orbital:
                            atom_type_orbital.append(f"L{m.group(1)}-dz{m.group(2)}")
                    if re.match(r"\s+number of atom", line) or (
                        "TOTAL ATOM NUMBER" in line
                    ):
                        if atom_type_orbital and len(atom_type_orbital) >= 2:
                            elem = atom_type_orbital[0]
                            zetas = " ".join(atom_type_orbital[1:])
                            if "orbital_info" not in rem:
                                rem["orbital_info"] = []
                            rem["orbital_info"].append(f"{elem}: {zetas}")
                        in_atom_type = False
                        atom_type_orbital = []

    return rem


def build_abacus_rem_lines(struct_name: str) -> list[str]:
    """Build REM block lines for a .res file from ABACUS output.

    Args:
        struct_name: Structure name (without extension).

    Returns:
        List of REM line strings (without "REM " prefix).
    """
    rem = extract_abacus_rem(struct_name)

    lines: list[str] = []
    lines.append("")

    # Functional
    if "functional" in rem:
        lines.append(rem["functional"])

    # Basis type + cutoff + kspacing
    parts = []
    if "basis_type" in rem:
        parts.append("Basis " + rem["basis_type"])
    if "cutoff" in rem:
        parts.append(f"Cut-off {rem['cutoff']} Ry")
    if "kspacing" in rem:
        parts.append(f"Spacing {rem['kspacing']}")
    if parts:
        lines.append(" ".join(parts))

    # Number of k-points
    if "nkpts" in rem:
        lines.append(f"No. kpts {rem['nkpts']}")

    lines.append("")

    # Pseudopotentials
    if "psps" in rem:
        for psp in rem["psps"]:
            lines.append(psp)

    # Orbital info (LCAO basis)
    if "orbital_info" in rem:
        lines.append("")
        for orb in rem["orbital_info"]:
            lines.append(orb)

    lines.append("")
    return lines


def compose_abacus_task_doc(struct_name: str) -> dict:
    """Extract results from a completed ABACUS calculation.

    Reads the ABACUS log and STRU_ION_D output, creates an ASE Atoms
    object, saves a ``.res`` file, and returns a dictionary suitable
    for constructing an ``AirssResultDoc``.

    Args:
        struct_name: Structure name (without extension).

    Returns:
        Dictionary with energy, structure, volume, formula, etc.
    """
    workdir = f"{struct_name}.abacus"
    input_path = f"{struct_name}.INPUT"

    # Detect and parse log file
    logfile = detect_logfile(workdir, input_path)
    if not logfile:
        raise RuntimeError(
            f"ABACUS log file not found for {struct_name}; refusing to collect result"
        )
    log_data = parse_abacus_log(logfile)
    if not log_data.get("scf_converged", False):
        raise RuntimeError(
            f"ABACUS SCF did not converge for {struct_name}; refusing to collect result"
        )

    from ase import Atoms
    from pymatgen.io.ase import AseAtomsAdaptor

    from .restools import save_airss_res

    energy = log_data.get("energy")
    logged_pressure = log_data.get("pressure")
    external_pressure = _parse_external_pressure_gpa(input_path)
    pressure = external_pressure if external_pressure is not None else logged_pressure
    volume = log_data.get("volume")

    # Read relaxed structure from STRU_ION_D
    stru_path = Path(workdir) / "OUT.ABACUS" / "STRU_ION_D"
    if not stru_path.is_file():
        stru_path = Path(workdir) / "STRU"

    elements = []
    positions = np.zeros((0, 3))
    cell = np.eye(3)

    if stru_path.is_file():
        elements, positions, cell = parse_abacus_stru(str(stru_path))

    atoms = Atoms(symbols=elements, positions=positions @ cell, cell=cell, pbc=True)

    # Compute symmetry via spglib
    try:
        import spglib

        sg = spglib.get_spacegroup(
            (
                atoms.get_cell().array,
                atoms.get_scaled_positions(),
                atoms.get_atomic_numbers(),
            ),
            symprec=0.1,
        )
        sym = sg.split()[0] if sg else "P1"
    except (ImportError, Exception):
        sym = "P1"

    # Compute enthalpy for pressure results
    enthalpy = energy
    if energy is not None and pressure is not None and volume is not None:
        enthalpy = energy + pressure * volume * GPA_TO_EV_PER_ANG3

    # Build REM lines from ABACUS metadata
    rem_lines = build_abacus_rem_lines(struct_name)

    info = {
        "uid": struct_name,
        "P": pressure if pressure is not None else 0.0,
        "V": atoms.get_volume(),
        "H": enthalpy if enthalpy is not None else 0.0,
        "nat": len(atoms),
        "sym": sym,
        "rem": rem_lines,
    }
    save_airss_res(atoms, info, fname=struct_name + ".res", force_write=True)

    structure = AseAtomsAdaptor.get_structure(atoms)

    # Read total time from ABACUS stdout capture
    total_time = None
    stdout_path = Path(workdir) / "abacus_out"
    if stdout_path.is_file():
        with open(stdout_path) as fh:
            for line in fh:
                m = re.search(r"TOTAL\s+Time\s*:\s*([0-9.]+)", line)
                if m:
                    total_time = float(m.group(1))

    return {
        "structure": structure,
        "volume": structure.volume,
        "reduced_formula": structure.reduced_formula,
        "formula": structure.composition.formula.replace(" ", ""),
        "natoms": len(atoms),
        "label": struct_name,
        "energy": energy,
        "energy_per_atom": energy / len(atoms) if energy else None,
        "pressure": pressure,
        "total_time": total_time,
        "res_content": Path(struct_name + ".res").read_text()
        if Path(struct_name + ".res").is_file()
        else None,
        "rem_lines": rem_lines,
    }

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
GPA_KBAR_TO_EV_PER_ANG3 = 0.006241510219780177


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

    # Pressure
    # LTS: #TOTAL-PRESSURE# (EXCLUDE KINETIC PART OF IONS): <val> GPa
    m = re.search(r"#TOTAL-PRESSURE#.*?([-eE+0-9.]+)\s*GPa", content, re.IGNORECASE)
    if m:
        result["pressure"] = float(m.group(1))
    else:
        # develop: TOTAL-PRESSURE: <val> KBAR
        m = re.search(r"TOTAL-PRESSURE:\s*([-eE+0-9.]+)\s*KBAR", content)
        if m:
            result["pressure"] = float(m.group(1)) / 10.0

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
                if not line:
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
                # Element name
                current_element = line
                i += 1
                # Magnetization (skip)
                while i < len(lines) and not lines[i].strip():
                    i += 1
                i += 1
                # Number of atoms
                while i < len(lines) and not lines[i].strip():
                    i += 1
                n_atoms = int(lines[i].strip())
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
        positions = positions @ inv_cell.T

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
    from ase import Atoms
    from pymatgen.io.ase import AseAtomsAdaptor

    from .restools import save_airss_res

    workdir = f"{struct_name}.abacus"
    input_path = f"{struct_name}.INPUT"

    # Detect and parse log file
    logfile = detect_logfile(workdir, input_path)
    log_data = {}
    if logfile:
        log_data = parse_abacus_log(logfile)

    energy = log_data.get("energy")
    pressure = log_data.get("pressure")
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

    # Compute enthalpy for pressure results
    enthalpy = energy
    if energy is not None and pressure is not None and volume is not None:
        enthalpy = energy + pressure * volume * GPA_KBAR_TO_EV_PER_ANG3

    info = {"uid": struct_name, "H": enthalpy}
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
    }

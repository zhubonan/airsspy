"""
Self-consistent CASTEP relaxation driver.

Replicates the behaviour of ``castep_relax.pl``, performing iterative
fixed-NPW geometry optimisations with restart capability and optional
alternating cell constraints.
"""

import fileinput
import json
import logging
import os
import re
import shutil
import subprocess
from typing import Optional

import numpy as np

from .casteptools import push_cell
from .common import RelaxError
from .utils import filter_out_stream

logger = logging.getLogger(__name__)

# CODATA 2002 physical constants for unit conversion
_units = {
    "hbar": 6.58211915e-16,  # eVs
    "Eh": 27.2113845,  # eV
    "kB": 8.617343e-5,  # eV/K
    "a0": 0.5291772108,  # Angstrom
    "c": 299792458,  # m/s
    "e": 1.60217653e-19,  # C
    "me": 5.4857990945e-4,  # u
}
_units["t0"] = _units["hbar"] / _units["Eh"]
_units["Pascal"] = _units["e"] * 1e30


class FullRelax:
    """
    Self-consistent CASTEP relaxation with restart and state tracking.

    Performs iterative fixed-NPW geometry optimisations. Initial cycles
    use a short iteration count, then full-length cycles until two
    consecutive successful runs confirm convergence.
    """

    def __init__(
        self,
        exe: str,
        struct_name: str,
        maxit: int,
        initial_cycle: int = 4,
        initial_length: int = 4,
        alter_cell_cons: bool = False,
    ) -> None:
        """
        Initialise a FullRelax instance.

        Args:
            exe: CASTEP executable name/path.
            struct_name: Name of the structure to relax (without extension).
            maxit: Maximum geometry iterations (0 for single-point).
            initial_cycle: Number of initial short cycles.
            initial_length: Number of iterations in initial short cycles.
            alter_cell_cons: Alternate cell constraints on/off between cycles.
        """
        self.exe = exe
        self.struct_name = struct_name
        self.maxit = maxit
        self.initial_length = initial_length
        if maxit == 0:
            self._init_relax = 0
            self.success = 2
        else:
            self._init_relax = initial_cycle
            self.success = 3

        self.alter_cell_cons = alter_cell_cons

    @property
    def dot_castep(self) -> str:
        return self.struct_name + ".castep"

    @property
    def dot_param(self) -> str:
        return self.struct_name + ".param"

    @property
    def dot_cell(self) -> str:
        return self.struct_name + ".cell"

    @property
    def dot_cell_out(self) -> str:
        return self.struct_name + "-out.cell"

    def _run_castep(self, timeout: Optional[float] = None) -> None:
        """Run the CASTEP executable."""
        args = " ".join([self.exe, self.struct_name])
        logger.debug("Calling subprocess '%s'", args)
        subprocess.call(args, shell=True, timeout=timeout)
        logger.debug("Checking %s for results", self.dot_castep)
        with open(self.dot_castep) as fh:
            line_container = []
            for line in fh:
                line_container.append(line)
                if len(line_container) > 20:
                    line_container.pop(0)
            if "Total time" not in "".join(line_container):
                raise RelaxError(
                    f"CASTEP finished with error. struct_name: {self.struct_name}"
                )

    def _set_short_geom_iter(self) -> None:
        """Set the geom iteration for initial rough relaxation."""
        with fileinput.input(self.dot_param, inplace=True) as f:
            for line in f:
                if "geom_max_iter" in line.lower():
                    print("#" + line, end="")
                else:
                    print(line, end="")
        with open(self.dot_param, "a") as f:
            f.write(f"\ngeom_max_iter: {self.initial_length}#MARKED\n")

    def _recover_geom_iter(self) -> None:
        """Recover the original geom_max_iter from the PARAM file."""
        with fileinput.input(self.dot_param, inplace=True) as f:
            for line in f:
                if line.startswith("#") and "geom_max_iter" in line.lower():
                    print(line[1:], end="")
                elif "#MARKED" in line:
                    continue
                else:
                    print(line, end="")

    def check_param(self) -> None:
        """Ensure write_cell_structure is set in the PARAM file."""
        cell_out = False
        with open(self.dot_param, "r+") as f:
            for line in f:
                if "write_cell_structure" in line.lower():
                    cell_out = True
            if not cell_out:
                f.write("write_cell_structure :  True")

        with fileinput.input(self.dot_param, inplace=True) as f:
            for line in f:
                if "#MARKED" in line:
                    pass
                elif line.startswith("#") and "geom_max_iter" in line.lower():
                    print(line[1:], end="")
                else:
                    print(line, end="")

    def save_state(self) -> None:
        """Save the relaxation state to a JSON file."""
        s_dict = {"_init_relax": self._init_relax, "success": self.success}
        sfile = "." + self.struct_name + "-fr.json"
        logger.debug("Save state in file: %s", sfile)
        with open(sfile, "w") as fh:
            json.dump(s_dict, fh)

    def load_state(self) -> None:
        """Load state from an unfinished run if present."""
        state_file = "." + self.struct_name + "-fr.json"
        if os.path.isfile(state_file):
            logger.debug("Recovering from %s", state_file)
            with open(state_file) as fh:
                s_dict = json.load(fh)
            for key, value in s_dict.items():
                setattr(self, key, value)
        self._cell_from_geom()

    def fixed_cell_off(self) -> None:
        """Remove fixed cell constraints and uncomment CELL_CONSTRAINTS block."""
        with open(self.struct_name + ".cell") as f:
            lines = f.readlines()

        in_block = False
        new_lines: list[str] = []
        for line in lines:
            if "FIX_ALL_CELL" in line.upper():
                continue
            if "%BLOCK CELL_CONSTRAINTS" in line.upper():
                in_block = True
            elif "%ENDBLOCK CELL_CONSTRAINTS" in line.upper():
                in_block = False
                if line[0] == "#":
                    new_lines.append(line[1:])
                else:
                    new_lines.append(line)
                continue
            if in_block:
                if line[0] == "#":
                    new_lines.append(line[1:])
                else:
                    new_lines.append(line)
            else:
                new_lines.append(line)

        with open(self.struct_name + ".cell", "w") as f:
            f.write("".join(new_lines))

    def fixed_cell_on(self) -> None:
        """Comment out CELL_CONSTRAINTS block and add FIX_ALL_CELL."""
        with open(self.struct_name + ".cell") as f:
            lines = f.readlines()

        in_block = False
        new_lines: list[str] = []
        for line in lines:
            if "%BLOCK CELL_CONSTRAINTS" in line.upper():
                in_block = True
            elif "%ENDBLOCK CELL_CONSTRAINTS" in line.upper():
                in_block = False
                new_lines.append("#" + line)
                continue
            if in_block:
                new_lines.append("#" + line)
            else:
                new_lines.append(line)

        new_lines.append("\nFIX_ALL_CELL: TRUE\n")

        with open(self.struct_name + ".cell", "w") as f:
            f.write("".join(new_lines))

    def _cell_from_geom(self) -> None:
        """Push cell and positions from the last .geom file into the .cell file."""
        struct_name = self.struct_name
        gfile = struct_name + ".geom"
        if not os.path.isfile(gfile):
            logger.debug("No geom file found for %s", struct_name)
            return

        logger.debug("Using structure from %s.geom", struct_name)
        with open(struct_name + ".cell") as cin:
            old_rest = filter_out_stream(cin, r"^%BLOCK lat", r"^%ENDBLOCK lat")
            old_rest = filter_out_stream(old_rest, r"^%BLOCK pos", r"^%ENDBLOCK pos")
            old_rest_content = old_rest.read()
        with open(struct_name + ".cell.tmp", "w") as f:
            cblock, pblock = geom_to_cell(gfile)
            f.write(cblock)
            f.write(pblock)
            f.write(old_rest_content)

        shutil.move(struct_name + ".cell.tmp", struct_name + ".cell")

    def run(self, timeout: Optional[float] = None) -> str:
        """
        Run the full relaxation procedure.

        Args:
            timeout: Walltime limit in seconds.

        Returns:
            ``"success"`` if converged, ``"limit reached"`` if maxit exceeded.
        """
        import time

        self.load_state()
        self.check_param()

        start = time.time()

        def remaining_timeout():
            return timeout - (time.time() - start)

        count = 1
        while self._init_relax > 0:
            logger.debug("Starting initial iteration %d", count)
            self._set_short_geom_iter()
            self._run_castep(remaining_timeout())
            self._recover_geom_iter()
            push_cell(self.dot_cell_out, self.dot_cell)
            self._init_relax -= 1
            self.save_state()
            logger.debug("Finished initial iteration %d", count)
            count += 1

        count = 1
        while self.success > 1:
            logger.debug("Starting iteration %d", count)

            if self.alter_cell_cons:
                if count % 2 == 1:
                    logger.debug("Remove fixed cell in %s", self.struct_name + ".cell")
                    self.fixed_cell_off()
                else:
                    logger.debug(
                        "Turn partial cell constraints back on in %s",
                        self.struct_name + ".cell",
                    )
                    self.fixed_cell_on()

            self._run_castep(remaining_timeout())
            push_cell(self.dot_cell_out, self.dot_cell)

            r_flag, i_count = check_relax_status(self.dot_castep)
            if r_flag:
                if self.success > 1:
                    self.success -= 1
            elif i_count > self.maxit:
                self.success = -1
            else:
                self.success = 3

            logger.debug("Finished iteration %d", count)
            logger.debug("success: %d", self.success)
            self.save_state()
            count += 1

        state_file = "." + self.struct_name + "-fr.json"
        if os.path.isfile(state_file):
            os.remove(state_file)
        if self.success == 1:
            return "success"
        return "limit reached"

    def __repr__(self) -> str:
        return f"<FullRelax(exe={self.exe!r}, struct_name={self.struct_name!r})>"


def check_relax_status(dot_castep: str) -> tuple[bool, int]:
    """
    Check whether the last geometry optimisation completed.

    Args:
        dot_castep: Path to the .castep file.

    Returns:
        Tuple of (relaxation_succeeded, total_iterations).
    """
    geom_line = re.compile(r"Geometry optimization (\w+)")
    count = 0
    g_res = ""

    with open(dot_castep) as fh:
        for line in fh:
            m = geom_line.search(line)
            if m:
                g_res = m.group(1)
            if ": finished iteration" in line:
                count += 1

    return g_res == "completed", count


def parse_geom_text_output(out_lines: list[str]) -> dict[str, np.ndarray]:
    """
    Parse a CASTEP .geom file and return trajectory data.

    Args:
        out_lines: Lines from the .geom file.

    Returns:
        Dictionary with keys ``cells``, ``positions``, ``forces``,
        ``geom_energy``, ``symbols``, and optionally ``velocities``.
        All arrays are in CASTEP atomic units, converted to
        Angstrom/eV.
    """
    hartree = _units["Eh"]
    bohr = _units["a0"]

    cell_list: list = []
    species_list: list = []
    geom_list: list = []
    forces_list: list = []
    energy_list: list = []
    velocity_list: list = []

    current_pos: list = []
    current_species: list = []
    current_forces: list = []
    current_velocity: list = []
    current_cell: list = []
    in_header = False

    for line in out_lines:
        if "begin header" in line.lower():
            in_header = True
            continue
        if "end header" in line.lower():
            in_header = False
            continue
        if in_header:
            continue

        sline = line.split()
        if "<-- E" in line:
            energy_list.append(float(sline[0]) * hartree)
        elif "<-- h" in line:
            current_cell.append(list(map(float, sline[:3])))
        elif "<-- R" in line:
            current_pos.append(list(map(float, sline[2:5])))
            current_species.append(sline[0])
        elif "<-- F" in line:
            current_forces.append(list(map(float, sline[2:5])))
        elif "<-- V" in line:
            current_velocity.append(list(map(float, sline[2:5])))
        elif "<-- T" in line:
            pass  # temperature not collected
        elif not line.strip() and current_cell:
            cell_list.append(current_cell)
            species_list.append(current_species)
            geom_list.append(current_pos)
            forces_list.append(current_forces)
            current_cell = []
            current_species = []
            current_pos = []
            current_forces = []
            if current_velocity:
                velocity_list.append(current_velocity)
                current_velocity = []

    if not species_list:
        raise RuntimeError("No data found in geom file")

    out: dict[str, object] = {
        "cells": np.array(cell_list) * bohr,
        "positions": np.array(geom_list) * bohr,
        "forces": np.array(forces_list) * hartree / bohr,
        "geom_energy": np.array(energy_list),
        "symbols": species_list[0],
    }
    if velocity_list:
        out["velocities"] = np.array(velocity_list) * bohr
    return out


def geom_to_cell(geom_file: str) -> tuple[str, str]:
    """
    Convert the last configuration in a .geom file to cell blocks.

    Args:
        geom_file: Path to the .geom file.

    Returns:
        Tuple of (cell_block_string, position_block_string) including
        ``%BLOCK``/``%ENDBLOCK`` enclosures.
    """
    with open(geom_file) as fh:
        glines = fh.read().split("\n")

    g_info = parse_geom_text_output(glines)

    last_cell = g_info["cells"][-1]
    last_pos = g_info["positions"][-1]

    pos_block = ["%BLOCK POSITIONS_ABS"]
    for specie, pos in zip(g_info["symbols"], last_pos):
        pos_block.append(f"{specie} {pos[0]:.7f} {pos[1]:.7f} {pos[2]:.7f}")
    pos_block.append("%ENDBLOCK POSITIONS_ABS\n")

    cell_block = ["%BLOCK LATTICE_CART"]
    for v in last_cell:
        cell_block.append(f"{v[0]:.7f} {v[1]:.7f} {v[2]:.7f}")
    cell_block.append("%ENDBLOCK LATTICE_CART\n")

    return "\n".join(cell_block), "\n".join(pos_block)

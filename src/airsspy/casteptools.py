"""
CASTEP-related utility functions.

Provides tools for parsing CASTEP output files, extracting results,
managing cell files, and running RASH (Random Ab initio Structure
Hunting) relaxation seeds.
"""

import io
import os
import platform
import re
import shutil
import time
import uuid

import ase.io
import numpy as np

from .utils import filter_out_stream, trim_stream


class CastepRunError(RuntimeError):
    """Error when CASTEP produced no result"""


class CastepSkip(RuntimeError):
    """Skip this CASTEP run"""


class CastepManualTimedout(RuntimeError):
    """CASTEP was manually timed out"""


def parse_param(seed: str) -> dict:
    """
    Parse a CASTEP .param file and return a dictionary.

    Keys are converted to lower case. Values are left as-is.

    Args:
        seed: Name of the seed (without extension).

    Returns:
        Dictionary of parameter key-value pairs.
    """
    dict_out = {}
    spliter = re.compile(r"[ :=]+")
    with open(seed + ".param") as seedfile:
        for line in seedfile:
            if "#" in line:
                continue
            pair = line.strip()
            pair = list(filter(None, spliter.split(pair, 1)))
            if pair:
                dict_out.update({pair[0].lower(): pair[1]})
    return dict_out


pattern_geom = r"""
^\ *
(?P<name>[A-Z]+):          # Optimisation name
\ +finished\ iteration
\ +(?P<number>[0-9]+)             # iteration number
\ +with\ enthalpy=
\ +(?P<H>[+-.E0-9]+)        # Enthalpy
\ +(?P<unit>[a-zA-Z]+)          # Enthalpy unit
\ *$
"""

pattern_time = r"""
^\ +
[0-9]+                      # Loop number
\ +
[+-.E0-9]+                   # Energy
\ +
[+-.E0-9]+                   # Fermi Energy
\ +
[+-.E0-9]+                   # Energy gain per atom
\ +
([.0-9]+)                      # Timer
\ +
<--\ SCF
$
"""


def parse_dot_castep(fb, aggregate=False):  # pylint: disable=too-many-locals
    """
    Parse geometry convergence information from a CASTEP output file.

    Args:
        fb: File handle or iterable of lines from a .castep file.
        aggregate: If True, aggregate timer values across continuation runs.

    Returns:
        Dictionary with keys: H (enthalpy list), iter_num, name, unit,
        time (timing), save_iter, save_times.
    """
    iter_start = re.compile(r" Starting \w+ iteration")
    sfc_line = re.compile(pattern_time, re.VERBOSE)
    match_geom = re.compile(pattern_geom, re.VERBOSE)

    eng = []
    iter_num = []
    iter_times = []
    save_times = []
    save_iter = []
    geom_name, unit = None, None

    last_time = 0
    current_iter = 0

    for line in fb:
        if iter_start.match(line):
            current_iter += 1
            continue

        if "Writing model" in line:
            save_times.append(last_time)
            save_iter.append(current_iter)

        scf_match = sfc_line.match(line)
        if scf_match:
            last_time = float(scf_match.group(1))
            continue

        geom_match = match_geom.match(line)
        if geom_match:
            eng.append(float(geom_match.group("H")))
            iter_num.append(int(geom_match.group("number")))
            iter_times.append(float(last_time))
            geom_name = geom_match.group("name")
            unit = geom_match.group("unit")
            continue

    out = {
        "H": eng,
        "iter_num": iter_num,
        "name": geom_name,
        "unit": unit,
        "time": iter_times,
        "save_iter": save_iter,
        "save_times": save_times,
    }
    if aggregate:
        timer_array = np.array(iter_times)
        last_record = 0
        for i, time_record in enumerate(iter_times):
            if time_record < last_record:
                timer_array[i:] = timer_array[i:] + last_record
            last_record = time_record
        out.update(time=timer_array)
    return out


def RASH_prepare_seed(seed: str, relaxed: str, amp: float) -> io.StringIO:
    """
    Prepare a seed file for RASH (Random Ab initio Structure Hunting) relaxation.

    Takes the lattice from the relaxed structure and the position block,
    then combines with the original seed's generation parameters.

    Args:
        seed: Name of the ORIGINAL seed.
        relaxed: Name (uid) of the relaxed structure.
        amp: Amplitude of shape perturbation (positive float).

    Returns:
        An StringIO object of the cell for RASH.
    """
    with open(relaxed + ".cell") as rlx_cell, open(seed + ".cell") as seed_cell:
        rlx_lattice = trim_stream(rlx_cell, r"^%BLOCK [Ll][Aa][Tt]", r"^%ENDBLOCK [Ll][Aa][Tt]", ["FIX_VOL", "ANG"]).read()
        rlx_lattice = rlx_lattice.replace("%ENDBLOCK", "#FIX\nENDBLOCK")
        rlx_pos = trim_stream(rlx_cell, r"^%BLOCK [Pp][Oo][Ss]", r"^%ENDBLOCK [Pp][Oo][Ss]").read()
        seed_rest = filter_out_stream(seed_cell, r"^%BLOCK [Ll][Aa][Tt]", r"^%ENDBLOCK [Ll][Aa][Tt]")
        seed_rest = filter_out_stream(
            seed_rest, r"^%BLOCK [Pp][Oo][Ss]", r"^%ENDBLOCK [Pp][Oo][Ss]", ["#POSAMP=", "#SYMMOPS=", "#NFORM=", "#SUPER"]
        ).read()
    out_string = io.StringIO()
    out_string.write(rlx_lattice + "\n")
    out_string.write(rlx_pos + "\n")
    out_string.write(seed_rest + "\n" + "#POSAMP=" + str(amp))
    out_string.seek(0)
    return out_string


def extract_REM(seed: str) -> dict:
    """
    Construct the REM lines from a CASTEP history file.

    Args:
        seed: Name of the seed (without extension).

    Returns:
        A dictionary of the REM information.
    """
    rem = {}
    with open(seed + ".history") as his:
        get_pspot = False
        pspot = []
        for line in his:
            if " Welcome to " in line:
                rem.update(version=line.strip(" |\n").split()[-1])
            if " from code version " in line:
                rem.update(fromcode=line.strip(" |\n"))
            if " Run started: " in line:
                rem.update(rundate=line.strip(" |\n"))
            if " using functional " in line:
                xc = line.split(":")[1].strip()
                rem.update(functional="Functional  " + xc)
            if " relativistic treatment " in line:
                rem.update(relativity=" Relativity " + line.split(":")[-1].strip())
            if " plane wave basis set cut-off " in line:
                rem.update(cutoff="Cut-off " + line.split(":")[-1].strip())
            if " size of standard grid " in line:
                rem.update(gridscale=" Grid scale " + line.split(":")[-1].strip())
            if " size of   fine   gmax " in line:
                rem.update(gmax=" Gmax " + line.split(":")[-1].strip())
            if " finite basis set correction  " in line:
                rem.update(fbsc=" FBSC" + line.split(":")[-1].strip())
            if " MP grid size for SCF calculation is " in line:
                rem.update(mpgrid="MP grid " + " ".join(line.strip().split()[-3:]))
            if " with an offset of  " in line:
                rem.update(offset=" Offset " + " ".join(line.strip().split()[-3:]))
            if " Number of kpoints used = " in line:
                rem.update(nkpts=" No. kpts " + line.split("=")[-1].strip())
            if not pspot and " Files used for pseudopotential" in line:
                get_pspot = True
            if get_pspot and "-------------------------------" in line:
                get_pspot = False
            if get_pspot:
                pspot.append("REM " + line.strip() + "\n")
    pspot_str = "".join(pspot).replace(":", " ")
    rem.update(pspot=pspot_str)

    with open(seed + ".cell") as cell:
        for line in cell:
            if "KPOINTS_MP_SPACING" in line:
                rem.update(spacing=" Spacing " + re.split(r"[:=\s]+", line.strip())[-1])
                break
        if "spacing" not in rem:
            rem.update(spacing="Spacing Unspecified")

    rem.update(projectdir=os.getcwd())
    rem.update(host=platform.node())
    return rem


def extract_result(seed: str) -> dict:
    """
    Extract and save results as attributes from a CASTEP history file.

    Args:
        seed: Name of the seed (without extension).

    Returns:
        Dictionary with P, H, V, sym, nat, chem_formula.
    """
    results = {}
    with open(seed + ".history") as res_file:
        for line in res_file:
            if "Pressure:" in line:
                results["pl"] = line
            if "Final free" in line or "Final Enthalpy" in line:
                results["hl"] = line
            if "Current cell volume =" in line:
                results["cl"] = line

    pl = results.get("pl")
    if pl is not None:
        P_str = pl.split(":")[1].strip("|*\n ")
    else:
        P_str = "0"

    if "hl" not in results or "cl" not in results:
        raise RuntimeError(f"Energy or cell volume not found in {seed}.history")
    H_str = results["hl"].split("=")[1].split()[0].strip()
    V_str = results["cl"].split("=")[1].split()[0].strip()

    P = float(P_str)
    H = float(H_str)
    V = float(V_str)

    atoms = ase.io.read(seed + ".cell")
    import spglib

    sg = spglib.get_spacegroup(atoms, 0.1)
    sg = sg.split()[0]
    sg = "(" + sg + ")"
    nat = atoms.get_number_of_atoms()
    chem_formula = atoms.get_chemical_formula()
    return {"P": P, "H": H, "V": V, "sym": sg, "nat": nat, "chem_formula": chem_formula}


def write_converge(seed: str, suffix: str = "castep") -> None:
    """
    Write convergence information to a .gconv file.

    Args:
        seed: Name of the seed (without extension).
        suffix: File suffix (default: ``"castep"``).
    """
    with open(seed + "." + suffix) as fb:
        conv_info = parse_dot_castep(fb, aggregate=True)
    with open(seed + ".gconv", "w") as out_file:
        out_file.write("# Geometry Optimisation convergence\n")
        out_file.write("# Method:  " + conv_info["name"] + "\n")
        out_file.write("# Energy Unit:  " + conv_info["unit"] + "\n")
        num = range(len(conv_info["H"]))
        for i, value, iter_number, timer in zip(num, conv_info["H"], conv_info["iter_num"], conv_info["time"]):
            out_file.write(f"{i:<5d}{value:<20.5f}{iter_number:<10d}{timer:.2f}\n")


def get_rand_cell_name(seed_name: str) -> str:
    """
    Return a string for naming a randomly generated cell file.

    Format: ``<seed_name>-<YYMMDD-HHMMSS>-<uuid6>.cell``

    Args:
        seed_name: Name of the seed.

    Returns:
        A unique cell file name.
    """
    timestamp = time.strftime("%y%m%d-%H%M%S")
    return seed_name + "-" + timestamp + "-" + str(uuid.uuid4())[:6] + ".cell"


def castep_finish_ok(dot_castep: str) -> bool:
    """
    Check for the hint that CASTEP finished OK.

    Args:
        dot_castep: Path to the .castep file.

    Returns:
        True if ``"Total time"`` appears in the last 20 lines.
    """
    if not os.path.isfile(dot_castep):
        return False

    with open(dot_castep) as fhandle:
        lines = fhandle.readlines()

    last_few = lines[-20:]
    for line in last_few:
        if "Total time" in line:
            return True
    return False


def castep_geom_count(dot_castep: str) -> int:
    """
    Count the number of geometry optimisation cycles in a .castep file.

    Args:
        dot_castep: Path to the .castep file.

    Returns:
        Number of geometry cycles.
    """
    count = 0
    with open(dot_castep) as fhandle:
        for line in fhandle:
            if "starting iteration" in line:
                count += 1
    return count


def push_cell(cellout: str, cell: str) -> None:
    """
    Copy the structure from one cell file to another.

    Moves the ``LATTICE`` and ``POSITIONS`` blocks from ``cellout`` into
    ``cell``, preserving everything else from ``cell``.

    Args:
        cellout: Path to the output cell file (source of new structure).
        cell: Path to the input cell file (to be updated).
    """
    with open(cellout) as fhout, open(cell) as fcell, open(str(cell) + ".tmp", "w") as ftmp:
        new_lattice = trim_stream(fhout, r"^%BLOCK [Ll][Aa][Tt]", r"^%ENDBLOCK [Ll][Aa][Tt]", ["FIX_VOL", "ANG"])
        new_pos = trim_stream(fhout, r"^%BLOCK [Pp][Oo][Ss]", r"^%ENDBLOCK [Pp][Oo][Ss]")
        ftmp.write(new_lattice.read())
        ftmp.write("\n")
        ftmp.write(new_pos.read())
        old_rest = filter_out_stream(fcell, r"^%BLOCK [Ll][Aa][Tt]", r"^%ENDBLOCK [Ll][Aa][Tt]")
        old_rest = filter_out_stream(old_rest, r"^%BLOCK [Pp][Oo][Ss]", r"^%ENDBLOCK [Pp][Oo][Ss]")
        ftmp.write(old_rest.read())

    shutil.move(str(cell) + ".tmp", cell)
    os.remove(cellout)


def gulp_relax_finish_ok(dot_castep: str) -> bool:
    """
    Check if a GULP relaxation finished successfully.

    Args:
        dot_castep: Path to the output file.

    Returns:
        True if ``"Final Enthalpy"`` is found.
    """
    if not os.path.isfile(dot_castep):
        return False

    with open(dot_castep) as fh:
        content = fh.read()

    match = re.search(r"Final Enthalpy += ([-e0-9\.]+)", content)
    return bool(match)

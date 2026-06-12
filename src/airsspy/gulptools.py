"""
GULP calculation utilities.

Provides tools for monitoring GULP geometry optimisation progress,
checking calculation health, and running GULP with automatic
termination on divergence.
"""

import re
import subprocess
import sys
from collections import namedtuple
from collections.abc import Iterable
from typing import Union

Ginfo = namedtuple("Ginfo", ["cycle", "energy", "gnorm", "cpu"])

GEOM_LINE_PATTERN = re.compile(
    r"^ +Cycle:([ 0-9]+)Energy:([-0-9e*. ]+)Gnorm:([-0-9e*. ]+)CPU:([-0-9e. ]+)"
)


def geom_opt_progress(
    gfile: Union[str, Iterable[str]],
) -> list[Ginfo]:
    """
    Parse geometry optimisation info from a GULP output file.

    Args:
        gfile: Path to the GULP output file, or an iterable of lines.

    Returns:
        A list of Ginfo namedtuples (cycle, energy, gnorm, cpu).
    """
    if isinstance(gfile, str):
        with open(gfile) as fhandle:
            fcontent = fhandle.readlines()
    else:
        fcontent = list(gfile)

    data = []
    for line in fcontent:
        match = GEOM_LINE_PATTERN.match(line)
        if match:
            data.append(Ginfo(*[match.group(i).strip() for i in range(1, 5)]))

    return data


def check_gulp(gfile: Union[str, Iterable[str]]) -> bool:
    """
    Check if a GULP calculation is progressing well.

    Returns False if energy or Gnorm fields contain overflow markers
    (``*``), or if the Gnorm diverges beyond 10 and keeps increasing.

    Args:
        gfile: Path to the GULP output file, or an iterable of lines.

    Returns:
        True if the calculation appears healthy.
    """
    gopt = geom_opt_progress(gfile)
    gnorms = []
    for entry in gopt:
        if "*" in entry.energy or "*" in entry.gnorm:
            return False
        gnorms.append(float(entry.gnorm))

        if len(gnorms) > 5 and gnorms[-1] > 10:
            if gnorms[-1] > gnorms[-4]:
                return False

    return True


def guarded_gulp(exe: str = "gulp") -> bool:
    """
    Run GULP with automatic termination on divergence.

    Monitors GULP stdout in real time and terminates the process if
    the energy or Gnorm overflows.

    Args:
        exe: Name or path of the GULP executable.

    Returns:
        True if GULP completed successfully, False if terminated.
    """
    proc = subprocess.Popen(
        [exe], stdin=sys.stdin, stdout=subprocess.PIPE, universal_newlines=True
    )

    run_ok = True
    data_lines: list[str] = []
    while proc.poll() is None:
        new_line = proc.stdout.readline()
        print(new_line, end="")
        data_lines.append(new_line)
        run_ok = check_gulp(data_lines)
        if not run_ok:
            proc.terminate()
    print(proc.stdout.read(), end="")
    sys.stdout.flush()

    return run_ok

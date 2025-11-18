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
General utility functions for AIRSS workflows
"""
import io
import re
from typing import List, Optional, TextIO, Tuple

import numpy as np
from ase import Atoms


def trim_stream(
    stream: TextIO, start: str, end: str, extra_remove: Optional[List[str]] = None
) -> io.StringIO:
    """
    Select a portion of a stream, return the portion without lines
    matched with extra_remove keywords.

    Args:
        stream: Input stream to trim
        start: Regular expression pattern to start recording
        end: Regular expression pattern to stop recording
        extra_remove: List of additional patterns to filter out

    Returns:
        StringIO object containing the trimmed content
    """
    out_stream = io.StringIO()
    record = False
    stream.seek(0)
    if extra_remove is None:
        extra_remove = []

    for line in stream:
        if re.match(start, line, re.IGNORECASE):
            record = True

        if record is True:
            # Apply filter, if there is the keyword then skip
            appear = False
            for pattern in extra_remove:
                if pattern in line:
                    appear = True
                    break
            # Only write when there is no match
            if not appear:
                out_stream.write(line)

        # Turn recording off at the end line
        if re.match(end, line, re.IGNORECASE):
            break

    # Reset streams
    stream.seek(0)
    out_stream.seek(0)
    return out_stream


def filter_out_stream(stream: TextIO, start: str, end: str) -> io.StringIO:
    """
    Opposite of trim_stream, only the portion outside start end is selected.

    Args:
        stream: Input stream to filter
        start: Regular expression pattern to start filtering out
        end: Regular expression pattern to stop filtering out

    Returns:
        StringIO object containing the filtered content
    """
    out_stream = io.StringIO()
    record = True
    stream.seek(0)

    for line in stream:
        if re.match(start, line, re.IGNORECASE):
            record = False

        if record is True:
            out_stream.write(line)

        # Turn recording back on at the end line
        if re.match(end, line, re.IGNORECASE):
            record = True

    # Reset streams
    stream.seek(0)
    out_stream.seek(0)
    return out_stream


def calc_kpt_tuple_recip(
    structure: Atoms, mp_spacing: float = 0.05, rounding: str = "up"
) -> Tuple[int, int, int]:
    """
    Calculate reciprocal-space sampling with real-space parameter.

    This function calculates k-point mesh based on a real-space spacing
    parameter, which is more intuitive than directly specifying k-points.

    Args:
        structure: ASE Atoms object containing the structure
        mp_spacing: Real-space spacing parameter in Å^-1
        rounding: Rounding method - "up" for ceiling, "down" for floor with +0.5

    Returns:
        Tuple of k-point mesh (kx, ky, kz)

    Note:
        Uses pymatgen for reciprocal lattice calculation.
    """
    # Import pymatgen for reciprocal lattice calculation
    from pymatgen.io.ase import AseAtomsAdaptor

    pmg_structure = AseAtomsAdaptor.get_structure(structure)
    # Get reciprocal lattice vectors with pymatgen. Note that pymatgen does include
    # the 2*pi factor used in many definitions of these vectors; hence it is divided by 2π for the
    # CASTEP convention
    recip_cell = pmg_structure.lattice.reciprocal_lattice.matrix / np.pi / 2

    # Get reciprocal cell vector magnitudes according to Pythagoras' theorem
    abc_recip = np.sqrt(np.sum(np.square(recip_cell), 1))
    k_samples = abc_recip / mp_spacing

    # Rounding
    if rounding == "up":
        k_samples = np.ceil(k_samples)
    else:
        k_samples = np.floor(k_samples + 0.5)

    return tuple(int(x) for x in k_samples)


def unique(items: List) -> List:
    """
    Get a list of ordered unique items.

    Args:
        items: List of items to deduplicate

    Returns:
        List containing unique items in original order
    """
    out = []
    for item in items:
        if item not in out:
            out.append(item)
    return out


def stream_to_list(stream: TextIO) -> List[str]:
    """
    Convert a stream to a list of strings.

    Args:
        stream: Input stream

    Returns:
        List of lines from the stream
    """
    stream.seek(0)
    return [line.rstrip() for line in stream]


def safe_cast_float(value: str, default: float = 0.0) -> float:
    """
    Safely convert a string to float with fallback value.

    Args:
        value: String to convert
        default: Default value if conversion fails

    Returns:
        Float value or default
    """
    try:
        return float(value)
    except (ValueError, TypeError):
        return default


def safe_cast_int(value: str, default: int = 0) -> int:
    """
    Safely convert a string to int with fallback value.

    Args:
        value: String to convert
        default: Default value if conversion fails

    Returns:
        Integer value or default
    """
    try:
        return int(value)
    except (ValueError, TypeError):
        return default


def extract_number_from_string(text: str, default: float = 0.0) -> float:
    """
    Extract the first number from a string.

    Args:
        text: String containing numbers
        default: Default value if no number found

    Returns:
        First number found or default
    """
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)
    if match:
        return float(match.group())
    return default


def find_pattern_in_file(
    filename: str, pattern: str, max_matches: int = 0
) -> List[str]:
    """
    Find all lines matching a pattern in a file.

    Args:
        filename: Path to file to search
        pattern: Regular expression pattern to match
        max_matches: Maximum number of matches to return (0 for unlimited)

    Returns:
        List of matching lines
    """
    matches = []
    compiled_pattern = re.compile(pattern)

    try:
        with open(filename) as f:
            for line in f:
                if compiled_pattern.search(line):
                    matches.append(line.rstrip())
                    if max_matches > 0 and len(matches) >= max_matches:
                        break
    except FileNotFoundError:
        pass

    return matches


def count_pattern_in_file(filename: str, pattern: str) -> int:
    """
    Count occurrences of a pattern in a file.

    Args:
        filename: Path to file to search
        pattern: Regular expression pattern to match

    Returns:
        Number of matches found
    """
    compiled_pattern = re.compile(pattern)
    count = 0

    try:
        with open(filename) as f:
            for line in f:
                if compiled_pattern.search(line):
                    count += 1
    except FileNotFoundError:
        pass

    return count


def format_time_elapsed(seconds: float) -> str:
    """
    Format elapsed time in human-readable format.

    Args:
        seconds: Time in seconds

    Returns:
        Formatted time string
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.0f}s"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        secs = seconds % 60
        return f"{hours}h {minutes}m {secs:.0f}s"

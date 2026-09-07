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
"""Compatibility exports for commonly used airsspy objects."""

from importlib import import_module
from typing import Final

__version__: Final = "0.1.4"

_LAZY_EXPORTS = {
    "Buildcell": ".build",
    "SeedAtoms": ".seed",
    "RESFile": ".restools",
    "TitlInfo": ".restools",
    "extract_res": ".restools",
    "format_minsep": ".restools",
    "get_minsep": ".restools",
    "get_spacegroup_atoms": ".restools",
    "parse_titl": ".restools",
    "read_res_atoms": ".restools",
    "read_res_pmg": ".restools",
    "save_airss_res": ".restools",
    "calc_kpt_tuple_recip": ".utils",
    "count_pattern_in_file": ".utils",
    "extract_number_from_string": ".utils",
    "filter_out_stream": ".utils",
    "find_pattern_in_file": ".utils",
    "format_time_elapsed": ".utils",
    "safe_cast_float": ".utils",
    "safe_cast_int": ".utils",
    "stream_to_list": ".utils",
    "trim_stream": ".utils",
    "unique": ".utils",
}


def __getattr__(name: str):
    if name not in _LAZY_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(_LAZY_EXPORTS[name], __name__)
    value = getattr(module, name)
    globals()[name] = value
    return value


__all__ = [
    # Core classes
    "SeedAtoms",
    "Buildcell",
    # RES file handling
    "RESFile",
    "TitlInfo",
    "extract_res",
    "save_airss_res",
    "parse_titl",
    "read_res_atoms",
    "read_res_pmg",
    "get_spacegroup_atoms",
    "get_minsep",
    "format_minsep",
    # Utilities
    "unique",
    "trim_stream",
    "filter_out_stream",
    "calc_kpt_tuple_recip",
    "stream_to_list",
    "safe_cast_float",
    "safe_cast_int",
    "extract_number_from_string",
    "find_pattern_in_file",
    "count_pattern_in_file",
    "format_time_elapsed",
]
